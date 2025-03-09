library(tidyverse)
library(splines)


simulate_genotypes <- function(N, M, min_maf, max_maf, 
                   seed = 1) {
  set.seed(seed)
  
  mafs <- runif(M, min_maf, max_maf)
  
  g_mat <- sapply(mafs, function(maf) {
    rbinom(N, 2, maf)
  }, simplify = TRUE)
  colnames(g_mat) <- paste0("g", seq(1, M))
  
  saveRDS(g_mat, "../data/processed/simulations/g_mat.rds")
  saveRDS(mafs, "../data/processed/simulations/maf_vec.rds")
}


simulate_exposure <- function(N, e_distr, e_icc, ge_var_tot, g_mat, 
                              sim_idx) {
  if (e_distr == "normal") {
    beta_ge_vec <- rep(sqrt(ge_var_tot / ncol(g_mat)), ncol(g_mat))
    e_means <- as.vector(g_mat %*% beta_ge_vec)
    var_e <- 1 - ge_var_tot
    beta_e_em <- sqrt(e_icc)  # Such that true E explains amount of variance equal to ICC
    var_em <- 1 - e_icc  # Such that measured E has total variance equal to 1
    e_df <- tibble(
      e = rnorm(N, e_means, sqrt(var_e)),
      e_m = rnorm(N, beta_e_em * e, sqrt(var_em))
    )
  } else if (e_distr == "normal10") {
    e_df <- tibble(
      e = rnorm(N, 10, 1),
      e_m = e
    )
  } else if (e_distr == "gamma") {
    beta_ge_vec <- rep(sqrt(ge_var_tot / ncol(g_mat)), ncol(g_mat))
    gamma_mult <- exp(as.vector(g_mat %*% beta_ge_vec))
    e_df <- tibble(
      e = rgamma(N, shape = 1, rate = 1) * gamma_mult
    ) %>%
      mutate(e = e / sd(e),
             e_m = e)
  }
  names(e_df) <- paste0(names(e_df), sim_idx)
  e_df
}


simulate_scenario <- function(ge_var_tot, e_distr, e_icc,
                              e_var_tot, g_var_tot, gxe_var_tot, nl_e,
                              n_sim, tag,
                              g_mat, maf_vec,
                              train_prop) {
  print(paste0("Simulating ", tag, "..."))
  
  N <- nrow(g_mat)
  M <- ncol(g_mat)
  
  sim_df <- tibble(
    id = seq(1, N),
    train = sample(c(0, 1), N, replace = TRUE, 
                   prob = c(1 - train_prop, train_prop))
  )
  
  e_df <- map(seq(1, n_sim), function(idx) {
    simulate_exposure(N, e_distr, e_icc, ge_var_tot, g_mat, idx)
  }) %>%
    bind_cols()
  
  beta_e <- sqrt(e_var_tot)
  
  g_mat_variances <- 2 * maf_vec * (1 - maf_vec)
  
  g_var_single <- g_var_tot / M
  beta_g_vars <- g_var_single / g_mat_variances
  beta_g_vec <- rnorm(M, 0, sqrt(beta_g_vars))
  # beta_g_vec <- sqrt(beta_g_vars)
  # beta_g_vec <- sample(c(-1, 1), M, replace = TRUE) * sqrt(beta_g_vars)

  gxe_var_single <- gxe_var_tot / M
  beta_gxe_vars <- gxe_var_single / g_mat_variances
  beta_gxe_vec <- rnorm(M, 0, sqrt(beta_gxe_vars))
  # beta_gxe_vec <- sqrt(beta_gxe_vars)
  # beta_gxe_vec <- sample(c(-1, 1), M, replace = TRUE) * sqrt(beta_gxe_vars)
  
  error_var <- 1 - e_var_tot - g_var_tot - gxe_var_tot
  
  y_df <- map(seq(1, n_sim), function(idx) {
    e <- e_df[[paste0("e", idx)]]
    # if (nl_e) e <- exp(e) / (exp(1) - (exp(1) - 1))  # Exponential function for nonlinearity, normalize to variance one
    # if (nl_e) e <- (e + e^2) / sqrt(3)  # Square for nonlinearity, normalize to mean zero, variance one
    # if (nl_e) e <- as.vector(scale(e + e^2))  # Square for nonlinearity, normalize to mean zero, variance one
    y_mean_vec <- g_mat %*% beta_g_vec
    if (nl_e & e_distr == "normal") {
      y_mean_vec <- y_mean_vec + sqrt(e + 10) * (beta_e / sd(sqrt(e + 10)))
    } else if (nl_e & e_distr == "gamma") {
      y_mean_vec <- y_mean_vec + sqrt(e) * (beta_e / sd(sqrt(e)))
    } else {
      y_mean_vec <- y_mean_vec + e * beta_e
    }
    y_mean_vec <- y_mean_vec + (g_mat * e) %*% beta_gxe_vec
    rnorm(N, y_mean_vec, sqrt(error_var))
  }) %>%
    setNames(paste0("y", seq(1, n_sim))) %>%
    bind_cols()
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  system(paste0("mkdir -p ", target_dir))
  bind_cols(sim_df, e_df, y_df) %>%
    write_csv(paste0(target_dir, "/phenos.csv"))
}


simulate_phenotypes <- function(scenario_df,
                                g_mat, maf_vec,
                                train_prop,
                                seed = 1) {
  
  set.seed(seed)
  
  scenario_df %>%
    rowwise() %>%
    group_walk(~ simulate_scenario(.$ge_var_tot, .$e_distr, .$e_icc,
                                   .$e_var_tot, .$g_var_tot, .$gxe_var_tot, 
                                   .$nl_e,
                                   .$n_sim, .$tag, 
                                   g_mat, maf_vec, train_prop))
}


test_g <- function(g, e, y, df, covars = NULL) {
  regression_df <- df %>%
    rename(g = {{ g }}, y = {{ y }})
  lm_form_str <- "y ~ g"
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  lm(as.formula(lm_form_str), data = regression_df) %>%
    broom::tidy() %>%
    filter(term == "g")
}


test_g_by_e <- function(g, e, y, df, covars = NULL) {
  regression_df <- df %>%
    rename(g = {{ g }}, e_test = {{ e }}, y = {{ y }})
  lm_form_str <- "y ~ g * e_test"
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  lm(as.formula(lm_form_str), data = regression_df) %>%
    broom::tidy() %>%
    filter(term == "g:e_test")
}


test_vqtl <- function(g, e, y, df, covars = NULL) {
  regression_df <- df %>%
    rename(g = {{ g }}, y = {{ y }})
  if (!is.null(covars)) {
    covar_lm_form_str <- paste0("y ~ ", paste(covars, collapse = " + "))
    covar_lm <- lm(as.formula(covar_lm_form_str), data = regression_df)
    regression_df$y <- resid(covar_lm)
  }
  regression_df <- regression_df %>%
    group_by(g) %>%
    mutate(z = abs(y - median(y)))
  lm_form_str <- "z ~ g"
  lm(as.formula(lm_form_str), data = regression_df) %>%
    broom::tidy() %>%
    filter(term == "g")
}


test_all_variants <- function(e_name, y_name, pheno_df, g_mat, regression_func) {
  map(colnames(g_mat), function(g_name) {
    pheno_df$g <- g_mat[, g_name]
    regression_func("g", e_name, y_name,
                    filter(pheno_df, train == 1))
  }) %>%
    setNames(colnames(g_mat)) %>%
    bind_rows(.id = "g")
}


test_associations <- function(tag, g_mat, n_cores = 1) {
  print(paste0("Testing variant associations: ", tag, "..."))
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  sim_idx_vec <- gsub("y", "", all_phenos)
  
  silent <- parallel::mclapply(sim_idx_vec, function(idx) {
    e_name <- paste0("e_m", idx)
    y_name <- paste0("y", idx)
    test_all_variants(e_name, y_name, pheno_df, g_mat, test_g) %>%
      write_csv(paste0(target_dir, "/gwas_res_", y_name, ".csv"))
    test_all_variants(e_name, y_name, pheno_df, g_mat, test_g_by_e) %>%
      write_csv(paste0(target_dir, "/igwas_res_", y_name, ".csv"))
    test_all_variants(e_name, y_name, pheno_df, g_mat, test_vqtl) %>%
      write_csv(paste0(target_dir, "/vgwas_res_", y_name, ".csv"))
  }, mc.cores = n_cores)
}


get_pgs_weights <- function(tag, p_thresh = 0.05) {
  print(paste0("Deriving PGS weights: ", tag, "..."))
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  
  walk(all_phenos, function(y_name) {
    read_csv(paste0(target_dir, "/gwas_res_", y_name, ".csv"),
             show_col_types = FALSE) %>%
      mutate(beta = ifelse(p.value < p_thresh, estimate, 0)) %>%
      select(g, beta) %>%
      write_csv(paste0(target_dir, "/mpgs_weights_", y_name, ".csv"))
    read_csv(paste0(target_dir, "/igwas_res_", y_name, ".csv"),
             show_col_types = FALSE) %>%
      mutate(beta = ifelse(p.value < p_thresh, estimate, 0)) %>%
      select(g, beta) %>%
      write_csv(paste0(target_dir, "/ipgs_weights_", y_name, ".csv"))
    read_csv(paste0(target_dir, "/vgwas_res_", y_name, ".csv"),
             show_col_types = FALSE) %>%
      mutate(beta = ifelse(p.value < p_thresh, estimate, 0)) %>%
      select(g, beta) %>%
      write_csv(paste0(target_dir, "/vpgs_weights_", y_name, ".csv"))
  })
}


calculate_pgs <- function(tag, g_mat) {
  print(paste0("Calculating PGS: ", tag, "..."))
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  sim_idx_vec <- gsub("y", "", all_phenos)
  
  walk(sim_idx_vec, function(idx) {
    e_name <- paste0("e", idx)
    e_m_name <- paste0("e_m", idx)
    y_name <- paste0("y", idx)
    mpgs_weights_df <- read_csv(paste0(target_dir, "/mpgs_weights_", y_name, ".csv"),
                                show_col_types = FALSE)
    pheno_df$mpgs <- drop(g_mat[, mpgs_weights_df$g] %*% mpgs_weights_df$beta)
    ipgs_weights_df <- read_csv(paste0(target_dir, "/ipgs_weights_", y_name, ".csv"),
                                show_col_types = FALSE)
    pheno_df$ipgs <- drop(g_mat[, ipgs_weights_df$g] %*% ipgs_weights_df$beta)
    vpgs_weights_df <- read_csv(paste0(target_dir, "/vpgs_weights_", y_name, ".csv"),
                                show_col_types = FALSE)
    pheno_df$vpgs <- drop(g_mat[, vpgs_weights_df$g] %*% vpgs_weights_df$beta)
    pheno_df %>%
      select(id, train, e = {{ e_name }}, e_m = {{ e_m_name }}, y = {{ y_name }}, mpgs, ipgs, vpgs) %>%
      write_csv(paste0(target_dir, "/all_pgs_", y_name, ".csv"))
  })
}


test_pgs <- function(tag, g_mat) {
  print(paste0("Testing PGS: ", tag, "..."))
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  
  test_pgs_by_e <- function(pgs_type, e, y, covars = NULL, 
                            test_nl_e = FALSE, robust_SE = FALSE) {
    regression_df <- read_csv(paste0(target_dir, "/all_pgs_", y, ".csv"),
                              show_col_types = FALSE) %>%
      rename(e_test = {{ e }})
    lm_form_str <- paste0("y ~ ", pgs_type, " * e_test")
    if (!is.null(covars)) {
      lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
    }
    if (test_nl_e) lm_form_str <- paste0(lm_form_str, " + ns(e_test, df = 10)")
    lm_fit <- lm(as.formula(lm_form_str), data = filter(regression_df, train == 0)) 
    if (robust_SE) {
      lm_fit <- lmtest::coeftest(lm_fit, vcov = sandwich::vcovHC, type = "HC0")
    }
    lm_fit %>%
      broom::tidy() %>%
      filter(grepl(pgs_type, term))
  }
  
  pgs_types <- c("mpgs", "ipgs", "vpgs")
  sim_res_df <- map(all_phenos, function(y_name) {
    map(pgs_types, function(pgs_type) {
      test_pgs_by_e(pgs_type, "e_m", y_name)
    }) %>%
      setNames(pgs_types) %>%
      bind_rows(.id = "pgs_type")
  }) %>%
    bind_rows()
  sim_res_df %>%
    write_csv(paste0(target_dir, "/test_res.csv"  ))
  
  # sim_res_df <- map(all_phenos, function(y_name) {
  #   map(pgs_types, function(pgs_type) {
  #     test_pgs_by_e(pgs_type, "e_m", y_name, test_nl_e = TRUE)
  #   }) %>%
  #     setNames(pgs_types) %>%
  #     bind_rows(.id = "pgs_type")
  # }) %>%
  #   bind_rows()
  # sim_res_df %>%
  #   write_csv(paste0(target_dir, "/test_res_nlE.csv"))
  
  # sim_res_df <- map(all_phenos, function(y_name) {
  #   map(pgs_types, function(pgs_type) {
  #     test_pgs_by_e(pgs_type, "e_m", y_name, robust_SE = TRUE)
  #   }) %>%
  #     setNames(pgs_types) %>%
  #     bind_rows(.id = "pgs_type")
  # }) %>%
  #   bind_rows()
  # sim_res_df %>%
  #   write_csv(paste0(target_dir, "/test_res_robustSE.csv"))
}
