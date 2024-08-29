library(tidyverse)


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


simulate_scenario <- function(e_var_tot, g_var_tot, gxe_var_tot, 
                              e_distr, e_icc,
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
  
  if (e_distr == "normal") {
    beta_e_em <- sqrt(e_icc)  # Such that true E explains amount of variance equal to ICC
    var_em <- 1 - e_icc  # Such that measured E has total variance equal to 1
    sim_df <- sim_df %>%
      mutate(
        e = rnorm(N, 0, 1),
        e_m = rnorm(N, beta_e_em * e, sqrt(var_em))
      )
  } else if (e_distr == "normal20") {
    sim_df <- sim_df %>%
      mutate(
        e = rnorm(N, 20, 1),
        e_m = e
      )
  } else if (e_distr == "gamma") {
    sim_df <- sim_df %>%
      mutate(
        e = rgamma(N, shape = 1, rate = 1),
        e_m = e
      )
  }
  
  beta_e <- sqrt(e_var_tot)
  
  g_mat_variances <- 2 * maf_vec * (1 - maf_vec)
  g_var_single <- g_var_tot / M
  beta_g_vec_unscaled <- rnorm(M, 0, sqrt(g_var_single))  # Assumes all genotype variances are 1
  beta_g_vec <- beta_g_vec_unscaled / sqrt(g_mat_variances)  # Scale because they're not

  gxe_var_single <- gxe_var_tot / M
  beta_gxe_vec_unscaled <- rnorm(M, 0, sqrt(gxe_var_single))
  beta_gxe_vec <- beta_gxe_vec_unscaled / sqrt(g_mat_variances)
  
  error_var <- 1 - e_var_tot - g_var_tot - gxe_var_tot
  
  y_df <- map(seq(1, n_sim), function(pheno_idx) {
    y_mean_vec <- sim_df$e * beta_e +
      g_mat %*% beta_g_vec +
      (g_mat * sim_df$e) %*% beta_gxe_vec
    rnorm(N, y_mean_vec, sqrt(error_var))
  }) %>%
    setNames(paste0("y", seq(1, n_sim))) %>%
    bind_cols()
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  system(paste0("mkdir -p ", target_dir))
  bind_cols(sim_df, y_df) %>%
    write_csv(paste0(target_dir, "/phenos.csv"))
}


simulate_phenotypes <- function(scenario_df,
                                g_mat, maf_vec,
                                train_prop,
                                seed = 1) {
  
  set.seed(seed)
  
  scenario_df %>%
    rowwise() %>%
    group_walk(~ simulate_scenario(.$e_var_tot, .$g_var_tot, .$gxe_var_tot, 
                                   .$e_distr, .$e_icc,
                                   .$n_sim, .$tag, g_mat, maf_vec, train_prop))
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


test_all_variants <- function(y_name, pheno_df, g_mat, regression_func) {
  map(colnames(g_mat), function(g_name) {
    pheno_df$g <- g_mat[, g_name]
    regression_func("g", "e_m", y_name,
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
  
  silent <- parallel::mclapply(all_phenos, function(y_name) {
    test_all_variants(y_name, pheno_df, g_mat, test_g) %>%
      write_csv(paste0(target_dir, "/gwas_res_", y_name, ".csv"))
    test_all_variants(y_name, pheno_df, g_mat, test_g_by_e) %>%
      write_csv(paste0(target_dir, "/igwas_res_", y_name, ".csv"))
    test_all_variants(y_name, pheno_df, g_mat, test_vqtl) %>%
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
  
  walk(all_phenos, function(y_name) {
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
      select(id, train, e, e_m, y = {{ y_name }}, mpgs, ipgs, vpgs) %>%
      write_csv(paste0(target_dir, "/all_pgs_", y_name, ".csv"))
  })
}


test_pgs <- function(tag, g_mat) {
  print(paste0("Testing PGS: ", tag, "..."))
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  
  test_pgs_by_e <- function(pgs_type, e, y, covars = NULL) {
    regression_df <- read_csv(paste0(target_dir, "/all_pgs_", y, ".csv"),
                              show_col_types = FALSE) %>%
      rename(e_test = {{ e }})
    lm_form_str <- paste0("y ~ ", pgs_type, " * e_test")
    if (!is.null(covars)) {
      lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
    }
    lm(as.formula(lm_form_str), data = filter(regression_df, train == 0)) %>%
      broom::tidy() %>%
      filter(grepl(pgs_type, term))
      # filter(term == paste0(pgs_type, ":e_test"))
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
    write_csv(paste0(target_dir, "/test_res.csv"))
}
