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
}


simulate_scenario <- function(e_var_tot, g_var_tot, gxe_var_tot, 
                              beta_prob_nonzero, 
                              n_sim, tag,
                              g_mat, train_prop) {
  print(paste0("Simulating ", tag, "..."))
  
  N <- nrow(g_mat)
  M <- ncol(g_mat)
  
  sim_df <- tibble(
    id = seq(1, N),
    train = sample(c(0, 1), N, replace = TRUE, 
                   prob = c(1 - train_prop, train_prop)),
    e = rnorm(N, 0, 1)
  )
  
  beta_e <- sqrt(e_var_tot)
  
  beta_g_var <- g_var_tot / M
  beta_g_vec <- rnorm(M, 0, sqrt(beta_g_var)) * (runif(M) < beta_prob_nonzero)
  
  beta_gxe_var <- gxe_var_tot / M
  beta_gxe_vec <- rnorm(M, 0, sqrt(beta_gxe_var)) * (runif(M) < beta_prob_nonzero)
  
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
                                g_mat,
                                train_prop,
                                seed = 1) {
  
  set.seed(seed)
  
  scenario_df %>%
    rowwise() %>%
    group_walk(~ simulate_scenario(.$e_var_tot, .$g_var_tot, .$gxe_var_tot, 
                                   .$prob_nonzero, 
                                   .$n_sim, .$tag, g_mat, train_prop))
}


test_g_by_e <- function(g, e, y, df, covars = NULL) {
  regression_df <- df %>%
    rename(g = {{ g }}, e = {{ e }}, y = {{ y }})
  lm_form_str <- "y ~ g * e"
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  lm(as.formula(lm_form_str), data = regression_df) %>%
    broom::tidy() %>%
    filter(term == "g:e")
}


test_g <- function(g, e, y, df, covars = NULL) {
  regression_df <- df %>%
    rename(g = {{ g }}, e = {{ e }}, y = {{ y }})
  lm_form_str <- "y ~ g + e"
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  lm(as.formula(lm_form_str), data = regression_df) %>%
    broom::tidy() %>%
    filter(term == "g")
}


test_all_variants <- function(y_name, pheno_df, g_mat, regression_func) {
  map(colnames(g_mat), function(g_name) {
    pheno_df$g <- g_mat[, g_name]
    regression_func("g", "e", y_name,
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
    print(paste0(y_name, "..."))
    test_all_variants(y_name, pheno_df, g_mat, test_g_by_e) %>%
      write_csv(paste0(target_dir, "/gwis_res_", y_name, ".csv"))
    test_all_variants(y_name, pheno_df, g_mat, test_g) %>%
      write_csv(paste0(target_dir, "/gwas_res_", y_name, ".csv"))
  }, mc.cores = n_cores)
}


get_pgs_weights <- function(tag) {
  print(paste0("Deriving PGS weights: ", tag, "..."))
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  
  walk(all_phenos, function(y_name) {
    read_csv(paste0(target_dir, "/gwis_res_", y_name, ".csv"),
             show_col_types = FALSE) %>%
      select(g, beta = estimate) %>%
      write_csv(paste0(target_dir, "/ipgs_weights_", y_name, ".csv"))
    read_csv(paste0(target_dir, "/gwas_res_", y_name, ".csv"),
             show_col_types = FALSE) %>%
      select(g, beta = estimate) %>%
      write_csv(paste0(target_dir, "/pgs_weights_", y_name, ".csv"))
  })
}


calculate_pgs <- function(tag, g_mat) {
  print(paste0("Calculating PGS: ", tag, "..."))
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  
  walk(all_phenos, function(y_name) {
    ipgs_weights_df <- read_csv(paste0(target_dir, "/ipgs_weights_", y_name, ".csv"),
                                show_col_types = FALSE)
    pheno_df$ipgs <- drop(g_mat[, ipgs_weights_df$g] %*% ipgs_weights_df$beta)
    pgs_weights_df <- read_csv(paste0(target_dir, "/pgs_weights_", y_name, ".csv"),
                               show_col_types = FALSE)
    pheno_df$pgs <- drop(g_mat[, pgs_weights_df$g] %*% pgs_weights_df$beta)
    pheno_df %>%
      select(id, train, e, y = {{ y_name }}, ipgs, pgs) %>%
      write_csv(paste0(target_dir, "/all_pgs_", y_name, ".csv"))
  })
}


test_pgs <- function(tag, g_mat) {
  print(paste0("Testing PGS: ", tag, "..."))
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  
  test_pgs_by_e <- function(y, pgs_type, covars = NULL) {
    regression_df <- read_csv(paste0(target_dir, "/all_pgs_", y, ".csv"),
                              show_col_types = FALSE)
    lm_form_str <- paste0("y ~ ", pgs_type, " * e")
    if (!is.null(covars)) {
      lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
    }
    lm(as.formula(lm_form_str), data = filter(regression_df, train == 0)) %>%
      broom::tidy() %>%
      filter(term == paste0(pgs_type, ":e"))
  }
  
  pgs_types <- c("ipgs", "pgs")
  sim_res_df <- map(all_phenos, function(y_name) {
    map(pgs_types, function(pgs_type) {
      test_pgs_by_e(y_name, pgs_type)
    }) %>%
      setNames(pgs_types) %>%
      bind_rows(.id = "pgs_type")
  }) %>%
    bind_rows()
  sim_res_df %>%
    write_csv(paste0(target_dir, "/test_res.csv"))
}
