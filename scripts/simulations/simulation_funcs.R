library(tidyverse)


simulate_genotypes <- function(N, M, min_maf, max_maf, 
                               seed = 1) {
  set.seed(seed)
  
  mafs <- runif(M, min_maf, max_maf)
  
  g_mat <- sapply(mafs, function(maf) {
    rbinom(N, 2, maf)
  }, simplify = TRUE)
  colnames(g_mat) <- paste0("g", seq(1, M))
  
  saveRDS(g_mat, paste0(sim_dir, "/g_mat.rds"))
  saveRDS(mafs, paste0(sim_dir, "/maf_vec.rds"))
}


simulate_exposure <- function(g_mat, dist_E = "normal", ge_var_tot = 0, icc = 1, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  
  N <- nrow(g_mat)
  
  e_mat <- matrix(rnorm(N), N, 1)  # Base exposure matrix (standard normal)
  colnames(e_mat) <- "e1"
  
  pgs <- as.vector(scale(scale(g_mat) %*% rnorm(ncol(g_mat))))  # "Polygenic score" with variance=1 (scale to prevent MAF-dependent contribution)
  pgs_mat <- matrix(pgs, nrow = N, ncol = 1, byrow = FALSE)  # Same PGS for each simulation replicate
  e_mat <- sqrt(1 - ge_var_tot) * e_mat + sqrt(ge_var_tot) * pgs_mat  # Add "PGS" effects based on G-E correlation
  
  noise <- matrix(rnorm(N), N, 1)    # Add any measurement error
  e_m_mat <- sqrt(icc) * e_mat + sqrt(1 - icc) * noise
  colnames(e_m_mat) <- "e_m1"

  if (dist_E == "normal10") {
    e_mat <- e_mat + 10
    e_m_mat <- e_m_mat + 10
  }
  
  if (dist_E == "gamma") {  # Use copula approach to generate gamma-distributed exposure
    U <- pnorm(e_mat)  # Push through normal CDF
    e_mat <- qgamma(U, shape = 1, rate = 1)  # Generate marginal gamma using quantile function
    U_m <- pnorm(e_m_mat)
    e_m_mat <- qgamma(U_m, shape = 1, rate = 1)
  }
  
  bind_cols(e_mat, e_m_mat)
}


simulate_scenario <- function(ge_var_tot, e_distr, e_icc,
                              e_var_tot, g_var_tot, gxe_var_tot, nl_e,
                              tag,
                              g_mat, maf_vec,
                              train_prop) {
  
  N <- nrow(g_mat)
  M <- ncol(g_mat)
  
  sim_df <- tibble(
    id = seq(1, N),
    train = sample(c(0, 1), N, replace = TRUE, 
                   prob = c(1 - train_prop, train_prop))
  )
  
  e_df <- simulate_exposure(g_mat, e_distr, ge_var_tot, e_icc)
  e <- e_df$e1
  

  
  g_mat_variances <- 2 * maf_vec * (1 - maf_vec)
  
  g_var_single <- g_var_tot / M
  beta_g_vars <- g_var_single / g_mat_variances
  beta_g_vec <- rnorm(M, 0, sqrt(beta_g_vars))
  y_mean_vec <- g_mat %*% beta_g_vec
  
  beta_e <- sqrt(e_var_tot)
  if (nl_e & e_distr == "normal") {
    y_mean_vec <- y_mean_vec + sqrt(e + 10) * (beta_e / sd(sqrt(e + 10)))
  } else if (nl_e & e_distr == "gamma") {
    y_mean_vec <- y_mean_vec + sqrt(e) * (beta_e / sd(sqrt(e)))
  } else {
    y_mean_vec <- y_mean_vec + e * beta_e
  }
  
  gxe_var_single <- gxe_var_tot / M
  beta_gxe_vars <- gxe_var_single / g_mat_variances
  beta_gxe_vec <- rnorm(M, 0, sqrt(beta_gxe_vars))
  y_mean_vec <- y_mean_vec + (g_mat * e) %*% beta_gxe_vec
  
  error_var <- 1 - e_var_tot - g_var_tot - gxe_var_tot
  y <- rnorm(N, y_mean_vec, sqrt(error_var))
  
  y_df <- tibble(y1 = y)
  
  bind_cols(sim_df, e_df, y_df)
}


test_g <- function(g_name, y_name, df, covars = NULL) {
  regression_df <- df %>%
    mutate(g = .data[[g_name]], y = .data[[y_name]])
  lm_form_str <- paste0("y ~ ", g_name)
  lm_form_str <- "y ~ g"
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  lm(as.formula(lm_form_str), data = regression_df) %>%
    broom::tidy() %>%
    filter(term == "g") %>%
    mutate(g = g_name,
           test = "g") %>%
    select(test, g, estimate, p.value)
}


test_gxe <- function(g_name, e_name, y_name, df, covars = NULL) {
  regression_df <- df %>%
    mutate(g = .data[[g_name]], e_test = .data[[e_name]], y = .data[[y_name]])
  lm_form_str <- "y ~ g * e_test"
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  lm(as.formula(lm_form_str), data = regression_df) %>%
    broom::tidy() %>%
    filter(term == "g:e_test") %>%
    mutate(g = g_name,
           test = "gxe") %>%
    select(test, g, estimate, p.value)
}


test_vqtl <- function(g_name, y_name, df, covars = NULL) {
  regression_df <- df %>%
    mutate(g = .data[[g_name]], y = .data[[y_name]])
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
    filter(term == "g") %>%
    mutate(g = g_name,
           test = "vqtl") %>%
    select(test, g, estimate, p.value)
}


test_variants <- function(e_name, y_name, pheno_df, g_mat, covars = NULL) {
  map_dfr(
    .x = colnames(g_mat),
    .f = function(g_name) {
      df <- pheno_df %>%
        bind_cols(g_mat[, g_name, drop = FALSE]) %>%
        filter(train == 1)
      
      g_res  <- test_g(g_name, y_name, df, covars)
      gxe_res <- test_gxe(g_name, e_name, y_name, df, covars)
      vqtl_res  <- test_vqtl(g_name, y_name, df, covars)
      
      bind_rows(g_res, gxe_res, vqtl_res)
    }
  )
}


get_pgs_weights <- function(assoc_df, p_thresh = 0.05) {
  assoc_df %>%
    mutate(beta = ifelse(p.value < p_thresh, estimate, 0))
}


calculate_pgs <- function(g_mat, weights_df) {
  stopifnot(names(g_mat) == weights_df$g)
  
  drop(g_mat %*% weights_df$beta)
}


test_pgs_by_e <- function(pgs_type, e_name, y_name, df, covars = NULL, 
                          test_nl_e = FALSE, robust_SE = FALSE) {
  df <- df %>%
    mutate(pgs = .data[[pgs_type]], e_test = .data[[e_name]], y = .data[[y_name]])
  if (all(df$pgs == 0)) return(tibble(estimate = 0, p.value = 1))
  lm_form_str <- "y ~ pgs * e_test"
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  if (test_nl_e) lm_form_str <- paste0(lm_form_str, " + splines::ns(e_test, df = 10)")
  lm_fit <- lm(as.formula(lm_form_str), data = filter(df, train == 0)) 
  if (robust_SE) {
    lm_fit <- lmtest::coeftest(lm_fit, vcov = sandwich::vcovHC, type = "HC0")
  }
  lm_fit %>%
    broom::tidy() %>%
    filter(term == "pgs:e_test") %>%
    select(estimate, p.value)
}


process_one_rep <- function(scn, rep, g_mat, maf_vec) {
  
  pheno_df <- simulate_scenario(  # Simulate E and Y
    scn$ge_var_tot, scn$e_distr, scn$e_icc,
    scn$e_var_tot, scn$g_var_tot, scn$gxe_var_tot, scn$nl_e,
    scn$tag,
    g_mat, maf_vec,
    train_prop = 0.7
  )
  
  assoc_tbl <- test_variants("e_m1", "y1", pheno_df, g_mat)  # Run each type of association test
  
  weights_df <- get_pgs_weights(assoc_tbl)  # Derive PGS weights
  
  pheno_df$mpgs <- calculate_pgs(g_mat, filter(weights_df, test == "g"))  # Calculate PGS
  pheno_df$ipgs <- calculate_pgs(g_mat, filter(weights_df, test == "gxe"))
  pheno_df$vpgs <- calculate_pgs(g_mat, filter(weights_df, test == "vqtl"))
  
  pgs_types <- c("mpgs", "ipgs", "vpgs")
  map_dfr(setNames(pgs_types, pgs_types), function(pgs_type) {  # Test PGSxE
    test_pgs_by_e(pgs_type, "e_m1", "y1", pheno_df)
  }, .id = "pgs_type")
}
