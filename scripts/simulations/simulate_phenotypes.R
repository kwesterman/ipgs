library(tidyverse)


train_prop <- 0.7


g_mat <- readRDS("../data/processed/simulations/g_mat.rds")
N <- nrow(g_mat)
M <- ncol(g_mat)

set.seed(1)

sim_df <- tibble(
  id = seq(1, N),
  train = sample(c(0, 1), N, replace = TRUE, 
                 prob = c(1 - train_prop, train_prop)),
  e = rnorm(N, 0, 1)
)

simulate_scenario <- function(beta_e, beta_g_var, beta_gxe_var, 
                              beta_prob_nonzero,
                              tag) {
  print(paste0("Simulating ", tag, "..."))
  
  beta_g_vec <- rnorm(M, 0, sqrt(beta_g_var)) * (runif(M) < beta_prob_nonzero)
  beta_gxe_vec <- rnorm(M, 0, sqrt(beta_gxe_var)) * (runif(M) < beta_prob_nonzero)
  
  y_df <- map(seq(1, M), function(pheno_idx) {
    y_mean_vec <- sim_df$e * beta_e +
      g_mat %*% beta_g_vec +
      (g_mat * sim_df$e) %*% beta_gxe_vec
    rnorm(N, y_mean_vec, 1)
  }) %>%
    setNames(paste0("y", seq(1, M))) %>%
    bind_cols()
  
  target_dir <- paste0("../data/processed/simulations/", tag)
  system(paste0("mkdir -p ", target_dir))
  bind_cols(sim_df, y_df) %>%
    write_csv(paste0(target_dir, "/phenos.csv"))
}

scenario_df <- read_csv("../data/processed/simulations/scenarios.csv") %>%
  rowwise() %>%
  group_walk(~ simulate_scenario(.$e, .$g, .$gxe, .$prob_nonzero, .$tag))
