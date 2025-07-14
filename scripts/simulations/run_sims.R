library(tidyverse)
library(furrr)


args <- commandArgs(trailingOnly = TRUE)
n_cpus <- as.numeric(args[1])


n_workers <- if (!is.na(n_cpus)) n_cpus - 1 else parallel::detectCores() - 1
plan(multisession, workers = n_workers)

seed <- 1

sim_dir <- "../data/processed/simulations"

source("simulations/simulation_funcs.R")


# Simulate genotypes

simulate_genotypes(N = 10000, M = 100, min_maf = 0.01, max_maf = 0.5, seed = seed)
g_mat <- readRDS(paste0(sim_dir, "/g_mat.rds"))
maf_vec <- readRDS(paste0(sim_dir, "/maf_vec.rds"))


# Run type I error simulations

t1e_scenario_df <- read_csv(paste0(sim_dir, "/t1e_scenarios.csv"),
                            show_col_types = FALSE)
n_sim_t1e <- 5000

set.seed(seed)

t1e_res_df <- map_dfr(
  split(distinct(t1e_scenario_df, tag, .keep_all = TRUE), 
        seq(1, length(unique(t1e_scenario_df$tag)))),
  function(scn) {
    message("Testing TIE scenario: ", scn$tag)
    sim_res_df <- future_map_dfr(seq_len(n_sim_t1e), function(rep) {
      process_one_rep(scn, rep, g_mat, maf_vec)
    }, .options = furrr_options(seed = TRUE))
    scenario_res_df <- sim_res_df %>%
      mutate(tag = scn$tag) %>%
      group_by(tag, pgs_type) %>%
      summarise(f = mean(p.value < 0.05),
                n_sim = n_sim_t1e,
                .groups = "drop")
    write_csv(scenario_res_df, 
              paste0("../data/processed/simulations/results_", scn$tag, ".csv"))
    scenario_res_df
  }
)

write_csv(t1e_res_df, "../data/processed/simulations/t1e_res.csv")


# Run power simulations

power_scenario_df <- read_csv(paste0(sim_dir, "/power_scenarios.csv"),
                              show_col_types = FALSE)
n_sim_power <- 500

set.seed(seed)

power_res_df <- map_dfr(
  split(distinct(power_scenario_df, tag, .keep_all = TRUE), 
        seq(1, length(unique(power_scenario_df$tag)))),
  function(scn) {
    message("Testing power scenario: ", scn$tag)
    sim_res_df <- future_map_dfr(seq_len(n_sim_power), function(rep) {
      process_one_rep(scn, rep, g_mat, maf_vec)
    }, .options = furrr_options(seed = TRUE))
    scenario_res_df <- sim_res_df %>%
      mutate(tag = scn$tag) %>%
      group_by(tag, pgs_type) %>%
      summarise(f = mean(p.value < 0.05),
                n_sim = n_sim_power,
                .groups = "drop")
    write_csv(scenario_res_df, 
              paste0("../data/processed/simulations/results_", scn$tag, ".csv"))
    scenario_res_df
  }
)

write_csv(power_res_df, "../data/processed/simulations/power_res.csv")