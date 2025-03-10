library(tidyverse)


source("simulations/simulation_funcs.R")


simulate_genotypes(N = 10000, M = 300, min_maf = 0.05, max_maf = 0.5, seed = 123)
g_mat <- readRDS("../data/processed/simulations/g_mat.rds")
maf_vec <- readRDS("../data/processed/simulations/maf_vec.rds")


# power_scenario_df <- read_csv("../data/processed/simulations/power_scenarios.csv",
#                               show_col_types = FALSE)
# power_tags <- power_scenario_df$tag
# 
# simulate_phenotypes(power_scenario_df, g_mat, maf_vec, train_prop = 0.7, seed = 1)
# 
# for (tag in power_tags) test_associations(tag, g_mat, n_cores = 16)
# 
# for (tag in power_tags) get_pgs_weights(tag, p_thresh = 0.05)
# 
# for (tag in power_tags) calculate_pgs(tag, g_mat)
# 
# for (tag in power_tags) test_pgs(tag, g_mat)


t1e_scenario_df <- read_csv("../data/processed/simulations/t1e_scenarios.csv",
                            show_col_types = FALSE)
# t1e_tags <- t1e_scenario_df$tag
t1e_tags <- filter(t1e_scenario_df, nl_e == TRUE)$tag

simulate_phenotypes(t1e_scenario_df, g_mat, maf_vec, train_prop = 0.7, seed = 1)

for (tag in t1e_tags) test_associations(tag, g_mat, n_cores = 16)

for (tag in t1e_tags) get_pgs_weights(tag, p_thresh = 0.05)

for (tag in t1e_tags) calculate_pgs(tag, g_mat)

for (tag in t1e_tags) test_pgs(tag, g_mat)
