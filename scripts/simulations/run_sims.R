library(tidyverse)


source("simulations/simulation_funcs.R")

# scenario_df <- expand_grid(
#   e_var_tot = c(0, 1),
#   g_var_tot = c(0, 1),
#   gxe_var_tot = seq(0, 1, 0.2),
#   prob_nonzero = 1,
#   n_sim = 200
# ) %>%
#   mutate(tag = paste0("e", e_var_tot, "_g", g_var_tot, "_gxe", gxe_var_tot))
# 
# scenario_df %>%
#   write_csv("../data/processed/simulations/scenarios.csv")
# tags <- scenario_df$tag
# write(tags, "../data/processed/simulations/scenario_tags.txt")


simulate_genotypes(N = 10000, M = 100, min_maf = 0.01, max_maf = 0.5)
g_mat <- readRDS("../data/processed/simulations/g_mat.rds")


power_scenario_df <- read_csv("../data/processed/simulations/power_scenarios.csv")
power_tags <- power_scenario_df$tag

simulate_phenotypes(power_scenario_df, g_mat, train_prop = 0.7)

for (tag in power_tags) test_associations(tag, g_mat, n_cores = 16)

for (tag in power_tags) get_pgs_weights(tag)

for (tag in power_tags) calculate_pgs(tag, g_mat)

for (tag in power_tags) test_pgs(tag, g_mat)


t1e_scenario_df <- read_csv("../data/processed/simulations/t1e_scenarios.csv")
t1e_tags <- t1e_scenario_df$tag

simulate_phenotypes(t1e_scenario_df, g_mat, train_prop = 0.7)

for (tag in t1e_tags) test_associations(tag, g_mat, n_cores = 16)

for (tag in t1e_tags) get_pgs_weights(tag)

for (tag in t1e_tags) calculate_pgs(tag, g_mat)

for (tag in t1e_tags) test_pgs(tag, g_mat)
