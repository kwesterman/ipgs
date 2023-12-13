library(tidyverse)


g_mat <- readRDS("../data/processed/simulations/g_mat.rds")

tags <- scan("../data/processed/simulations/scenario_tags.txt", what = character())

for (tag in tags) {
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

