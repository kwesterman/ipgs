library(tidyverse)


tags <- scan("../data/processed/simulations/scenario_tags.txt", what = character())

for (tag in tags) {
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
