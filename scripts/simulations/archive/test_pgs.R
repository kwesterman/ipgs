library(tidyverse)


g_mat <- readRDS("../data/processed/simulations/g_mat.rds")


tags <- scan("../data/processed/simulations/scenario_tags.txt", what = character())

for (tag in tags) {
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
    print(paste0(y_name, "..."))
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
