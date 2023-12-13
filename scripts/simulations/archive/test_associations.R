library(tidyverse)


# args <- commandArgs(trailingOnly = TRUE)
# tag <- args[[1]]

g_mat <- readRDS("../data/processed/simulations/g_mat.rds")


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

test_all_variants <- function(g_mat, y_name, regression_func) {
  map(colnames(g_mat), function(g_name) {
    pheno_df$g <- g_mat[, g_name]
    regression_func("g", "e", y_name,
                    filter(pheno_df, train == 1))
  }) %>%
    setNames(colnames(g_mat)) %>%
    bind_rows(.id = "g")
}

tags <- scan("../data/processed/simulations/scenario_tags.txt", what = character())

for (tag in tags) {
  target_dir <- paste0("../data/processed/simulations/", tag)
  pheno_df <- read_csv(paste0(target_dir, "/phenos.csv"),
                       show_col_types = FALSE)
  all_phenos <- grep("^y\\d*", names(pheno_df), value = TRUE)
  
  silent <- parallel::mclapply(all_phenos, function(y_name) {
    print(paste0(y_name, "..."))
    test_all_variants(g_mat, y_name, test_g_by_e) %>%
      write_csv(paste0(target_dir, "/gwis_res_", y_name, ".csv"))
    test_all_variants(g_mat, y_name, test_g) %>%
      write_csv(paste0(target_dir, "/gwas_res_", y_name, ".csv"))
  }, mc.cores = 6)
}
