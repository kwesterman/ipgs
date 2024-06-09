export_mpgs <- function(export_tag) {
  opt_param_df <- read_csv("../data/processed/pgs/optimized_mpgs_params.csv")
  thresh_val <- opt_param_df$threshold[opt_param_df$tag == export_tag]
  snp_df <- read_tsv(paste0("../data/processed/pgs/", export_tag, ".snp"), 
                     show_col_types = FALSE) %>%
    filter(P < thresh_val) %>%
    select(SNP)
  pgs_weights_df <- read_tsv(paste0("../data/processed/pgs/", export_tag, "_pgsInput"), 
                             show_col_types = FALSE) %>%
    inner_join(snp_df, by = "SNP")
  pgs_weights_df
}

export_ipgs <- function(export_tag) {
  opt_param_df <- read_csv("../data/processed/pgs/optimized_ipgs_params.csv")
  thresh_val <- opt_param_df$threshold[opt_param_df$tag == export_tag]
  snp_df <- read_tsv(paste0("../data/processed/pgs/", export_tag, ".snp"), 
                     show_col_types = FALSE) %>%
    filter(P < thresh_val) %>%
    select(SNP)
  pgs_weights_df <- read_tsv(paste0("../data/processed/pgs/", export_tag, "_pgsInput"), 
                             show_col_types = FALSE) %>%
    inner_join(snp_df, by = "SNP")
  pgs_weights_df
}

export_vpgs <- function(export_tag) {
  opt_param_df <- read_csv("../data/processed/pgs/optimized_vpgs_params.csv")
  thresh_val <- opt_param_df$threshold[opt_param_df$tag == export_tag]
  snp_df <- read_tsv(paste0("../data/processed/pgs/", export_tag, ".snp"), 
                     show_col_types = FALSE) %>%
    filter(P < thresh_val) %>%
    select(SNP)
  pgs_weights_df <- read_tsv(paste0("../data/processed/pgs/", export_tag, "_pgsInput"), 
                             show_col_types = FALSE) %>%
    inner_join(snp_df, by = "SNP")
  pgs_weights_df
}