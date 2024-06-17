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

export_ppgs <- function(export_tag) {
  opt_param_df <- read_csv("../data/processed/pgs/optimized_ppgs_params.csv")
  pathway_weights <- filter(opt_param_df, tag == export_tag) %>%
    select(-tag) %>%
    unlist(use.names = TRUE)
  prset_snp_df <- read_tsv(paste0("../data/processed/pgs/", export_tag, ".snp"),
                     show_col_types = FALSE) %>%
    select(SNP, contains("HALLMARK"))
  snp_weights <- as.matrix(select(prset_snp_df, -SNP)) %*% pathway_weights
  snp_df <- bind_cols(prset_snp_df["SNP"], 
                      tibble(opt_weights = as.vector(snp_weights)))
  pgs_weights_df <- read_tsv(paste0("../data/processed/pgs/", 
                                    gsub("_prset", "", export_tag), "_pgsInput"),
                             show_col_types = FALSE) %>%
    inner_join(snp_df, by = "SNP") %>%
    mutate(beta = beta * opt_weights) %>%
    select(-opt_weights)
  pgs_weights_df
}