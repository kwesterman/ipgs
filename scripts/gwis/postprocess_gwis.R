library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
filepath <- args[1]
e <- args[2]
y <- args[3]

gwis_dir <- dirname(filepath)


### Define functions

calc_lambda <- function(x, p=0.5){
  # Calculate genomic inflation lambda value
  x = x[!is.na(x)]
  x.quantile <- quantile(x, p)
  round(qchisq(1 - x.quantile, 1) / qchisq(p, 1), 2)
}

make_qq <- function(data, pval_col, main=""){
  # Make a quantile-quantile plot
  data <- filter(data, data[[pval_col]] > 0)  # In case extremely low p-values are stored as zero
  
  # Process p-values
  y_vals <- sort(-log10(data[[pval_col]]))
  x_vals <- -log10(rev(ppoints(length(y_vals))))  # ppoints generates a uniform probability distribution
  
  # Trim points at higher p-values (credit to RaMWAS package for code snippet)
  levels = as.integer((x_vals - x_vals[1]) / (tail(x_vals, 1) - x_vals[1]) * 2000)
  keep = c(TRUE, diff(levels) != 0)
  levels = as.integer((y_vals - y_vals[1])/(tail(y_vals, 1) - y_vals[1]) * 2000)
  keep = keep | c(TRUE, diff(levels) != 0)
  keep = which(keep)
  
  par(ps = 18)
  plot(x = x_vals[keep], y = y_vals[keep], 
       xlab = expression(-log[10](italic(p)) * " (Expected)"), 
       ylab = expression(-log[10](italic(p)) * " (Observed)"),
       main = main, cex = 0.8, 
       cex.lab = 0.8, cex.main = 0.9, 
       pch = 16, ylim = c(0, ceiling(max(y_vals))))
  abline(0, 1, lty = 2)
  legend(x = 'topleft', y = 'topleft',
         bquote(lambda == .(calc_lambda(data[[pval_col]]))), 
         cex = 0.9, bty = "n")
}


### Read in summary stats and subset to columns of interest

ss_cols <- c(
  CHR = "CHR", SNP = "RSID", POS = "POS",
  EA = "Effect_Allele", NEA = "Non_Effect_Allele", AF = "AF",
  N = "N_Samples", P_int = "P_Value_Interaction", P_joint = "P_Value_Joint", 
  robust_P_int = "robust_P_Value_Interaction", robust_P_joint = "robust_P_Value_Joint"
)

high_qual_variants <- read_tsv("../data/processed/ukb_rsIDs_maf0.005_info0.5.txt", 
                               col_names=F, col_types="c")[[1]]

ss_df <- fread(filepath, stringsAsFactors=F, data.table=F) %>%
  select(all_of(ss_cols), matches("^Beta"), matches("^Var_Beta")) %>%
  mutate(across(contains("P_"), ~ as.numeric(.))) %>%
  filter(SNP %in% high_qual_variants)


### Prepare files for downstream analysis

beta_col = paste0("Beta_G-", e)
prs_dir <- paste0(gwis_dir, "/../prs")
ss_df %>%
  filter(robust_P_int < 0.05) %>%
  select(SNP, CHR, POS, EA, NEA, beta_int = {{beta_col}}, P_int, robust_P_int) %>% 
  group_by(SNP) %>%
  arrange(robust_P_int) %>%
  slice(1) %>%
  ungroup() %>%
  write_tsv(paste0(prs_dir, "/", e, "_", y, "_prsInput"))
  #filter(robust_P_int == min(robust_P_int)) %>%

ldsc_dir <- paste0(gwis_dir, "/../ldsc")
ss_df %>%
  select(SNP, CHR, POS, EA, NEA, AF, N, Beta = {{beta_col}}, P = robust_P_int) %>% 
  group_by(SNP) %>%
  arrange(P) %>%
  slice(1) %>%
  ungroup() %>%
  write_tsv(paste0(ldsc_dir, "/", e, "_", y, "_ldscInput"))
  #filter(P == min(P)) %>%


### Create Q-Q plot
  
qq_dir <- paste0(dirname(filepath), "/qq_plots/")
system(paste0("mkdir -p ", qq_dir))
write(calc_lambda(ss_df$robust_P_int), paste0(qq_dir, gsub("_merged|\\.tbl", "_lambda", basename(filepath))))
plot_filepath <- paste0(qq_dir, gsub("_merged|\\.tbl", "_QQ.pdf", basename(filepath)))
pdf(file = plot_filepath)
make_qq(ss_df, "P_int")
dev.off()
robust_plot_filepath <- paste0(qq_dir, gsub("_merged|\\.tbl", "_robust_QQ.pdf", basename(filepath)))
pdf(file = robust_plot_filepath)
make_qq(ss_df, "robust_P_int")
dev.off()
