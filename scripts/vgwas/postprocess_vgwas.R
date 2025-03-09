library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
filepath <- args[1]
pheno <- args[2]
#tag <- args[3]

vgwas_dir <- dirname(filepath)
MAF <- 0.005

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
  CHR = "CHR", SNP = "SNP", POS = "BP",
  EA = "A1", NEA = "A2",
  N = "N", P = "P", beta="BETA", SE="SE")

high_qual_variants <- read_tsv("../data/processed/ukb_rsIDs_maf0.005_info0.5.txt", col_names=F, col_types="c")[[1]]

ss_df <- fread(filepath, stringsAsFactors=F, data.table=F) %>%
  select(all_of(ss_cols)) %>%
  mutate_at("P", ~ as.numeric(.)) %>%
  filter(SNP %in% high_qual_variants)


### Prepare files for downstream analysis

pgs_dir <- paste0(vgwas_dir, "/../pgs")
system(paste0("mkdir -p ", pgs_dir))
ss_df %>%
  filter(P < 0.05) %>%
  select(SNP, CHR, POS, EA, NEA, beta, P) %>% 
  group_by(SNP) %>% 
  arrange(P) %>%
  slice(1) %>% 
  ungroup() %>% 
  write_tsv(paste0(pgs_dir, "/", pheno, "_vQTL_pgsInput"))
 

ldsc_dir <- paste0(vgwas_dir, "/../ldsc")
system(paste0("mkdir -p ", ldsc_dir))
ss_df %>%
  select(SNP, CHR, POS, EA, NEA, N, beta, P) %>% 
  group_by(SNP) %>%
  arrange(P) %>%
  slice(1) %>%
  ungroup() %>% 
  write_tsv(paste0(ldsc_dir, "/", pheno, "_vQTL_ldscInput"))



### Create Q-Q plot
 
qq_dir <- paste0(dirname(filepath), "/qq_plots/")
system(paste0("mkdir -p ", qq_dir))

write(calc_lambda(ss_df$P), paste0(qq_dir, gsub("_merged|\\.tbl", "_lambda", basename(filepath))))
plot_filepath <- paste0(qq_dir, gsub("_merged|\\.tbl", "_QQ.pdf", basename(filepath)))
pdf(file = plot_filepath)
make_qq(ss_df, "P")
dev.off()


### Delete geno file
#system(paste0("rm ../data/processed/", pheno, "_ukb_geno_1to22_maf", MAF, "_info0.5.txt"))


##EOF

