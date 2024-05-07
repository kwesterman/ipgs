library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
filepath <- args[1]
y <- args[2]

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
  N = "N", P_vgwas = "P", BETA="BETA", SE="SE")

high_qual_variants <- as.data.frame(read_tsv(paste0(vgwas_dir, "/", y, "_ukb_geno_1to22_maf", MAF, "_info0.5.txt"), col_names=F, col_types="c"))$X3

ss_df <- fread(filepath, stringsAsFactors=F, data.table=F) %>%
  select(all_of(ss_cols)) %>%
  mutate(across(contains("P_"), ~ as.numeric(.))) %>%
  filter(SNP %in% high_qual_variants)


### Prepare files for downstream analysis

beta_col = "BETA"
prs_dir <- paste0(vgwas_dir, "/../prs")
system(paste0("mkdir -p ", prs_dir))
ss_df %>%
  filter(P_vgwas < 0.05) %>%
  select(SNP, CHR, POS, EA, NEA, BETA, P_vgwas) %>% 
  group_by(SNP) %>% 
  arrange(P_vgwas) %>%
  slice(1) %>% 
  ungroup() %>% 
  write_tsv(paste0(prs_dir, "/", y, "_prsInput"))
 

ldsc_dir <- paste0(vgwas_dir, "/../../ldsc")
system(paste0("mkdir -p ", ldsc_dir))
ss_df %>%
  select(SNP, CHR, POS, EA, NEA, N, BETA, P=P_vgwas) %>% 
  group_by(SNP) %>%
  arrange(P) %>%
  slice(1) %>%
  ungroup() %>% 
  write_tsv(paste0(ldsc_dir, "/", y, "_ldscInput"))



### Create Q-Q plot
 
qq_dir <- paste0(dirname(filepath), "/qq_plots/")
system(paste0("mkdir -p ", qq_dir))

write(calc_lambda(ss_df$P_vgwas), paste0(qq_dir, gsub("_merged|\\.tbl", "_lambda", basename(filepath))))
plot_filepath <- paste0(qq_dir, gsub("_merged|\\.tbl", "_QQ.pdf", basename(filepath)))
pdf(file = plot_filepath)
make_qq(ss_df, "P_vgwas")
dev.off()


### Delete geno file
system(paste0("rm ", vgwas_dir, "/", y, "_ukb_geno_1to22_maf", MAF, "_info0.5.txt"))


##EOF
