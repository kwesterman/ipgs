# Format phenotype file for QUAIL


# load packages
library(tidyverse)
library(data.table)

# system arguments
args = commandArgs(trailingOnly=TRUE)
phenoFile=args[1]
phenoName=args[2]

covars=c("sex", "age", "age_squared", "ageBySex", "gPC1", "gPC2", "gPC3", "gPC4", "gPC5", "gPC6", "gPC7", "gPC8", "gPC9", "gPC10")
print(covars)


# gather phenos
phenos=fread(paste0(phenoFile)) %>%
 select(FID=id, IID=id, all_of(phenoName), all_of(covars)) %>%
 filter(complete.cases(.))


# write files 
phenos %>%
 select(FID, IID, all_of(phenoName)) %>%
 fwrite(file=paste0("../data/vgwas/quail/ukb_", phenoName, "_quail.txt"), row.names=F, col.names=T, sep=' ')


phenos %>%
 select(FID, IID, all_of(covars)) %>%
 fwrite(file=paste0("../data/vgwas/quail/ukb_covars_quail.txt"), row.names=F, col.names=T, sep=' ')



#EOF
