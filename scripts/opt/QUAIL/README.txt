README for QUAIL scripts (JEG)

* Minor edits made to original Step-2 QUAIL script on line 19: Add heading for "ERROR" in plink glm.linear output: 
** Original: colnames(df) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "TEST", "N", "BETA", "SE", "Z", "P", "A2")
** _JEG: colnames(df) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "TEST", "N", "BETA", "SE", "Z", "P", "ERROR", "A2")
** The resulting file was renamed as: Step2_QUAIL_vQTL_JEG.R


* Preliminary tests for running QUAIL with the covariates listed below errored out, as follows: 
* covariates: age, sex, age_squared, ageBySex, gPC1, gPC2, gPC3, gPC4, gPC5, gPC6, gPC7, gPC8. gPC9, gPC10

** Error 1-produced when using 14 raw covariates*: Error: Cannot proceed with --glm regression on phenotype 'int_rank_score', since genotype/covariate scales vary too widely for numerical stability of the current implementation. Try rescaling your covariates 	with e.g. --covar-variance-standardize.
*** In response, '--covar-variance-standardie ' was added to the plink2 commands.

** Error 2-produced when using 14 raw covariates with added flag for '--covar-variance-standardize '*: Error: Cannot proceed with --glm 	regression on phenotype 'int_rank_score', since variance inflation factor for covariate 'age' is too high (VIF_TOO_HIGH). You may want to remove redundant covariates and try again.
*** In response, '--vif 300 ' was added to the plink2 commands (NOTE: the default is set to 50; 100 and 200 produced the same error 	messages).

** Additional plink2 flags were added directly to the Step-2 QUAIL script on line 88:
*** Original: job_vqtl <- paste0(plink_path, " --bfile ", genotype, " --rm-dup 'exclude-all' 'list' ", " --pheno ", pheno_rank_score , " --covar ", covariate, " --out ", output, "_QUAIL_vQTL --linear --no-psam-pheno")
*** _JEG: job_vqtl <- paste0(plink_path, " --bfile ", genotype, " --rm-dup 'exclude-all' 'list' ",  " --pheno ", pheno_rank_score , " --covar ", covariate, " --out ", output, "_QUAIL_vQTL --linear --covar-variance-standardize --vif '300' ", " --no-psam-pheno ", " --memory '6000' ")
** The resulting file was renamed as: Step2_QUAIL_vQTL_JEG.R





#END OF README




