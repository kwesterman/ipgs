library(tidyverse)
library(data.table)


### Basic variables ------------------------------------------------------------

base_pheno_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz",
                       data.table = FALSE, stringsAsFactors = FALSE)

ac_fields <- c(ac = 54, ac_date = 53)
centers <- c(
  "11012" = "Barts", "11021" = "Birmingham", "11011" = "Bristol", "11008" = "Bury",
  "11003" = "Cardiff", "11024" = "Cheadle_revisit", "11020" = "Croydon",
  "11005" = "Edinburgh", "11004" = "Glasgow", "11018" = "Hounslow", "11010" = "Leeds",
  "11016" = "Liverpool", "11001" = "Manchester", "11017" = "Middlesborough",
  "11009" = "Newcastle", "11013" = "Nottingham", "11002" = "Oxford",
  "11007" = "Reading", "11014" = "Sheffield", "10003" = "Stockport_pilot",
  "11006" = "Stoke", "11022" = "Swansea", "11023" = "Wrexham",
  "11025" = "Cheadle_pilot", "11027" = "Newcastle_pilot"
)
ac_df <- base_pheno_df %>%
  select(id = f.eid, ac = f.54.0.0, ac_date = f.53.0.0) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(id = as.integer(id),
         ac = centers[ac])  # Recode numbers to location names

basic_fields <- c(sex = "f.31.0.0", age = "f.21003.0.0",  
                  bmi = "f.21001.0.0", wc = "f.48.0.0", hc = "f.49.0.0",
                  sbp = "f.4080.0.0", dbp = "f.4079.0.0", 
                  fasting_hrs = "f.74.0.0")
basic_phenos_df <- base_pheno_df %>%
  select(id = f.eid, all_of(basic_fields)) %>%
  mutate(age_squared = age^2,
         ageBySex = age * sex,
         whr = wc / hc)
# mutate(whrAdjBmi = residuals(lm(whr ~ bmi, data=., na.action=na.exclude)))

saveRDS(ac_df, "ac_df.rds")
saveRDS(basic_phenos_df, "basic_phenos_df.rds")

### Covariates and confounders -------------------------------------------------

income_fields <- c(income = "f.738.0.0")
income_coding <- c(
  "1" = "Less than 18,000",
  "2" = "18,000 to 30,999",
  "3" = "31,000 to 51,999",
  "4" = "52,000 to 100,000",
  "5" = "Greater than 100,000",
  "-1" = "Do not know",
  "-3" = "Prefer not to answer"
)
income_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_may_2023/ukb672750.tab.gz",
                   data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, all_of(income_fields)) %>%
  mutate(income = income_coding[as.character(income)])

education_fields <- c(education = "f.6138.0.0")
education_coding <- c(
  "1" = "College or University degree",
  "2" = "A levels/AS levels or equivalent",
  "3" = "O levels/GCSEs or equivalent",
  "4" = "CSEs or equivalent",
  "5" = "NVQ or HND or HNC or equivalent",
  "6" = "Other professional qualifications eg: nursing, teaching",
  "-7" = "None of the above",
  "-3" = "Prefer not to answer"
)
education_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_aug_2022/ukb669148.tab.gz",
                      data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, all_of(education_fields)) %>%
  mutate(education = education_coding[as.character(education)])

lifestyle_fields <- c(smoking = "f.20116.0.0")
smoking_coding <- c(
  "-3" = "Prefer not to answer",
  "0" =	"Never",
  "1" =	"Previous",
  "2" =	"Current"
)
lifestyle_df <- base_pheno_df %>%
  select(id = f.eid, all_of(lifestyle_fields)) %>%
  mutate(smoking = smoking_coding[as.character(smoking)])

alcohol_fields <- c(alcohol = "f.1558.0.0")
alcohol_coding <- c(
  "1" = "Daily or almost daily",
  "2" = "Three or four times a week",
  "3" = "Once or twice a week",
  "4" = "One to three times a month",
  "5" = "Special occasions only",
  "6" = "Never",
  "-3" = "Prefer not to answer"
)
alcohol_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb40167.tab.gz",
                    data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, all_of(alcohol_fields)) %>%
  mutate(alcohol = alcohol_coding[as.character(alcohol)])

covariate_df <- income_df %>%
  full_join(education_df, by = "id") %>%
  full_join(lifestyle_df, by = "id") %>%
  full_join(alcohol_df, by = "id")

saveRDS(covariate_df, "covariate_df.rds")

### Biomarkers -----------------------------------------------------------------

bm_fields <- c(
  alt = 30620, alb = 30600, apoB = 30640, hscrp = 30710, chol = 30690, glu = 30740, 
  hba1c = 30750, hdl = 30760, ldl = 30780, shbg = 30830, tg = 30870, vitD = 30890
) %>%
  sapply(function(f) paste0("f.", f, ".0.0"))

biomarker_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz", 
                      data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id = f.eid, all_of(bm_fields))

# Collect medication data for biomarker adjustments
drug_df <- base_pheno_df %>%
  select(id=f.eid, contains("f.20003.0"), contains("f.6177.0"))

statin_ids <- c(
  1140861958, 1140888594, 1140888648, 1141146234, 1141192410, 1140861922, 1141146138
)
drug_df$num_statins <- 0
for (statin in statin_ids) {
  drug_df$num_statins <- drug_df$num_statins + rowSums(drug_df[, grep("20003", names(drug_df), value=TRUE)] == statin, na.rm=TRUE)
}
drug_df$bp_med <- rowSums(drug_df[, (
  grep("6177", names(drug_df), value=TRUE)
)] == 2, na.rm=TRUE)
drug_df <- select(drug_df, id, num_statins, bp_med)

# Make medication-based biomarker adjustments
biomarker_df <- left_join(biomarker_df, 
                          select(drug_df, id, num_statins),
                          by="id")  # Not currently joining by instance!
statin_adj_bms <- c("chol", "ldl", "apoB")
statin_adj_factors <- c(
  chol = 0.749,
  ldl = 0.684,
  apoB = 0.719
)
for (bm in statin_adj_bms) {
  adj_factor <- ifelse(biomarker_df$num_statins > 0, statin_adj_factors[bm], 1)
  biomarker_df[[paste0(bm, "_statinadj")]] <- biomarker_df[[bm]] / adj_factor
}

basic_phenos_df <- left_join(basic_phenos_df, 
                             select(drug_df, id, bp_med), 
                             by = "id")
bp_adj_factors <- c(sbp = 15, dbp = 10)
for (bm in c("sbp", "dbp")) {
  adj_factor <- ifelse(basic_phenos_df$bp_med, bp_adj_factors[bm], 0)
  basic_phenos_df[[paste0(bm, "_medsadj")]] <- basic_phenos_df[[bm]] + adj_factor
}

saveRDS(biomarker_df, "biomarker_df.rds")
saveRDS(basic_phenos_df, "basic_phenos_df.rds")

### Outcomes for sample exclusion ----------------------------------------------

medical_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb38040.tab.gz",
                    data.table = FALSE, stringsAsFactors = FALSE) %>%
  select(id=f.eid, everything())
medical_df$diabetes <- rowSums(medical_df[, grepl("f\\.2443\\.0\\.", names(medical_df)), drop=FALSE] == 1, na.rm=T) > 0
medical_df$MI <- rowSums(medical_df[, grepl("f\\.6150\\.0\\.", names(medical_df)), drop=FALSE] == 1, na.rm=T) > 0
medical_df$angina <- rowSums(medical_df[, grepl("f\\.6150\\.0\\.", names(medical_df)), drop=FALSE] == 2, na.rm=T) > 0
cirrhosis_codes <- c(paste0("K70", 2:4), "K717", paste0("K74", 0:6))
cirrhosis_primary_ids <- c()
for (f in grep("f\\.41202\\.0\\.", names(medical_df), value=T)) {
  cirrhosis_primary_ids <- c(cirrhosis_primary_ids, medical_df$id[medical_df[[f]] %in% cirrhosis_codes])
}
cirrhosis_secondary_ids <- c()
for (f in grep("f\\.41204\\.0\\.", names(medical_df), value=T)) {
  cirrhosis_secondary_ids <- c(cirrhosis_secondary_ids, medical_df$id[medical_df[[f]] %in% cirrhosis_codes])
}
medical_df$pregnant <- rowSums(medical_df[, grepl("f\\.3140\\.0\\.", names(medical_df)), drop=FALSE] == 1, na.rm=T) > 0
cancer_tmp <- select(medical_df, 1, contains("f.40005."))
cancer <- cancer_tmp[, 1:7] %>%
  inner_join(select(ac_df, id, ac_date), by="id")  # Add assessment center dates
cancer$ac_date = as.Date(cancer$ac_date)
for (i in 2:7) {
  x <- ifelse(abs(difftime(cancer[, i, drop=TRUE], cancer$ac_date, units="days")) <= 365, TRUE, FALSE)  # TRUE if cancer diagnosis within a year of assessment center visit
  cancer <- cbind(cancer, x)
}
cancer$cancer_within_1yearac = apply(cancer[, 9:14], 1, function(x) {
  ifelse(any(x == TRUE, na.rm=TRUE), TRUE, FALSE)
})
cancer[names(cancer) == "x"] <- NULL

medical_df <- medical_df %>%
  left_join(select(cancer, id, cancer_within_1yearac), by="id") %>%
  mutate(CHD = MI | angina,
         cirrhosis = id %in% c(cirrhosis_primary_ids, cirrhosis_secondary_ids)) %>%
  select(id, diabetes, CHD, cirrhosis, pregnant, cancer_within_1yearac)

saveRDS(medical_df, "medical_df.rds")

### Add genetic PCs and relatedness --------------------------------------------

gPC_df <- base_pheno_df %>%
  select(id = f.eid, contains("f.22009.0"), used_in_PCA = f.22020.0.0) %>%
  rename_with(.fn = ~gsub("f.22009.0.", "gPC", .)) %>%
  mutate(unrelated = (used_in_PCA == 1)) %>%  # An unrelated subset was used in the central PCA
  select(id, all_of(paste0("gPC", 1:20)), unrelated)

saveRDS(gPC_df, "gPC_df.rds")

### Merge and write "raw" phenotypes -------------------------------------------

# ac_df <- readRDS("ac_df.rds")
# basic_phenos_df <- readRDS("basic_phenos_df.rds")
# covariate_df <- readRDS("covariate_df.rds")
# biomarker_df <- readRDS("biomarker_df.rds")
# medical_df <- readRDS("medical_df.rds")
# full_diet_df <- readRDS("full_diet_df.rds")
# gPC_df <- readRDS("gPC_df.rds")

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw27892_232_14_Nov_2022.txt", 
                          what = character())

phenos <- ac_df %>%
  left_join(basic_phenos_df, by = "id") %>%
  left_join(covariate_df, by = "id") %>%
  left_join(biomarker_df, by = "id") %>%
  left_join(medical_df, by = "id") %>%
  left_join(gPC_df, by = "id") %>%
  filter(!(id %in% withdrawn_consent)) %>%
  mutate(id = format(id, scientific = FALSE))

write_csv(phenos, "../data/processed/ukb_phenos_raw.csv")

### Phenotype processing and exclusions ----------------------------------------

logged_risk_factors <- c("alt", "tg", "hscrp")
risk_factors <- c(
  "bmi",
  "sbp", "dbp", "sbp_medsadj", "dbp_medsadj",
  "alt_log", 
  "chol", "ldl", "hdl", "apoB",
  "tg_log",
  "hba1c", "glu",
  "hscrp_log",
  "vitD",
  "chol_statinadj", "ldl_statinadj", "apoB_statinadj"
)

processed_phenos <- phenos %>%
  filter(!diabetes & !CHD & !cirrhosis & !cancer_within_1yearac & !pregnant) %>%
  mutate(across(one_of(logged_risk_factors), list(log = log))) %>%
  mutate(across(one_of(risk_factors), 
                ~ifelse(findInterval(., mean(., na.rm = TRUE) + c(-5, 5) * sd(., na.rm = TRUE)) != 1, 
                        as.numeric(NA), .)))

calc_resid_product <- function(v1, v2, covars, df) {
  adj_v1_lm_str <- paste0(v1, " ~ ", paste(covars, collapse = " + "))
  adj_v1_lm_fit <- lm(as.formula(adj_v1_lm_str), data = df, na.action = na.exclude)
  adj_v2_lm_str <- paste0(v2, " ~ ", paste(covars, collapse = " + "))
  adj_v2_lm_fit <- lm(as.formula(adj_v2_lm_str), data = df, na.action = na.exclude)
  prod_pheno <- as.vector(scale(resid(adj_v1_lm_fit))) * 
    as.vector(scale(resid(adj_v2_lm_fit)))
  prod_pheno <- ifelse(
    findInterval(prod_pheno, mean(prod_pheno, na.rm = T) + 
                   5 * c(-1, 1) * sd(prod_pheno, na.rm = T)) == 1,
    prod_pheno, NA
  )
  as.vector(scale(prod_pheno))
}
covars <- c("sex", "age", "age_squared", "ageBySex", paste0("gPC", 1:10))
processed_phenos$whr_hscrp_log_prod <- calc_resid_product("whr", "hscrp_log", covars, processed_phenos)
processed_phenos$whr_hba1c_prod <- calc_resid_product("whr", "hba1c", covars, processed_phenos)

### Write processed phenotypes -------------------------------------------------

write_csv(processed_phenos, "../data/processed/ukb_phenos.csv")

processed_phenos %>%
  filter(unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_unrelated.csv")

### Add Pan-UKBB data to generate European subset ------------------------------

anc_rel_df <- fread("/humgen/florezlab/UKBB_app27892/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_app27892.csv",
                    data.table = FALSE, stringsAsFactors = FALSE) %>%
  mutate(f.eid = as.character(f.eid),
         unrelated = !related_return2442) %>%
  select(id = f.eid, ancestry = pop_return2442, unrelated,
         one_of(paste0("PC", 1:10, "_return2442"))) %>%
  rename_at(vars(contains("PC")), ~ gsub("_return2442", "", .)) %>%
  rename_at(vars(contains("PC")), ~ gsub("PC", "gPC", .))

processed_phenos_panUKBB <- processed_phenos %>%
  select(-contains("gPC"), -unrelated) %>%
  inner_join(anc_rel_df, by = "id")

processed_phenos_panUKBB %>%
  write_csv("../data/processed/ukb_phenos_panUKBB.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "EUR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_EUR_unrelated.csv")

### Generate training-validation-test split ------------------------------------

set.seed(123)

processed_phenos_EUR <- filter(processed_phenos_panUKBB, ancestry == "EUR")

training_ids <- sample(processed_phenos_EUR$id, 
                       round(0.5 * nrow(processed_phenos_EUR)), 
                       replace = FALSE)
processed_phenos_EUR %>%
  filter(ancestry == "EUR",
         id %in% training_ids) %>%
  write_csv("../data/processed/ukb_training_set.csv")

validation_ids <- sample(setdiff(processed_phenos_EUR$id, training_ids), 
                         round(0.25 * nrow(processed_phenos_EUR)), 
                         replace = FALSE)
processed_phenos_EUR %>%
  filter(id %in% validation_ids) %>%
  write_csv("../data/processed/ukb_validation_set.csv")

processed_phenos_EUR %>%
  filter(!(id %in% c(training_ids, validation_ids))) %>%
  write_csv("../data/processed/ukb_testing_set.csv")
