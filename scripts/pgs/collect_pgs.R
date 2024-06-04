library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)
tag <- args[1]

prefix <- paste0("../data/processed/pgs/", tag)


pt_range_df <- read_table(paste0(prefix, "_pt_range_list.txt"),
                          col_names = c("threshold", "min", "max"),
                          col_types = c("cdd"))

all_chromosomes <- 1:22

sum_chromosomes <- function(pattern) {
  map(all_chromosomes, function(chr) {
    f <- str_replace(pattern, "chr[^\\.]+", paste0("chr", chr))
    read_tsv(f, show_col_types = FALSE, na = "nan") %>%
      select(id = `#IID`, score = SCORE1_SUM)
  }) %>%
    bind_rows() %>%
    group_by(id) %>%
    summarise(score = sum(score))
}

pgs_df <- map(pt_range_df$threshold, function(thresh) {
  pattern <- paste0(prefix, "_chrCHR.", thresh, ".sscore")
  score_name <- paste0("thresh", thresh)
  sum_chromosomes(pattern) %>%
    rename(!! score_name := score)
}) %>%
  reduce(~ inner_join(.x, .y, by = "id"))
write_csv(pgs_df, paste0(prefix, "_all_pgs.csv"))
