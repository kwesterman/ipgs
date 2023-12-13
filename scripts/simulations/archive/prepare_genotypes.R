library(tidyverse)


N <- 10000
M <- 100
min_maf <- 0.01
max_maf <- 0.5

set.seed(1)

mafs <- runif(M, min_maf, max_maf)

g_mat <- sapply(mafs, function(maf) {
  rbinom(N, 2, maf)
}, simplify = TRUE)
colnames(g_mat) <- paste0("g", seq(1, M))

saveRDS(g_mat, "../data/processed/simulations/g_mat.rds")