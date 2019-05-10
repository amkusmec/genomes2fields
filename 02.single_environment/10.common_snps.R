library(tidyverse)


# Find SNPs that are segregating across all sites -------------------------
blue <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
sites <- names(blue)
sites <- sites[!str_detect(sites, "2017$")]

assoc <- seq_along(sites) %>%
  map(function(s) {
    paste0("data/gemma/output/single", s, ".assoc.txt") %>%
      read_tsv() %>%
      pull(rs)
  })
common <- reduce(assoc, intersect)
write_lines(common, "data/gemma/common_snps.txt")
