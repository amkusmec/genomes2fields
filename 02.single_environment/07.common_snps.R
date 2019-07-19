library(tidyverse)


# Find SNPs that are segregating across all sites -------------------------
blue <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
blue <- blue[!str_detect(names(blue), "2017$")]
taxa <- lapply(blue, function(x) x$blue$PedigreeNew)

GT <- read_rds("data/gbs/add_snps.rds")$GD

seg <- lapply(taxa, function(tx) {
  which(apply(GT[rownames(GT) %in% tx, ], 2, function(x) length(unique(x)) > 1))
})
seg <- reduce(seg, intersect)
write_lines(colnames(GT)[seg], "data/single_site_gwas/common_snps.txt")
