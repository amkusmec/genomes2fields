### Examines the NAM RIL GBS data from the NCBI SRA. All samples have too much 
### missing data to use.

library(tidyverse)
library(parallel)


# Move and uncompress SNPs ------------------------------------------------
# system(paste0("cp /media/analyses/SRP009896.GBS.Sequencing.of.NAM/SNPs/", 
#               "SRP009896.GBS.Sequencing.of.NAM.SNPs_Matrix.Recalled.", 
#               "genotypes.vcf.gz data/gbs/recalled_NAM_gbs.vcf.gz"))
# system("gunzip data/gbs/recalled_NAM_gbs.vcf.gz")

system(paste0("cp /media/analyses/SRP009896.GBS.Sequencing.of.NAM/", 
              "SNPs.at.least.3.reads/SRP009896.GBS.Sequencing.of.NAM.", 
              "SNPs_Matrix_N3.Recalled.genotypes.vcf.gz data/gbs/recalled_NAM3", 
              "_gbs.vcf.gz"))
system("gunzip data/gbs/recalled_NAM3_gbs.vcf.gz")


# Split the information from the SNP calls --------------------------------
vcf <- read_tsv("data/gbs/recalled_NAM3_gbs.vcf", comment = "##") %>%
  filter(!(`#CHROM` %in% c("chrMt", "chrPt", "chrUNKNOWN")))

snp_info <- vcf[, 1:9] %>%
  rename(CHROM = `#CHROM`) %>%
  mutate(CHROM = str_replace(CHROM, "chr", "") %>% as.integer(), 
         ID = str_replace(ID, "chr", "") %>%
           str_replace(., "-", "_") %>%
           paste0("X", .)) %>%
  select(ID, CHROM, POS, REF, ALT)

vcf <- vcf[, -c(1:9)]


# Keep only samples that match inbreds in the yield files -----------------
inbreds <- read_rds("data/phenotype/yield_munged.rds") %>%
  pull(Pedigree) %>%
  str_split(., "/") %>%
  unlist(use.names = FALSE) %>%
  unique() %>%
  sort() %>%
  str_to_upper() %>%
  make.names() %>%
  str_replace_all(., "\\.", "") %>%
  str_replace_all(., "_", "")

vcf_names <- names(vcf) %>%
  str_replace(., "_DUP[0-9]{1,2}", "") %>%
  str_to_upper() %>% make.names() %>%
  str_replace_all(., "\\.", "") %>%
  str_replace_all(., "_", "")

idx_vcf <- which(vcf_names %in% inbreds)

# Exclude inbreds already in the G2F and Ames files
g2f_lines <- read_lines("data/gbs/G2F_gbs_keep.txt") %>%
  str_split(., ",") %>% unlist(use.names = FALSE) %>%
  str_replace(., "_DUP[0-9]{1,2}", "") %>%
  unique() %>% sort() %>% make.names()
ames_lines <- read_lines("data/gbs/Ames_gbs_keep.txt") %>%
  str_split(., ",") %>% unlist(use.names = FALSE) %>%
  str_replace(., "_DUP[0-9]{1,2}", "") %>%
  unique() %>% sort() %>% make.names()

idx_g2f <- which(vcf_names %in% g2f_lines)
idx_ames <- which(vcf_names %in% ames_lines)
idx_final <- setdiff(idx_vcf, idx_g2f)

vcf <- vcf[, idx_final]


# Parse the VCF file ------------------------------------------------------
cl <- makeCluster(10)
clusterEvalQ(cl, library(tidyverse))
vcf <- split(vcf, snp_info$CHROM)

parsed <- parLapply(cl, vcf, function(v) {
  apply(v, 2, function(x) {
    str_split(x, ":") %>%
      sapply(., function(y) y[1])
  })
})

stopCluster(cl)
parsed <- do.call("rbind", parsed)
write_rds(parsed, "data/gbs/recalled_NAM3_gbs_parsed.rds")
write_rds(snp_info, "data/gbs/recalled_NAM3_gbs_snp_info.rds")

rm(cl, vcf)
gc()


# Convert parsed VCF into ABH calls ---------------------------------------
parsed <- apply(parsed, 2, function(x) {
  if_else(x == "./.", "N", 
          if_else(x == "0/0", "A", 
                  if_else(x == "1/1", "B", 
                          if_else(x == "0/1" | x == "1/0", "H", as.character(NA)))))
})


# Require at least 20% non-missing calls per sample -----------------------
missing <- apply(parsed, 2, function(x) sum(x == "N")/length(x))
tibble(Missing = missing) %>%
  ggplot(., aes(x = Missing)) + theme_classic() +
    geom_histogram(binwidth = 0.01, colour = "black", fill = "red", alpha = 0.8) +
    geom_vline(xintercept = 0.8, linetype = 2) +
    labs(x = "% Missing Calls", y = "Count") +
    scale_x_continuous(labels = scales::percent) +
    ggtitle("NAM GBS (3 reads)")
ggsave("figures/munge/NAM3_gbs_missing.pdf", width = 6, height = 4, units = "in", 
       dpi = 300)
parsed <- parsed[, missing <= 0.8]

### No samples have less than 80% missing calls.
