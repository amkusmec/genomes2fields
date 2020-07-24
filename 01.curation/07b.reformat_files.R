### Reformat SNP files into VCF format or changing sample names.

setwd("data/gbs")

library(tidyverse)
library(readxl)


# Find RIL families that are used in G2F ----------------------------------
families <- read_rds("../phenotype/yield_munged.rds") %>%
  pull(Pedigree) %>% str_split(., "/") %>% unlist(use.names = FALSE) %>%
  unique()
families <- families[str_detect(families, "Z[0-9]{3}")]
write_lines(families, "rils_to_keep.txt")
families <- str_replace(families, "E[0-9]{4}", "") %>% unique()

# We need data for families Z013 (KI3) and Z022 (OH43)


# Reformat VCF file -------------------------------------------------------
info <- read_lines("NAM_rna.vcf", n_max = 16)
vcf <- read_tsv("NAM_rna.vcf", comment = "##") %>%
  mutate(`#CHROM` = str_replace(`#CHROM`, "chr", ""), 
         ID = str_replace(ID, "chr", "") %>%
           str_replace(., "-", "_"))
names(vcf) <- str_to_upper(names(vcf))
taxa <- names(vcf)[10:36]
# vcf <- vcf %>% select(`#CHROM`:B73, KI3, OH43)
# taxa <- names(vcf)[10:12]

write_lines(info, "NAM_rna.vcf")
write_tsv(vcf, "NAM_rna.vcf", append = TRUE, col_names = TRUE)

rm(vcf, info)
gc()


# Reformat NAM markers ----------------------------------------------------
hmp <- map_df(1:10, function(i) {
  paste0("NAM/NAM_map_and_genos-121025/hapmap/NAM_SNP_genos_raw_20090921_chr", 
         i, ".hmp") %>%
    read_tsv(., comment = "")
}) %>%
  mutate(`rs#` = paste(chrom, pos, sep = "_"))
names(hmp)[12:4737] <- str_to_upper(names(hmp)[12:4737])
hmp <- hmp[, c(1:11, 39:4737)]
hmp <- hmp %>% select(`rs#`:QCcode, contains("Z013"), contains("Z022"))

write_tsv(hmp, "NAM/NAM_markers.hmp.txt")


# Reformat HapMap1 --------------------------------------------------------
hmp <- read_tsv("HapMap1/maizeHapMapV1_B73RefGenV2_20110309.hmp.txt", comment = "") %>%
  filter(chrom != 0, chrom != 99) %>%
  mutate(`rs#` = paste(chrom, pos, sep = "_")) %>%
  filter(pos != 256527)
taxa2 <- names(hmp)[12:38] %>%
  str_split(., ":") %>%
  map_chr(function(x) x[1])
names(hmp)[12:38] <- taxa2
# hmp <- hmp %>% select(`rs#`:B73, KI3, OH43)
write_tsv(hmp, "HapMap1/maizeHapMapV1.hmp.txt")

rm(hmp, taxa2)
gc()


# Reformat HapMap2 --------------------------------------------------------
hmp <- map_df(1:10, function(i) {
  hmp <- paste0("HapMap2/maizeHapMapV2_B73RefGenV2_201203028_chr", i, ".hmp.txt") %>%
    read_tsv(., comment = "") %>%
    filter(chrom != 0)
  taxa3 <- names(hmp)[12:116] %>% str_split(., ":") %>%
    map_chr(function(x) x[1])
  names(hmp)[12:116] <- taxa3
  idx <- which(taxa3 %in% taxa)
  hmp[, c(1:11, idx + 11)]
}) %>%
  mutate(`rs#` = paste(chrom, pos, sep = "_"))
write_tsv(hmp, "HapMap2/maizeHapMapV2.hmp.txt")
