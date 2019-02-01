setwd("data/gbs")

library(tidyverse)
library(readxl)


# Reformat VCF file -------------------------------------------------------
info <- read_lines("NAM_rna.vcf", n_max = 16)
vcf <- read_tsv("NAM_rna.vcf", comment = "##") %>%
  mutate(`#CHROM` = str_replace(`#CHROM`, "chr", ""), 
         ID = str_replace(ID, "chr", "") %>%
           str_replace(., "-", "_"))
names(vcf) <- str_to_upper(names(vcf))
taxa <- names(vcf)[10:36]

write_lines(taxa, "NAM_to_keep.txt")
write_lines(info, "NAM_rna.vcf")
write_tsv(vcf, "NAM_rna.vcf", append = TRUE, col_names = TRUE)

rm(vcf, info)
gc()


# Reformat NAM markers ----------------------------------------------------
nam_map <- read_xls("NAM/NAM_map_and_genos-121025/NAM_map_20080419.xls") %>%
  select(marker, ch, cumulative)
nam_pos <- read_xlsx("NAM/NAM_map_and_genos-121025/NAM_1144SNPs_AGPv2_positions.xlsx") %>%
  select(SNP_NAME:AGPv2_pos) %>%
  mutate(SNP_NAME = as.character(SNP_NAME),
         AGPv2_Chr = as.numeric(as.character(AGPv2_Chr)),
         AGPv2_pos = as.numeric(as.character(AGPv2_pos)))
nam_map <- inner_join(nam_map, nam_pos, by = c("marker" = "SNP_NAME")) %>%
  select(marker, AGPv2_Chr, AGPv2_pos, cumulative) %>%
  rename(Marker = marker, Chr = AGPv2_Chr, Physical = AGPv2_pos, 
         Genetic = cumulative) %>%
  filter(!is.na(Physical)) %>%
  arrange(Chr, Physical)
nam_map <- nam_map[-c(243, 279, 280, 319, 450, 460, 503, 518, 533, 698, 702, 746, 
                      748, 838, 887, 961), ]
nam_map <- nam_map[-873, ]
nam_map <- select(nam_map, Chr, Marker, Genetic, Physical)

hmp <- map_df(1:10, function(i) {
  paste0("NAM/NAM_map_and_genos-121025/hapmap/NAM_SNP_genos_raw_20090921_chr", 
         i, ".hmp") %>%
    read_tsv(., comment = "") %>%
    filter(`rs#` %in% nam_map$Marker)
}) %>%
  mutate(`rs#` = paste(chrom, pos, sep = "_"))
names(hmp)[12:4737] <- str_to_upper(names(hmp)[12:4737])

nam_map <- mutate(nam_map, Marker = paste(Chr, Physical, sep = "_"))

write_tsv(hmp, "NAM/NAM_markers.hmp.txt")
write_tsv(nam_map, "NAM/NAM_markers.map", col_names = FALSE)


# Reformat HapMap1 --------------------------------------------------------
hmp <- read_tsv("HapMap1/maizeHapMapV1_B73RefGenV2_20110309.hmp.txt", comment = "") %>%
  filter(chrom != 0, chrom != 99) %>%
  mutate(`rs#` = paste(chrom, pos, sep = "_"))
taxa2 <- names(hmp)[12:38] %>%
  str_split(., ":") %>%
  map_chr(function(x) x[1])
names(hmp)[12:38] <- taxa2
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
