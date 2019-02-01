setwd("data/gbs")

library(tidyverse)


# Reformat VCF file -------------------------------------------------------
info <- read_lines("NAM_rna.vcf", n_max = 16)
vcf <- read_tsv("NAM_rna.vcf", comment = "##") %>%
  mutate(`#CHROM` = str_replace(`#CHROM`, "chr", ""), 
         ID = str_replace(ID, "chr", ""))
names(vcf) <- str_to_upper(names(vcf))
taxa <- names(vcf)[10:36]

write_lines(taxa, "NAM_to_keep.txt")
write_lines(info, "NAM_rna.vcf")
write_tsv(vcf, "NAM_rna.vcf", append = TRUE, col_names = TRUE)

rm(vcf, header, info)
gc()


# Reformat HapMap1 --------------------------------------------------------
hmp <- read_tsv("HapMap1/maizeHapMapV1_B73RefGenV2_20110309.hmp.txt", comment = "") %>%
  filter(chrom != 0)
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
  hmp[, c(1:12, idx + 11)]
})
write_tsv(hmp, "HapMap2/maizeHapMapV2.hmp.txt")
