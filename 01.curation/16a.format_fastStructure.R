library(tidyverse)


snps <- read_rds("data/gbs/add_snps.rds")
header <- read_lines("data/gbs/GAN_reference_chrom1.vcf", n_max = 9)
info <- snps$GM %>%
  select(Chromosome, Position, SNP, Alleles) %>%
  separate(Alleles, c("REF", "ALT"), sep = ",", remove = TRUE) %>%
  mutate(QUAL = ".", 
         FILTER = "PASS", 
         INFO = ".", 
         FORMAT = "GT") %>%
  rename(`#CHROM` = Chromosome, POS = Position, ID = SNP)

GD <- t(snps$GD) %>%
  apply(., 2, function(x) { 
    if_else(x == 0, "0|0", 
            if_else(x == 1, "0|1", "1|1"))
  })

vcf <- bind_cols(info, as_tibble(GD))

write_lines(header, "data/gbs/g2f_hyb.vcf")
write_tsv(vcf, "data/gbs/g2f_hyb.vcf", append = TRUE, col_names = TRUE)

system("src/plink --vcf data/gbs/g2f_hyb.vcf --double-id --out data/gbs/g2f_hyb_bin")


write_lines(c("#!/bin/bash", ""), "01.curation/16b.fastStructure.sh")
for (i in 2:15) {
  for (j in 1:19) {
    paste0("python ~/trs_selection/src/fastStructure/structure.py -K ", 
           i, " --input=data/gbs/g2f_hyb_bin --output=data/structure/g2f_hyb_rep", j, " &") %>%
      write_lines(., "01.curation/16b.fastStructure.sh", append = TRUE)
  }
  
  paste0("python ~/trs_selection/src/fastStructure/structure.py -K ", 
         i, " --input=data/gbs/g2f_hyb_bin --output=data/structure/g2f_hyb_rep20") %>%
    write_lines(., "01.curation/16b.fastStructure.sh", append = TRUE)
}
