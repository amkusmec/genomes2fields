### Standardize the reference allele across the NAM parent genotypes.

library(tidyverse)
options(readr.num_columns = 0)

setwd("data/gbs") # Make referencing files easier


for (i in 1:10) {
  cat("+ Chromosome", i, "\n")
  reference <- paste0("NAM_reference_chrom", i, ".vcf.gz") %>%
    read_tsv(., comment = "##") %>%
    select(`#CHROM`:ALT)
  target <- paste0("NAM/NAM_markers_chrom", i, ".vcf") %>%
    read_tsv(., comment = "##")
  
  # Identify the common markers
  common <- intersect(reference$ID, target$ID)
  cat("-|", length(common), "common markers\n")
  
  # Check concordance of REF/ALT alleles
  swapped <- 0
  multi <- 0
  for (j in seq_along(common)) {
    idr <- which(reference$ID == common[j])
    idt <- which(target$ID == common[j])
    
    if ((reference$REF[idr] == target$REF[idt]) & (reference$ALT[idr] != target$ALT[idt])) {
      multi <- multi + 1
    } else if ((reference$REF[idr] != target$REF[idt]) & (reference$ALT[idr] == target$ALT[idt])) {
      multi <- multi + 1
    } else if ((reference$REF[idr] == target$ALT[idt]) & (reference$ALT[idr] == target$REF[idt])) {
      swapped <- swapped + 1
      
      temp <- target$REF[idt]
      target$REF[idt] <- target$ALT[idt]
      target$ALT[idt] <- temp
      
      idx0 <- which(target[idt, ] == "0/0")
      idx1 <- which(target[idt, ] == "1/1")
      target[idt, idx0] <- "1/1"
      target[idt, idx1] <- "0/0"
    }
  }
  cat("-|", multi, "potentially multiallelic markers\n")
  cat("-| Swapped alleles for", swapped, "markers\n")
  
  # Write the swapped tables back to disk
  header <- paste0("NAM/NAM_markers_chrom", i, ".vcf") %>%
    read_lines(., n_max = 10) %>%
    write_lines(., paste0("NAM/NAM_markers_chrom", i, ".vcf"))
  
  write_tsv(target, paste0("NAM/NAM_markers_chrom", i, ".vcf"), append = TRUE, 
            col_names = TRUE)
}
