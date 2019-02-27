library(tidyverse)
library(parallel)
options(readr.num_columns = 0)

setwd("data/gbs") # Make referencing files easier
cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(tidyverse))


for (i in 2:10) {
  cat("+ Chromosome", i, "\n")
  system(paste0("gunzip NAM_imputed_chrom", i, ".vcf.gz"))
  reference <- paste0("NAM_imputed_chrom", i, ".vcf") %>%
    read_tsv(., comment = "##")
  target <- paste0("G2F_Ames_combined_chrom", i, ".vcf") %>%
    read_tsv(., comment = "##")
  
  # Identify the common markers
  common <- intersect(reference$ID, target$ID)
  cat("-|", length(common), "common markers\n")
  
  # Subsets for parallel processing
  ref_temp <- filter(reference, ID %in% common)
  tar_temp <- filter(target, ID %in% common)
  
  split_idx <- clusterSplit(cl, seq_along(common))
  split_idx <- rep(1:detectCores(), times = sapply(split_idx, length))
  
  ref_temp <- split(ref_temp, split_idx)
  tar_temp <- split(tar_temp, split_idx)
  
  # Check concordance of REF/ALT alleles
  res <- clusterMap(cl, function(reference, target) {
    concord <- 0
    swapped <- 0
    multi <- 0
    remove <- c()
    
    for (j in 1:nrow(reference)) {
      if ((reference$REF[j] == target$REF[j]) & 
          (reference$ALT[j] != target$ALT[j])) {
        multi <- multi + 1
      } else if (reference$REF[j] != target$REF[j]) {
        # Check for the presence of reference$REF in target$ALT
        talt <- unlist(str_split(target$ALT[j], ","))
        ralt <- unlist(str_split(reference$ALT[j], ","))
        if (reference$REF[j] %in% talt) { # Allele is present; swap target
          swapped <- swapped + 1
          
          t_idx <- which(talt == reference$REF[j])
          temp <- target$REF[j]
          target$REF[j] <- talt[t_idx]
          talt[t_idx] <- temp
          target$ALT[j] <- paste(talt, collapse = ",")
          
          temp <- as.character(length(talt) + 1)
          target[j, 10:ncol(target)] <- str_replace(target[j, 10:ncol(target)], 
                                                    "0", temp)
          target[j, 10:ncol(target)] <- str_replace(target[j, 10:ncol(target)], 
                                                    as.character(t_idx), "0")
          target[j, 10:ncol(target)] <- str_replace(target[j, 10:ncol(target)], 
                                                    temp, as.character(t_idx))
        } else if (any(ralt %in% talt)) { # Matching alternate alleles; swap both
          swapped <- swapped + 1
          
          r_idx <- which(ralt %in% talt)[1]
          t_idx <- which(talt == ralt[r_idx])
          
          # Swap reference alleles
          temp <- reference$REF[j]
          reference$REF[j] <- ralt[r_idx]
          ralt[r_idx] <- temp
          reference$ALT[j] <- paste(ralt, collapse = ",")
          
          temp <- as.character(length(ralt) + 1)
          reference[j, 10:ncol(reference)] <- str_replace(reference[j, 10:ncol(reference)], 
                                                          "0", temp)
          reference[j, 10:ncol(reference)] <- str_replace(reference[j, 10:ncol(reference)], 
                                                          as.character(r_idx), "0")
          reference[j, 10:ncol(reference)] <- str_replace(reference[j, 10:ncol(reference)], 
                                                          temp, as.character(r_idx))
          
          # Swap target alleles
          temp <- target$REF[j]
          target$REF[j] <- talt[t_idx]
          talt[t_idx] <- temp
          target$ALT[j] <- paste(talt, collapse = ",")
          
          temp <- as.character(length(talt) + 1)
          target[j, 10:ncol(target)] <- str_replace(target[j, 10:ncol(target)], 
                                                    "0", temp)
          target[j, 10:ncol(target)] <- str_replace(target[j, 10:ncol(target)], 
                                                    as.character(t_idx), "0")
          target[j, 10:ncol(target)] <- str_replace(target[j, 10:ncol(target)], 
                                                    temp, as.character(t_idx))
        } else { # No matching alleles; remove site
          remove <- c(remove, reference$ID[j])
        }
      } else { # Alleles are concordant 
        concord <- concord + 1
      }
    }
    
    return(list(reference = reference, target = target, swapped = swapped, 
                multi = multi, concord = concord, remove = remove))
  }, reference = ref_temp, target = tar_temp)
  
  cat("-|", sum(sapply(res, function(x) x$concord)), "concordant sites\n")
  cat("-|", sum(sapply(res, function(x) x$multi)), "multiallelic sites\n")
  cat("-| Swapped alleles for", sum(sapply(res, function(x) x$swapped)), "sites\n")
  cat("-|", sum(sapply(res, function(x) length(x$remove))), "sites to remove\n")
  
  # Replace common sites
  ref_common <- do.call("rbind", lapply(res, function(x) x$reference)) %>%
    as_tibble()
  tar_common <- do.call("rbind", lapply(res, function(x) x$target)) %>%
    as_tibble()
  
  idx <- which(reference$ID %in% ref_common$ID)
  reference[idx, ] <- ref_common
  
  idx <- which(target$ID %in% tar_common$ID)
  target[idx, ] <- tar_common
  
  # Remove sites with no matching alleles
  remove <- sapply(res, function(x) x$remove) %>% unlist(use.names = FALSE)
  reference <- filter(reference, !(ID %in% remove))
  target <- filter(target, !(ID %in% remove))
  
  ### Write the swapped tables back to disk
  # G2F/Ames
  paste0("G2F_Ames_combined_chrom", i, ".vcf") %>%
    read_lines(., n_max = 16) %>%
    write_lines(., paste0("G2F_Ames_combined_chrom", i, ".vcf"))
  write_tsv(target, paste0("G2F_Ames_combined_chrom", i, ".vcf"), append = TRUE, 
            col_names = TRUE)
  
  # NAM
  paste0("NAM_imputed_chrom", i, ".vcf") %>%
    read_lines(., n_max = 9) %>%
    write_lines(., paste0("NAM_imputed_chrom", i, ".vcf"))
  write_tsv(reference, paste0("NAM_imputed_chrom", i, ".vcf"), append = TRUE, 
            col_names = TRUE)
}

stopCluster(cl)
