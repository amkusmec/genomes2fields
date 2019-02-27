library(tidyverse)
library(parallel)


# Get the combined GBS SNP file -------------------------------------------
setwd("data/gbs")
system(paste0("cp /media/analyses/SRP106367+SRP021921.Combined.Aaron.SNPs/", 
              "SNPs/SRP106367+SRP021921.Combined.Aaron.SNPs_Matrix-20190222.", 
              "Recalled.genotypes.vcf.gz ./G2F_Ames_combined.vcf.gz"))
system("gunzip G2F_Ames_combined.vcf.gz")
setwd("../..")


# Load and parse the SNPs -------------------------------------------------
header <- read_lines("data/gbs/G2F_Ames_combined.vcf", n_max = 16)
vcf <- read_tsv("data/gbs/G2F_Ames_combined.vcf", comment = "##") %>%
  filter(!(`#CHROM` %in% c("chrMt", "chrPt", "chrUNKNOWN"))) %>%
  mutate(`#CHROM` = str_replace(`#CHROM`, "chr", "") %>% as.integer(), 
         ID = str_replace(ID, "chr", "") %>% str_replace(., "-", "_"))

snp_info <- vcf[, 1:9]
vcf <- vcf[, -c(1:9)]

cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(tidyverse))

split_idx <- clusterSplit(cl, 1:nrow(vcf))
split_idx <- rep(1:detectCores(), times = sapply(split_idx, length))
vcf <- split(vcf, split_idx)

parsed <- parLapply(cl, vcf, function(v) {
  apply(v, 2, function(x) {
    str_split(x, ":") %>%
      sapply(., function(y) y[1])
  })
})

stopCluster(cl)
parsed <- do.call("rbind", parsed)
rm(cl, vcf, split_idx)
gc()


# Identify and handle duplicates ------------------------------------------
parsed <- as_tibble(parsed) %>%
  bind_cols(snp_info, .) %>%
  arrange(`#CHROM`, POS)

# Only keep "true" duplicates; i.e., some samples have "DUP" in the name but
# are the only remaining sample for that taxon
names(parsed)[str_detect(names(parsed), "TR")] <- paste0("TR9_1_1_6_DUP", 1:2)
duplicates <- which(str_detect(colnames(parsed), "DUP"))
dup_table <- tibble(Taxa = str_replace(colnames(parsed)[duplicates], 
                                       "_DUP[0-9]{1,2}", ""),
                    Label = colnames(parsed)[duplicates], 
                    Column = duplicates) %>%
  filter(duplicated(Taxa) | duplicated(Taxa, fromLast = TRUE))

# Across all duplicates, set discordant calls to missing
combined <- dup_table %>%
  split(., .$Taxa) %>%
  map(function(df) {
    samples <- parsed[, df$Column]
    consensus <- apply(samples, 1, function(x) {
      calls <- unique(x[x != "./."])
      if (length(calls) == 0) {
        "./."
      } else if (length(calls) == 1) {
        calls
      } else if (length(calls) > 1) {
        if (any(str_detect(calls, "0/1") & str_detect(calls, "1/0"))) {
          "0/1"
        } else {
          "./."
        }
      }
    })
  }) %>%
  bind_cols()

# Combine consensus duplicates with single samples and remove the "DUP" identifiers
parsed <- parsed[, -dup_table$Column] %>%
  bind_cols(., combined)
names(parsed) <- str_replace(names(parsed), "_DUP[0-9]{1,2}", "")


# Split by chromosome for easier imputation -------------------------------
parsed <- split(parsed, parsed$`#CHROM`)
map(parsed, function(df) {
  write_lines(header, paste0("data/gbs/G2F_Ames_combined_chrom", unique(df$`#CHROM`), ".vcf"))
  write_tsv(df, paste0("data/gbs/G2F_Ames_combined_chrom", unique(df$`#CHROM`), ".vcf"), 
            col_names = TRUE, append = TRUE)
  NULL
})
