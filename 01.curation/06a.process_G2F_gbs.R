library(tidyverse)
library(parallel)


# Get G2F SNP metadata ----------------------------------------------------
### See 01.munge_metadata.R for information about configuring icommands.

### This line will need to change based on where you installed icommands.
Sys.setenv("PATH" = paste("/mnt/01/amkusmec/bin/icommands", Sys.getenv("PATH"), sep = ":"))

setwd("data/gbs")
system("iinit")
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Nov_2016_V.3/b._2014_gbs_data/", 
              "g2f_2014_gbs_data.csv"))
system("iexit")
setwd("../..")


# Move and uncompress SNPs ------------------------------------------------
system(paste0("cp /media/analyses/SRP106367.PRJNA385022.GBS.Sequencing.of.", 
              "Maize.G2F.GxE.Inbreds/SNPs/SRP106367.PRJNA385022.", 
              "SNPs_Matrix-20181202.Recalled.genotypes.vcf.gz ", 
              "data/gbs/recalled_G2F_gbs.vcf.gz"))
system("gunzip data/gbs/recalled_G2F_gbs.vcf.gz")


# Split the information from the SNP calls --------------------------------
vcf <- read_tsv("data/gbs/recalled_G2F_gbs.vcf", comment = "##") %>%
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

meta <- read_csv("data/gbs/g2f_2014_gbs_data.csv") %>%
  filter(`<DNAName>` != "BLANK") %>%
  mutate(DNA = `<DNAName>` %>% str_to_upper() %>%
           make.names() %>% str_replace_all(., "\\.", "") %>%
           str_replace_all(., "_", ""), 
         GBS = `<GBSSampleName>` %>% str_split(., ":") %>%
           sapply(., function(x) x[1]) %>% str_to_upper() %>%
           make.names() %>% str_replace_all(., "\\.", "") %>%
           str_replace_all(., "_", "")) %>%
  select(DNA, GBS)

vcf_names <- names(vcf) %>%
  str_replace(., "_DUP[0-9]{1,2}", "") %>%
  str_to_upper() %>% make.names() %>%
  str_replace_all(., "\\.", "") %>%
  str_replace_all(., "_", "")

idx_inbred <- which(meta$DNA %in% inbreds | meta$GBS %in% inbreds)
idx_vcf <- which(vcf_names %in% unlist(meta[idx_inbred, ], use.names = FALSE))

vcf <- vcf[, idx_vcf]


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
write_rds(parsed, "data/gbs/recalled_G2F_gbs_parsed.rds")
write_rds(snp_info, "data/gbs/recalled_G2F_gbs_snp_info.rds")

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
    ggtitle("G2F GBS")
ggsave("figures/munge/G2F_gbs_missing.pdf", width = 6, height = 4, units = "in", 
       dpi = 300)
parsed <- parsed[, missing <= 0.8]


# Identify duplicate samples ----------------------------------------------
duplicates <- which(str_detect(colnames(parsed), "DUP"))
dup_table <- tibble(Taxa = str_replace(colnames(parsed)[duplicates], 
                                       "_DUP[0-9]{1,2}", ""), 
                    Label = colnames(parsed)[duplicates], 
                    Column = duplicates)


# Compute IBS between all pairs of duplicates -----------------------------
ibs <- dup_table %>%
  split(., .$Taxa) %>%
  map_df(function(df) {
    if (nrow(df) == 1) {
      return(tibble(Taxa = NA, Col1 = NA, Col2 = NA, IBS = NA))
    }
    
    combos <- combn(1:nrow(df), 2)
    temp <- tibble(Taxa = df$Taxa[1], 
                   Col1 = df$Label[as.vector(combos[1, ])], 
                   Col2 = df$Label[as.vector(combos[2, ])], 
                   IBS = 0)
    for (i in 1:ncol(combos)) {
      use <- parsed[, temp$Col1[i]] != "N" & parsed[, temp$Col2[i]] != "N"
      temp$IBS[i] <- sum(parsed[use, temp$Col1[i]] == parsed[use, temp$Col2[i]])/sum(use)
    }
    
    temp
  }) %>%
  filter(!is.na(Col1))

ggplot(ibs, aes(x = IBS)) + theme_classic() +
  geom_histogram(binwidth = 0.01, colour = "black", fill = "blue", alpha = 0.8) +
  geom_vline(xintercept = 0.99, linetype = 2) +
  labs(x = "Identitiy by State", y = "Count") +
  ggtitle("G2F GBS")
ggsave("figures/munge/G2F_gbs_ibs.pdf", width = 6, height = 4, units = "in", 
       dpi = 300)


# Identify samples that can be merged -------------------------------------
to_merge <- ibs %>%
  filter(IBS >= 0.99) %>%
  group_by(Taxa) %>%
  summarise(Columns = c(Col1, Col2) %>% unique() %>%
              sort() %>% paste(., collapse = ","))


# Samples to call but not merge -------------------------------------------
no_merge <- tibble(Columns = colnames(parsed)) %>%
  mutate(Taxa = str_replace(Columns, "_DUP[0-9]{1,2}", "")) %>%
  filter(!(Taxa %in% to_merge$Taxa)) %>%
  select(Taxa, Columns)

bind_rows(to_merge, no_merge) %>%
  pull(Columns) %>%
  write_lines(., "data/gbs/G2F_gbs_keep.txt")
