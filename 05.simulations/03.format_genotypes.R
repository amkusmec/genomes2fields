library(tidyverse)
library(lubridate)
library(purrrlyr)
library(parallel)
library(readxl)


# Remove plots with extreme flowering times -------------------------------
weather <- read_rds("data/weather/env_variables.rds") %>%
  mutate(Date = paste(Year, Month, Day, sep = "-") %>% ymd())
yield <- read_rds("data/phenotype/yield_munged.rds") %>%
  by_row(function(r) {
    a <- weather %>%
      filter(Environment == r$Environment, Date >= r$Planted, 
             Date <= r$Anthesis) %>%
      summarise(A = sum(CHU, na.rm = TRUE)) %>%
      unlist()
    
    tibble(CHUA = a)
  }, .collate = "rows") %>%
  select(-.row) %>%
  mutate(Site = paste(Environment, Year, sep = "_"))

cutoffs <- yield %>%
  group_by(Site) %>%
  summarise(Lower = quantile(CHUA, 0.025), 
            Upper = quantile(CHUA, 0.975)) %>%
  ungroup()

sites <- weather %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  pull(Site) %>%
  unique()
yield <- inner_join(yield, cutoffs, by = "Site") %>%
  filter(Site %in% sites) %>%
  filter(CHUA >= Lower, CHUA <= Upper) %>%
  select(-Lower, -Upper)


# Identify inbred parents -------------------------------------------------
yield <- yield %>%
  mutate(Pedigree = str_replace(Pedigree, "\\(1\\)", "_1") %>%
           str_replace(., "DK3iiH6", "3IIH6") %>%
           str_replace(., "ICI 441", "ICI441") %>%
           str_replace(., "ILLHY", "ILL_HY") %>%
           str_replace(., "TR 9-1-1-6", "TR9_1_1_6") %>%
           str_replace(., "\\?\\?\\?TX205", "TX205") %>%
           str_replace(., "MBNIL B", "MBNILB") %>%
           str_replace(., "LH195EXPVP", "LH195") %>%
           str_replace(., "LH195EX-PVP", "LH195"),
         Female = str_split(Pedigree, "/") %>% sapply(., function(x) x[1]) %>%
           str_to_upper() %>%
           str_replace_all(., "\\.", "") %>% str_replace_all(., "-", "_"), 
         Male = str_split(Pedigree, "/") %>% sapply(., function(x) x[2]) %>%
           str_to_upper() %>%
           str_replace_all(., "\\.", "") %>% str_replace_all(., "-", "_"), 
         PedigreeNew = paste(Female, Male, sep = "/"))

# Change some reciprocal hybrids and remove one self cross
yield <- yield %>%
  mutate(PedigreeNew2 = paste(Male, Female, sep = "/"), 
         Test = PedigreeNew2 %in% intersect(PedigreeNew, PedigreeNew2),
         PedigreeNew = if_else(Test, 
                               if_else(Male >= Female, paste(Male, Female, collapse = "/"), 
                                       paste(Female, Male, collapse = "/")), 
                               PedigreeNew), 
         Female = str_split(PedigreeNew, "/") %>% sapply(., function(x) x[1]),
         Male = str_split(PedigreeNew, "/") %>% sapply(., function(x) x[2])) %>%
  select(-PedigreeNew2, -Test) %>%
  filter(PedigreeNew != "PHB47/PHB47")

# Keep the inbred parents
inbreds <- c(pull(yield, Female), pull(yield, Male)) %>%
  unique() %>% sort()


# Load the GBS files ------------------------------------------------------
gbs <- list.files("data/gbs", "GAN_reference_chrom[0-9]{1,2}.vcf", 
                  full.names = TRUE) %>%
  map_df(read_tsv, comment = "##") %>%
  arrange(`#CHROM`, POS)
taxa <- names(gbs)[-c(1:9)] %>%
  str_replace(., "MBNIL_B", "MBNILB")

meta <- read_xlsx("data/gbs/AmesUSInbreds_AllZeaGBSv1.0_SampleNameKey.xlsx") %>% 
  mutate(DNA = `Inbred name` %>% str_to_upper() %>%
           str_replace_all(., "\\.", "") %>% str_replace_all(., "_", ""),
         GBS = `GBS name (Sample:Flowcell:Lane:Well)` %>%
           str_split(., ":") %>% sapply(., function(x) x[1]) %>%
           str_to_upper() %>% str_replace_all(., "\\.", "") %>% 
           str_replace_all(., "_", "")) %>%
  select(DNA, GBS) %>%
  filter(GBS %in% taxa & (str_detect(GBS, "AMES") | str_detect(GBS, "PI"))) %>%
  distinct(GBS, .keep_all = TRUE) %>%
  filter(!(DNA %in% taxa)) %>%
  filter(!(DNA %in% c("ILL HY", "TZI 8")))
idxt <- which(taxa %in% meta$GBS)
taxa[idxt] <- meta$DNA

# Keep SNPs used in the hybrid panel
snps <- read_rds("data/gbs/add_snps.rds")
GM <- snps$GM

gbs <- gbs %>%
  mutate(ID = paste0("X", ID)) %>%
  filter(ID %in% GM$SNP)


# Create synthetic hybrid genotypes ---------------------------------------
hybrids <- yield %>%
  select(PedigreeNew, Female, Male) %>%
  distinct(PedigreeNew, .keep_all = TRUE) %>%
  filter(Female %in% taxa, Male %in% taxa) %>%
  mutate(FemaleRow = match(Female, taxa), 
         MaleRow = match(Male, taxa))

gbs1 <- gbs[, -c(1:9)] %>%
  apply(., 2, str_remove, pattern = "\\|[012]")
gbs2 <- gbs[, -c(1:9)] %>%
  apply(., 2, str_remove, pattern = "[012]\\|")
gbs_list <- list(gbs1, gbs2)

cl <- makeCluster(10)
clusterExport(cl, list("gbs_list"))

hyb <- clusterMap(cl, function(x, y) {
  g <- sample(2, 2, replace = TRUE)
  paste(gbs_list[[g[1]]][, x], gbs_list[[g[2]]][, y], sep = "|")
}, x = hybrids$FemaleRow, y = hybrids$MaleRow, 
RECYCLE = FALSE, .scheduling = "static")
hyb <- t(matrix(unlist(hyb, use.names = FALSE), ncol = length(hyb)))
dimnames(hyb) <- list(hybrids$PedigreeNew, GM$SNP)

stopCluster(cl)
rm(cl, gbs, gbs_list, gbs1, gbs2)
gc()

hyb <- hyb[rownames(hyb) %in% rownames(snps$GD), ]

write_rds(hyb, "data/gbs/synthetic_hybrids_phased.rds")


# Reformat haplotypes for simulation --------------------------------------
haplotypes <- list()
for (i in 1:10) {
  haplotypes[[i]] <- hyb[, which(map$Chr == i)] %>%
    apply(., 2, function(x) {
      temp <- str_split(x, "\\|") %>% unlist(use.names = FALSE) %>% as.numeric()
      if_else(temp > 1, 1, temp) %>% 
        matrix(., ncol = 1, nrow = length(x)*2)
    })
}
write_rds(haplotypes, "data/gbs/synthetic_hybrids_haplotypes.rds")
