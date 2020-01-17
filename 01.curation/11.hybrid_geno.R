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
                               if_else(Male >= Female, paste(Male, Female, sep = "/"), 
                                       paste(Female, Male, sep = "/")), 
                               PedigreeNew)) %>%
  separate(PedigreeNew, c("Male", "Female"), sep = "/", remove = FALSE) %>%
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

# Keep SNP information
GM <- gbs[, 1:9] %>%
  rename(Chromosome = `#CHROM`, Position = POS, SNP = ID) %>%
  mutate(Alleles = paste(REF, ALT, sep = ","), 
         SNP = make.names(SNP)) %>%
  select(SNP, Chromosome, Position, Alleles)

# Swap dimensions so that samples are rows and SNPs are columns
gbs <- t(gbs[, -c(1:9)])
colnames(gbs) <- GM$SNP
rownames(gbs) <- taxa

# Convert to numeric genotypes
gbs <- apply(gbs, 2, function(x) {
  if_else(x == "0|0", 0, 
          if_else(x == "0|1" | x == "1|0", 1, 2))
})
rownames(gbs) <- taxa


# Keep only SNPs that are segregating in both panels ----------------------
idx <- which(str_detect(rownames(gbs), "Z[0-9]{3}E[0-9]{4}"))

seg_nam <- apply(gbs[idx, ], 2, function(x) length(unique(x)))
seg_nam <- names(seg_nam)[seg_nam > 1]

seg_g2f <- apply(gbs[-idx, ], 2, function(x) length(unique(x)))
seg_g2f <- names(seg_g2f)[seg_g2f > 1]

seg <- intersect(seg_nam, seg_g2f)
gbs <- gbs[, colnames(gbs) %in% seg]
GM <- filter(GM, SNP %in% seg)


# Create synthetic hybrid genotypes ---------------------------------------
hybrids <- yield %>%
  select(PedigreeNew, Female, Male) %>%
  distinct(PedigreeNew, .keep_all = TRUE) %>%
  filter(Female %in% taxa, Male %in% taxa) %>%
  mutate(FemaleRow = match(Female, taxa), 
         MaleRow = match(Male, taxa))

cl <- makeCluster(10)
clusterExport(cl, list("gbs"))

hyb <- clusterMap(cl, function(x, y) {
  apply(gbs[c(x, y), ], 2, mean, na.rm = FALSE)
}, x = hybrids$FemaleRow, y = hybrids$MaleRow, 
  RECYCLE = FALSE, .scheduling = "static")
hyb <- t(matrix(unlist(hyb, use.names = FALSE), ncol = length(hyb)))
dimnames(hyb) <- list(hybrids$PedigreeNew, GM$SNP)

stopCluster(cl)
rm(cl, gbs)
gc()

write_rds(hyb, "data/gbs/synthetic_hybrids.rds")


# Split hybrid genotypes --------------------------------------------------
hyb1416 <- yield %>% 
  filter(Year <= 2016) %>%
  pull(PedigreeNew) %>%
  unique()
hyb17 <- yield %>%
  filter(Year == 2017) %>%
  pull(PedigreeNew) %>%
  unique() %>%
  setdiff(., hyb1416)

# Hybrids only grown in 2017
idx17 <- which(rownames(hyb) %in% hyb17)
geno17 <- hyb[idx17, ]

# Hybrids grown in 2014-2016
hyb <- hyb[-idx17, ]


# Identify hybrids for GxE analysis ---------------------------------------
ped_counts <- yield %>%
  filter(Year <= 2016, PedigreeNew %in% rownames(hyb)) %>%
  distinct(Site, PedigreeNew) %>%
  count(PedigreeNew)

count(ped_counts, n) %>%
  arrange(desc(n)) %>%
  mutate(nn = cumsum(nn)) %>%
  ggplot(., aes(x = n, y = nn)) + theme_bw() +
    geom_point(size = 3) + geom_vline(xintercept = 6, linetype = 2, colour = "red") +
    labs(x = "Number of Location-Years", y = "Number of Hybrids")
ggsave("figures/munge/hybrid_count.pdf", width = 6, height = 4, units = "in", dpi = 300)

yield <- yield %>%
  select(-Male, -Female) %>%
  mutate(Only17 = PedigreeNew %in% hyb17, 
         Geno = PedigreeNew %in% c(rownames(hyb), rownames(geno17)), 
         Obs = PedigreeNew %in% (filter(ped_counts, n >= 6) %>% pull(PedigreeNew)))

gbs_hyb <- filter(yield, Obs) %>% pull(PedigreeNew) %>% unique()
hyb <- hyb[rownames(hyb) %in% gbs_hyb, ]


# Apply QC filters --------------------------------------------------------
# Remove monomorphic SNPs
n_alleles <- apply(hyb, 2, function(x) length(unique(x)))
sum(n_alleles <= 1) # 6
hyb <- hyb[, n_alleles > 1]
geno17 <- geno17[, n_alleles > 1]
GM <- GM[n_alleles > 1, ]

# Remove SNPs with MAF < 0.025
maf <- apply(hyb, 2, function(x) sum(x)/(2*length(x)))
sum(maf < 0.025 | maf > 1 - 0.025) # 109,415
hyb <- hyb[, maf >= 0.025 & maf <= 1 - 0.025]
geno17 <- geno17[, maf >= 0.025 & maf <= 1 - 0.025]
GM <- GM[maf >= 0.025 & maf <= 1 - 0.025, ]

# Recode the minor allele based on synthetic hybrid genotypes
maf <- maf[names(maf) %in% GM$SNP]
for (i in seq_along(maf)) {
  if (maf[i] > 0.5) {
    idx0 <- which(hyb[, i] == 0)
    idx0.5 <- which(hyb[, i] == 0.5)
    idx1.5 <- which(hyb[, i] == 1.5)
    idx2 <- which(hyb[, i] == 2)
    
    hyb[idx0, i] <- 2
    hyb[idx0.5, i] <- 1.5
    hyb[idx1.5, i] <- 0.5
    hyb[idx2, i] <- 0
    
    maf[i] <- 1 - maf[i]
    GM$Alleles[i] <- str_split(GM$Alleles[i], ",") %>%
      unlist(use.names = FALSE) %>%
      rev() %>% paste(., collapse = ",")
  }
}

# Augment GM
GM$maf <- maf
ggplot(GM, aes(x = maf)) + theme_classic() +
  geom_histogram(binwidth = 0.01, fill = "orange", alpha = 0.8, colour = "black") +
  labs(x = "Minor Allele Frequency", y = "Count")
ggsave("figures/munge/maf.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Save the final tables ---------------------------------------------------
write_rds(list(GD = hyb, GD17 = geno17, GM = GM), "data/gbs/add_snps.rds")
write_rds(yield, "data/phenotype/yield_augment.rds")
