### Various exploratory analyses of the PCA results to examine the population 
### structure. It is easy to see structure induced by the inclusion of the 
### NAM RILs and by the crossing structure of the experiments.

library(tidyverse)


# Separate into groups based on PC1 < 0 -----------------------------------
pca <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(rownames = "Pedigree") %>%
  separate(Pedigree, c("Parent1", "Parent2"), sep = "/", remove = FALSE) %>%
  mutate(Group = if_else(PC1 < 125, "1", "2"))

table(pca$Group)

v_exp <- summary(read_rds("data/gbs/pca.rds"))$importance


# Scatter plots -----------------------------------------------------------
ggplot(pca, aes(x = PC1, y = PC2, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC1 (17.3%)", y = "PC2 (5.6%)")
ggsave("figures/munge/pc1_2.pdf", width = 6, height = 4, units = "in", dpi = 300)

ggplot(pca, aes(x = PC1, y = PC3, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC1 (17.3%)", y = "PC3 (5.0%)")

ggplot(pca, aes(x = PC1, y = PC4, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC1 (17.3%)", y = "PC4 (4.1%)")

ggplot(pca, aes(x = PC2, y = PC3, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC2 (5.6%)", y = "PC3 (5.0%)")

ggplot(pca, aes(x = PC2, y = PC4, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC2 (5.6%)", y = "PC4 (4.1%)")
ggsave("figures/munge/pc2_4.pdf", width = 6, height = 4, units = "in", dpi = 300)

ggplot(pca, aes(x = PC3, y = PC4, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC3 (5.0%)", y = "PC4 (4.1%)")
ggsave("figures/munge/pc3_4.pdf", width = 6, height = 4, units = "in", dpi = 300)


# Parents in the two groups -----------------------------------------------
parents1 <- c(filter(pca, Group == "1") %>% pull(Parent1), 
              filter(pca, Group == "1") %>% pull(Parent2)) %>%
  unique()
parents2 <- c(filter(pca, Group == "2") %>% pull(Parent1), 
              filter(pca, Group == "2") %>% pull(Parent2)) %>%
  unique()
all_parents <- c(pull(pca, Parent1), pull(pca, Parent2)) %>% unique()

png("figures/munge/pca_parents.png", width = 8, height = 6, units = "in", res = 300)
VennDiagram::draw.pairwise.venn(area1 = length(parents1), area2 = length(parents2), 
                                cross.area = length(intersect(parents1, parents2)), 
                                lwd = rep(3, 2), fill = c("red", "blue"), 
                                fontfamily = rep("sanserif", 3), cex = rep(1.5, 3), 
                                cat.fontfamily = rep("sanserif", 2), cat.pos = rep(0, 2), 
                                category = c("Group 1", "Group 2"), 
                                cat.cex = rep(2, 2), cat.dist = rep(0.05, 2))
dev.off()

intersect(parents1, parents2)
setdiff(parents2, parents1)

all_parents <- c(pull(pca, Parent1), pull(pca, Parent2))
all_parents <- all_parents[all_parents %in% intersect(parents2, parents1)]
table(all_parents)

# These are the RILs from Z013 (KI3) and Z022 (OH43)
# KI3 is not included as a parent for other hybrids; OH43 is.


# Check minor allele frequencies without NAM RILs -------------------------
gbs <- read_rds("data/gbs/add_snps.rds")$GD
nam <- which(str_detect(rownames(gbs), "Z[0-9]{3}E[0-9]{4}"))

ngbs <- gbs[nam, ]
ogbs <- gbs[-nam, ]
rm(gbs); gc()

nmaf <- apply(ngbs, 2, function(x) sum(x)/(2*length(x)))
omaf <- apply(ogbs, 2, function(x) sum(x)/(2*length(x)))
cor(nmaf, omaf)

tibble(Type = rep(c("NAM", "Other"), each = length(nmaf)), 
       MAF = c(nmaf, omaf)) %>%
  ggplot(., aes(x = MAF, group = Type, fill = Type)) + theme_classic() +
    geom_histogram(alpha = 0.6, binwidth = 0.01, colour = "black") + 
    labs(x = "Minor Allele Frequency", y = "Count", fill = "") +
    scale_fill_manual(values = c("NAM" = "blue", "Other" = "red")) +
    facet_wrap(~ Type, ncol = 1, strip.position = "right")
ggsave("figures/munge/structure_maf.pdf", width = 6, height = 4, units = "in", dpi = 300)

summary(nmaf)
summary(omaf)

sum(omaf == 0)
sum(omaf == 0)/length(omaf)

sum(nmaf == 0)
sum(nmaf == 0)/length(nmaf)

hmp2_snps <- read_tsv("data/gbs/HapMap2/maizeHapMapV2_recode.vcf", comment = "##") %>%
  mutate(ID = paste0("X", ID)) %>%
  pull(ID)
hmp2_snps <- hmp2_snps[!duplicated(hmp2_snps)]

sum(names(omaf)[omaf == 0] %in% hmp2_snps)


# Inter-SNP distance ------------------------------------------------------
GM <- read_rds("data/gbs/add_snps.rds")$GM

GM <- GM %>%
  split(., .$Chromosome) %>%
  map_df(function(df) {
    mutate(df, Distance = c(0, Position[-1] - Position[-length(Position)]))
  })

GM %>%
  mutate(Chromosome = factor(Chromosome), 
         Distance = log10(Distance)) %>%
  ggplot(., aes(x = Distance)) + theme_classic() +
    geom_density(aes(colour = Chromosome, group = Chromosome)) +
    labs(y = "Density", x = expression(paste(log[10], "(Distance)")))
ggsave("figures/munge/snp_distance.pdf", width = 6, height = 4, units = "in", dpi = 300)

GM %>% split(., .$Chromosome) %>%
  map_dbl(function(df) quantile(df$Distance, probs = 0.95))


# Non-segregating SNPs and PC loadings ------------------------------------
# Identify non-segregating sites
nseg_nam <- which(nmaf == 0)
nseg_oth <- which(omaf == 0)
nseg <- union(nseg_nam, nseg_oth) %>% sort()

# Genomic distribution of non-segregating SNPs
chrom <- GM %>%
  group_by(Chromosome) %>%
  summarise(Length = max(Position)/1e6) %>%
  ungroup()

GM <- GM %>%
  mutate(NSeg = if_else(row_number() %in% nseg_nam, "NAM", 
                        if_else(row_number() %in% nseg_oth, "G2F", "Segregating")))

GM %>% mutate(Position = Position/1e6) %>%
  filter(NSeg != "Segregating") %>%
  ggplot() + theme_classic() + labs(x = "Chromosome", y = "Position (Mbp)", colour = "") +
    geom_rect(aes(xmin = Chromosome - 0.25, xmax = Chromosome + 0.25, 
                  ymax = Length + 100/1e6), data = chrom, ymin = 0, colour = "black", 
              fill = "white") +
    geom_segment(aes(x = Chromosome - 0.2, xend = Chromosome + 0.2, 
                     y = Position, yend = Position, colour = NSeg), size = 0.1) +
    scale_colour_manual(values = c("NAM" = "#1F78B4", "G2F" = "#B2DF8A")) +
    scale_x_continuous(breaks = 1:10)
ggsave("figures/munge/nonsegregating_distribution.pdf", width = 8, height = 6, 
       units = "in", dpi = 300)

# Density quantiles for sub-populations without non-segregating SNPs
# NAM
GM %>% 
  filter(NSeg != "NAM") %>%
  split(., .$Chromosome) %>%
  map_dbl(function(df) quantile(df$Distance, probs = 0.95))

# G2F
GM %>%
  filter(NSeg != "G2F") %>%
  split(., .$Chromosome) %>%
  map_dbl(function(df) quantile(df$Distance, probs = 0.95))


pc <- read_rds("data/gbs/pca.rds")
pc1_loadings <- pc$rotation[, 1]

# Fun but relatively uninformative plot
GM %>%
  mutate(PC1 = pc1_loadings, 
         Position = Position/1e6) %>%
  ggplot(.) + theme_classic() + 
    geom_segment(aes(x = Position, y = 0, xend = Position, yend = PC1, 
                     colour = NSeg, size = NSeg)) +
    scale_colour_manual(values = c("NAM" = "#1F78B4", "G2F" = "#B2DF8A", 
                                   "Segregating" = "grey80")) +
    scale_size_manual(values = c("NAM" = 0.3, "G2F" = 0.3, "Segregating" = 0.1)) +
    labs(x = "Position (Mbp)", y = "PC1 Loading") + facet_wrap(~ Chromosome)

GM %>% mutate(PC1 = pc1_loadings) %>%
  ggplot(., aes(x = PC1, group = NSeg, colour = NSeg)) + theme_classic() +
    geom_density() + facet_wrap(~ Chromosome, scales = "free") +
    scale_colour_manual(values = c("NAM" = "#1F78B4", "G2F" = "#B2DF8A", 
                                   "Segregating" = "grey80")) +
    labs(x = "PC1 Loading", y = "Density", colour = "")
ggsave("figures/munge/pc1_loadings_density.pdf", width = 10, height = 6, 
       units = "in", dpi = 300)


# k-means clustering ------------------------------------------------------
group1 <- pca %>%
  filter(Group == "1") %>%
  select(Pedigree:Parent2, PC2, PC4)

km <- sapply(2:10, function(x) kmeans(group1[, 4:5], centers = x, 
                                      iter.max = 50, nstart = 10), 
             simplify = FALSE)
opt_clust <- which.min(map_dbl(km, function(x) x$tot.withinss))

# Optimal number of clusters is 10 (at least), but 6 is easier to illustrate 
# the relevant structure in the population.
group1 <- mutate(group1, Cluster = km[[5]]$cluster)

group1 %>%
  mutate(Cluster = as.character(Cluster)) %>%
  ggplot(., aes(x = PC2, y = PC4, colour = Cluster)) + theme_classic() +
    geom_point() + scale_color_brewer(type = "qual", palette = "Dark2") +
    labs(x = "PC2 (5.6%)", y = "PC4 (4.1%)")
ggsave("figures/munge/clustered_pc2_4.pdf", width = 6, height = 4, units = "in", 
       dpi = 300)

p1 <- with(group1, table(Parent1, Cluster))
p2 <- with(group1, table(Parent2, Cluster))
