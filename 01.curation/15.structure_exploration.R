library(tidyverse)


# Separate into groups based on PC1 < 0 -----------------------------------
pca <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(rownames = "Pedigree") %>%
  separate(Pedigree, c("Parent1", "Parent2"), sep = "/", remove = FALSE) %>%
  mutate(Group = if_else(PC1 < 0, "1", "2"))

table(pca$Group)

v_exp <- summary(read_rds("data/gbs/pca.rds"))$importance


# Scatter plots -----------------------------------------------------------
ggplot(pca, aes(x = PC1, y = PC2, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC1 (30.4%)", y = "PC2 (4.9%)")
ggsave("figures/munge/pc1_2.pdf", width = 6, height = 4, units = "in", dpi = 300)

ggplot(pca, aes(x = PC1, y = PC3, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC1 (30.4%)", y = "PC3 (3.8%)")

ggplot(pca, aes(x = PC1, y = PC4, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC1 (30.4%)", y = "PC4 (3.6%)")

ggplot(pca, aes(x = PC2, y = PC3, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC2 (4.9%)", y = "PC3 (3.8%)")

ggplot(pca, aes(x = PC2, y = PC4, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC2 (4.9%)", y = "PC4 (3.6%)")

ggplot(pca, aes(x = PC3, y = PC4, colour = Group)) + theme_classic() +
  geom_point(alpha = 0.8) + geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + guides(colour = FALSE) +
  scale_colour_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "PC3 (3.8%)", y = "PC4 (3.6%)")
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
    facet_wrap(~ Type, ncol = 1)
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
