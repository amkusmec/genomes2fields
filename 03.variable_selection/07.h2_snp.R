library(tidyverse)
library(sommer)

rxn <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value)

snps <- read_rds("data/gbs/add_snps.rds")
K <- A.mat(snps$GD - 1)
D <- D.mat(snps$GD - 1)

idx <- match(rownames(K), rxn$PedigreeNew)
rxn <- rxn[idx, ]
pca <- read_rds("data/gbs/pca_covariates.rds")
rxn <- rxn %>%
  inner_join(., as_tibble(pca, rownames = "PedigreeNew"), by = "PedigreeNew") %>%
  mutate(PedigreeNewD = PedigreeNew)


p <- colnames(rxn)[2:8]
mA <- lapply(p, function(x) {
  form <- as.formula(paste0(x, " ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9"))
  mmer(form, random = ~ vs(PedigreeNew, Gu = K), rcov = ~ units, 
       data = rxn, date.warning = FALSE)
})
h2A <- lapply(mA, function(m2) pin(m2, h2 ~ V1/(V1 + V2)))
h2A <- do.call("rbind", h2A) %>%
  as_tibble() %>%
  mutate(Phenotype = p, 
         Model = "A")

mD <- lapply(p, function(x) {
  form <- as.formula(paste0(x, " ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9"))
  mmer(form, random = ~ vs(PedigreeNew, Gu = K) + 
         vs(PedigreeNewD, Gu = D), rcov = ~ units, 
       data = rxn, date.warning = FALSE)
})
h2D <- lapply(mD, function(m2) pin(m2, h2 ~ (V1 + V2)/(V1 + V2 + V3)))
h2D <- do.call("rbind", h2D) %>%
  as_tibble() %>%
  mutate(Phenotype = p, 
         Model = "A+D")

h2 <- bind_rows(h2A, h2D)

ggplot(h2, aes(x = Phenotype, y = Estimate, group = Model, colour = Model)) +
  geom_pointrange(aes(ymin = Estimate - SE, ymax = Estimate + SE), 
                  position = position_dodge(width = 0.4)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), 
        panel.grid.major.y = element_line(colour = "grey80", size = 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_colour_manual(values = c("A" = "blue", "A+D" = "orange")) +
  labs(x = "", y = expression(h[SNP]^2))
ggsave("figures/select/rxn_h2.pdf", width = 6, height = 4, units = "in", dpi = 300)
