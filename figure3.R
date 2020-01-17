library(tidyverse)
library(coda)
library(BGLR)
library(abind)
library(grid)
library(gridExtra)


# A: phenotypic distributions ---------------------------------------------
rxn <- read_rds("data/phenotype/glmm_rxn_norm_parameters.rds") %>%
  mutate(Parameter = str_replace(Parameter, "\\(Intercept\\)", "Hybrid"))
pA <- ggplot(rxn, aes(x = Value)) + theme_bw() + 
  geom_density() + facet_wrap(~ Parameter, scales = "free") + 
  labs(x = "", y = "Density")


# B: genetic correlations -------------------------------------------------
output <- list.files("data/phenotype", "mcmc_chain[1-8]\\.rds", full.names = TRUE) %>%
  map(function(f) {
    temp <- read_rds(f)
    idx1 <- which(!str_detect(colnames(temp$Sol), "meta"))
    idx2 <- which(!str_detect(colnames(temp$VCV), "meta"))
    mcmc(cbind(temp$Sol[, idx1], temp$VCV[, idx2]), 
         start = 500001, end = 2500000, thin = 2000)
  }) %>%
  as.mcmc.list()

gibbs <- lapply(output, as.matrix)
gibbs <- do.call("rbind", gibbs)

vcov_gibbs <- gibbs[, str_detect(colnames(gibbs), ":")]
cor_list <- lapply(1:nrow(vcov_gibbs), function(j) {
  x <- vcov_gibbs[j, ]
  m <- matrix(x, nrow = sqrt(length(x)), ncol = sqrt(length(x)), byrow = TRUE)
  temp <- str_remove(names(x)[1:5], "\\.PedigreeNew") %>%
    str_split(., ":") %>% sapply(., function(i) i[1])
  dimnames(m) <- list(temp, temp)
  m <- m[order(rownames(m)), order(colnames(m))]
  cov2cor(m)
})
cor_list <- abind(cor_list, along = 3)

ql <- apply(cor_list, c(1, 2), quantile, probs = 0.025)
ql[upper.tri(ql, diag = TRUE)] <- NA
ql <- as_tibble(ql, rownames = "Phenotype1") %>%
  gather(Phenotype2, Lower, -Phenotype1) %>%
  filter(!is.na(Lower)) %>%
  mutate(Phenotype2 = str_replace(Phenotype2, "\\(Intercept\\)", "Hybrid"))

qu <- apply(cor_list, c(1, 2), quantile, probs = 0.975)
qu[upper.tri(qu, diag = TRUE)] <- NA
qu <- as_tibble(qu, rownames = "Phenotype1") %>%
  gather(Phenotype2, Upper, -Phenotype1) %>%
  filter(!is.na(Upper)) %>%
  mutate(Phenotype2 = str_replace(Phenotype2, "\\(Intercept\\)", "Hybrid"))

limits <- inner_join(ql, qu, by = c("Phenotype1", "Phenotype2")) %>%
  mutate(Lower = round(Lower, 2), 
         Upper = round(Upper, 2), 
         Label = paste0("(", Lower, ",", Upper, ")")) %>%
  arrange(Phenotype1, Phenotype2)

vcov <- enframe(colMeans(vcov_gibbs)) %>%
  mutate(name = str_remove(name, "\\.PedigreeNew")) %>%
  separate(name, c("Phenotype1", "Phenotype2"), sep = ":", remove = TRUE) %>%
  spread(Phenotype2, value)
vcov <- as.matrix(vcov[, -1])
rownames(vcov) <- colnames(vcov)

corg <- cov2cor(vcov)
corg[upper.tri(corg, diag = TRUE)] <- NA
corg <- as_tibble(corg, rownames = "Phenotype1") %>%
  gather(Phenotype2, Mean, -Phenotype1) %>%
  filter(!is.na(Mean)) %>%
  mutate(Phenotype2 = str_replace(Phenotype2, "\\(Intercept\\)", "Hybrid")) %>%
  arrange(Phenotype1, Phenotype2)

corg <- inner_join(corg, limits, by = c("Phenotype1", "Phenotype2")) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE))

pB <- ggplot(corg, aes(x = Phenotype1, y = Phenotype2)) + theme_classic() + 
  geom_tile(aes(fill = Mean)) + labs(x = "", y = "", fill = expression(r[G]^2)) +
  geom_text(aes(label = Label), size = 2.5) + 
  scale_fill_gradient2(low = "blue", high = "red") + 
  scale_x_discrete(position = "top") + 
  theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))


# C: heritabilities -------------------------------------------------------
rxn <- spread(rxn, Parameter, Value) %>%
  rename(`(Intercept)` = Hybrid)
est <- map_df(2:6, function(p) {
  varE_A <- scan(paste0("data/bglr/h2/eigA_", names(rxn)[p], "_varE.dat"))[-c(1:400)]
  varU_A <- scan(paste0("data/bglr/h2/eigA_", names(rxn)[p], "_ETA_A_varU.dat"))[-c(1:400)]
  
  varE_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_varE.dat"))[-c(1:400)]
  varUA_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_ETA_A_varU.dat"))[-c(1:400)]
  varUD_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn)[p], "_ETA_D_varU.dat"))[-c(1:400)]
  
  tibble(Phenotype = names(rxn)[p], 
         Model = rep(c("A", "A+D"), each = length(varE_A)), 
         h2 = c(varU_A/(varU_A + varE_A), 
                (varUA_D + varUD_D)/(varUA_D + varUD_D + varE_D)))
}) %>% 
  mutate(Phenotype = str_replace(Phenotype, "\\(Intercept\\)", "Hybrid")) %>%
  group_by(Phenotype, Model) %>%
  summarise(Mean = mean(h2), 
            Lower = quantile(h2, probs = 0.025), 
            Upper = quantile(h2, probs = 0.975)) %>%
  ungroup() %>%
  mutate(Model = factor(Model), 
         Phenotype = factor(Phenotype, levels = unique(Phenotype), ordered = TRUE))

pC <- ggplot(est, aes(x = Phenotype, y = Mean, group = Model)) + theme_bw() +
  geom_pointrange(aes(ymin = Lower, ymax = Upper, colour = Model), 
                  position = position_dodge(0.5)) +
  scale_colour_manual(values = c("A" = "skyblue", "A+D" = "palegreen")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) + 
  labs(x = "", y = expression(h[SNP]^2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank())


# Final figure ------------------------------------------------------------
lay <- matrix(c(1, 1, 1, 1, 2, 2, 3, 3), nrow = 2, byrow = TRUE)
gA <- grobTree(ggplotGrob(pA), 
               textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gB <- grobTree(ggplotGrob(pB), 
               textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gC <- grobTree(ggplotGrob(pC), 
               textGrob("C", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                        hjust = "left", vjust = "top", 
                        gp = gpar(fontface = "bold", fontsize = 14)))
gp <- arrangeGrob(gA, gB, gC, layout_matrix = lay)
ggsave("figures/select/rxn_norm_plots.pdf", plot = gp, width = 14, 
       height = 10, units = "in", dpi = 300)
