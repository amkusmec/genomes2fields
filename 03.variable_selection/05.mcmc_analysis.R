### Process and summarize the results of the random regression model.

library(tidyverse)
library(MCMCglmm)
library(abind)


# Collect and process the output of each chain ----------------------------
# g2 <- read_rds("data/weather/ga_final_model.rds")
sel <- read_lines("data/phenotype/selected_variables.txt")
output <- list.files("data/phenotype", "mcmc_chain[1-8]\\.rds", full.names = TRUE) %>%
  map(function(f) {
    temp <- read_rds(f)
    idx1 <- which(!str_detect(colnames(temp$Sol), "meta"))
    idx2 <- which(!str_detect(colnames(temp$VCV), "meta"))
    mcmc(cbind(temp$Sol[, idx1], temp$VCV[, idx2]), 
         start = 500001, end = 2500000, thin = 2000)
  }) %>%
  as.mcmc.list()


# Convergence analysis ----------------------------------------------------
# 3,754 parameters total
gdiag <- geweke.diag(output)
idx <- lapply(gdiag, function(g) which(abs(g$z) < 1.96))
length(reduce(idx, intersect)) # 2,277 (60.7%) converged in all chains
length(reduce(idx, union)) # 3,754 (100%) converged in at least one chain
unlist(idx) %>% enframe() %>% count(value) %>% count(n)
# All parameters converged in at least 4 chains
# < 100 parameters converged in < 6 chains


# Format posterior means for GWAS -----------------------------------------
eff_size <- effectiveSize(output)
# Minimum ESS is 7,117

gibbs <- lapply(output, as.matrix)
gibbs <- do.call("rbind", gibbs)
params <- enframe(colMeans(gibbs))

# Extract the fixed effects
fixed <- filter(params, name %in% sel)

# Extract the random effects
rand <- filter(params, str_detect(name, "PedigreeNew"), 
               !str_detect(name, "units"), !str_detect(name, ":")) %>%
  mutate(name = str_replace(name, "\\.PedigreeNew\\.", ":")) %>%
  separate(name, c("Parameter", "PedigreeNew"), sep = ":", remove = TRUE) %>%
  full_join(., fixed, by = c("Parameter" = "name")) %>%
  mutate(Value = if_else(is.na(value.y), value.x, value.x + value.y)) %>%
  select(-value.x, -value.y)

write_rds(rand, "data/phenotype/glmm_rxn_norm_parameters.rds")

pA <- ggplot(rand, aes(x = Value)) + theme_bw() + 
  geom_density() + facet_wrap(~ Parameter, scales = "free") + 
  labs(x = "", y = "Density")
ggsave("figures/select/rxn_norm_density.pdf", plot = pA, width = 8, 
       height = 5, units = "in", dpi = 300)


# Genetic covariance matrix -----------------------------------------------
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

vcov <- filter(params, str_detect(name, ":")) %>%
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
ggsave("figures/select/gen_cor.pdf", plot = pB, width = 6, 
       height = 4, units = "in", dpi = 300)

# library(grid)
# library(gridExtra)
# 
# lay <- matrix(c(1, 1, 1, 2, 2), nrow = 1)
# grob1 <- grobTree(ggplotGrob(pA), 
#                   textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
#                            hjust = "left", vjust = "top", 
#                            gp = gpar(fontface = "bold", fontsize = 14)))
# grob2 <- grobTree(ggplotGrob(pB), 
#                   textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
#                            hjust = "left", vjust = "top", 
#                            gp = gpar(fontface = "bold", fontsize = 14)))
# gp <- arrangeGrob(grob1, grob2, layout_matrix = lay)
# ggsave("figures/select/rxn_norm_plots.pdf", plot = gp, width = 15, height = 5, 
#        units = "in", dpi = 300)
