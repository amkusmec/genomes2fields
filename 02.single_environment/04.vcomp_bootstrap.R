### Bootstrap analysis of variability in expressed genetic variance for each 
### location-year.

library(tidyverse)
library(sommer)
library(parallel)


yield <- read_rds("data/phenotype/yield_agron0.rds") %>%
  filter(Year != 2017) %>%
  filter(Replicate != 0 | is.na(Replicate)) %>%
  mutate(Replicate = as.character(Replicate), 
         Block = paste(Replicate, Block, sep = "_"), 
         Rowf = as.character(Row), 
         Colf = as.character(Column)) %>%
  split(., .$Site)

vcomp <- read_rds("data/phenotype/variance_components.rds")
vcomp <- vcomp[!str_detect(names(vcomp), "2017$")]


cl <- makeCluster(32)
clusterEvalQ(cl, { library(tidyverse); library(sommer) })

boot <- clusterMap(cl, function(temp, v) {
  ped <- unique(temp$PedigreeNew)
  n_ped <- length(ped)
  ped_list <- map(ped, function(p) which(temp$PedigreeNew == p))
  names(ped_list) <- ped
  
  components <- matrix(0, nrow = 1000, ncol = length(v$sigma))
  colnames(components) <- names(v$sigma)
  
  pev <- numeric(1000)
  
  for (i in 1:1000) {
    idx <- unlist(ped_list[sample(ped, n_ped, replace = TRUE)], use.names = FALSE)
    m <- mmer(v$call$fixed, random = v$call$random, rcov = v$call$rcov, 
              data = temp[idx, ], date.warning = FALSE, verbose = FALSE)
    components[i, names(m$sigma)] <- m$sigma %>% unlist(use.names = FALSE)
    pev[i] <- mean(diag(m$PevU$`u:PedigreeNew`$Yield))
  }
  
  list(components = components, pev = pev)
}, temp = yield, v = vcomp, RECYCLE = FALSE, .scheduling = "static")
stopCluster(cl)
write_rds(boot, "data/phenotype/vcomp_bootstrap.rds")


# Analysis ----------------------------------------------------------------
# Heterogeneity in scaled genetic variance
gvar <- boot %>%
  map_df(function(l) {
    totals <- rowSums(l$components)
    tibble(VarG = median(l$components[, "u:PedigreeNew"]/totals, na.rm = TRUE),
           Lower = quantile(l$components[, "u:PedigreeNew"]/totals, probs = 0.025, na.rm = TRUE),
           Upper = quantile(l$components[, "u:PedigreeNew"]/totals, probs = 0.975, na.rm = TRUE))
  }) %>%
  mutate(Site = names(boot)) %>%
  dplyr::select(Site, everything())
ggplot(gvar, aes(x = Site)) + theme_classic() +
  geom_point(aes(y = VarG)) +
  geom_segment(aes(x = Site, xend = Site, y = Lower, yend = Upper)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.9)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = expression(paste("% ", sigma[g]^2)))
ggsave("figures/single/vg_het.pdf", width = 8, height = 5, units = "in", dpi = 300)

# # Heterogeneity in scaled residual variance
# rvar <- boot %>%
#   map_df(function(l) {
#     totals <- rowSums(l$components)
#     tibble(VarE = median(l$components[, "u:units"]/totals, na.rm = TRUE),
#            Lower = quantile(l$components[, "u:units"]/totals, probs = 0.025, na.rm = TRUE),
#            Upper = quantile(l$components[, "u:units"]/totals, probs = 0.975, na.rm = TRUE))
#   }) %>%
#   mutate(Site = names(boot)) %>%
#   dplyr::select(Site, everything())
# ggplot(rvar, aes(x = Site)) + theme_classic() +
#   geom_point(aes(y = VarE)) +
#   geom_segment(aes(x = Site, xend = Site, y = Lower, yend = Upper)) +
#   scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(x = "", y = expression(paste("% ", sigma[e]^2)))
# ggsave("figures/single/ve_het.pdf", width = 8, height = 5, units = "in", dpi = 300)

# # Heterogeneity in generalized heritability
# gh2 <- jack %>%
#   map_df(function(l) {
#     pev <- rowMeans(l$pev, na.rm = TRUE)
#     g <- 1 - pev/l$components[, "u:PedigreeNew"]
#     tibble(GenH2 = mean(g),
#            Lower = quantile(g, probs = 0.025, na.rm = TRUE),
#            Upper = quantile(g, probs = 0.975, na.rm = TRUE))
#   }) %>%
#   mutate(Site = names(jack)) %>%
#   dplyr::select(Site, everything()) %>%
#   filter(!is.nan(GenH2))
# ggplot(gh2, aes(x = Site)) + theme_classic() +
#   geom_point(aes(y = GenH2)) +
#   geom_segment(aes(x = Site, xend = Site, y = Lower, yend = Upper)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(x = "", y = expression(H[gen]^2))
# ggsave("figures/single/genH2.pdf", width = 8, height = 5, units = "in", dpi = 300)
