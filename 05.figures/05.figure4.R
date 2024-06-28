library(tidyverse)
library(GGally)
library(coda)
library(BGLR)
library(abind)
library(grid)
library(gridExtra)


# A: phenotypic distributions ---------------------------------------------
pop_params <- list.files("data/phenotype", "mcmc_chain[1-8]\\.rds", full.names = TRUE) %>%
  map(function(f) read_rds(f)$Sol[, c(1, 46:49)])
pop_params <- do.call("rbind", pop_params)

sel <- read_lines("data/phenotype/selected_variables.txt")
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  filter(Site != "NEH3_2015")
data <- data[, sel]
ec_means <- colMeans(data)
ec_means <- cbind(matrix(rep(1, nrow(pop_params)), ncol = 1), 
                  matrix(ec_means, ncol = 4, nrow = nrow(pop_params), byrow = TRUE))
# int <- mean(rowSums(pop_params*ec_means))
int <- mean(pop_params[, 1])

ec_ranges <- apply(data, 2, function(x) range(x) - mean(x)) %>% 
  t() %>% 
  as_tibble(rownames = "Parameter") %>% 
  rename(Min = V1, Max = V2) %>% 
  mutate(Parameter2 = str_replace(Parameter, "\\(Intercept\\)", "Intercept") %>%
           str_replace("TMAX_X0.825_X0.95", "Pre-anthesis max. temp.") %>%
           str_replace("TMAX_X0.875_X1.425", "Post-anthesis max. temp.") %>%
           str_replace("SR_X0.65_X1.425", "Mid-season solar radiation") %>%
           str_replace("NET_X0.175_X1.4", "Whole season net evapotranspiration"))

# EC responses have already had the fixed effects added to them
# Intercepts require adjustment
rxn <- read_rds("data/phenotype/glmm_rxn_norm_parameters.rds") %>%
  mutate(Parameter2 = str_replace(Parameter, "\\(Intercept\\)", "Intercept") %>%
           str_replace("TMAX_X0.825_X0.95", "Pre-anthesis max. temp.") %>%
           str_replace("TMAX_X0.875_X1.425", "Post-anthesis max. temp.") %>%
           str_replace("SR_X0.65_X1.425", "Mid-season solar radiation") %>%
           str_replace("NET_X0.175_X1.4", "Whole season net evapotranspiration"), 
         Value = if_else(Parameter2 == "Intercept", Value + int, Value), 
         Value = Value*62.77) # Convert bu/ac => kg/ha

rxn3 <- read_rds("data/phenotype/glmm_rxn_norm_parameters.rds") %>%
  mutate(Parameter2 = str_replace(Parameter, "\\(Intercept\\)", "Average yield") %>%
           str_replace("TMAX_X0.825_X0.95", "Pre-anthesis\nmax. temp.") %>%
           str_replace("TMAX_X0.875_X1.425", "Post-anthesis\nmax. temp.") %>%
           str_replace("SR_X0.65_X1.425", "Mid-season\nsolar radiation") %>%
           str_replace("NET_X0.175_X1.4", "Whole season net\nevapotranspiration"), 
         Value = if_else(Parameter2 == "Intercept", Value + int, Value), 
         Value = Value*62.77) %>% # Convert bu/ac => kg/ha 
  select(-Parameter) %>% 
  pivot_wider(names_from = "Parameter2", values_from = "Value")

pp <- ggpairs(rxn3, columns = 2:6, progress = FALSE, 
        upper = list(continuous = "blank"), 
        lower = list(continuous = wrap(ggally_points, alpha = 0.5)), 
        diag = list(continuous = wrap(ggally_barDiag, bins = 30, colour = "black", 
                                      fill = "grey60"))) + 
  theme_classic()
for (i in 2:5) {
  for (j in 1:(i - 1)) {
    pp[i, j] <- pp[i, j] + geom_hline(yintercept = 0, linetype = 2, colour = "red") + 
      geom_vline(xintercept = 0, linetype = 2, colour = "red")
    
    if (i < 5) {
      pp[i, j] <- pp[i, j] + annotate("rect", xmin = c(0, -Inf, -Inf, 0), 
                                      ymin = c(0, 0, -Inf, -Inf), 
                                      xmax = c(Inf, 0, 0, Inf), 
                                      ymax = c(Inf, Inf, 0, 0), 
                                      fill = c("green", "orange", "red", "orange"), 
                                      alpha = 0.2)
    } else {
      pp[i, j] <- pp[i, j] + annotate("rect", xmin = c(0, -Inf, -Inf, 0), 
                                      ymin = c(0, 0, -Inf, -Inf), 
                                      xmax = c(Inf, 0, 0, Inf), 
                                      ymax = c(Inf, Inf, 0, 0), 
                                      fill = c("orange", "red", "orange", "green"), 
                                      alpha = 0.2)
    }
  }
}
png("figures/select/rxn_pairs.png", width = 10.5, height = 7, units = "in", res = 300)
pp
dev.off()

rxn2 <- rxn %>% 
  filter(Parameter2 != "Intercept") %>% 
  select(-Parameter) %>% 
  group_by(Parameter2) %>% 
  arrange(desc(Value)) %>% 
  slice(c(1:5, (n() - 4):n())) %>% 
  ungroup() %>% 
  inner_join(rxn %>% filter(Parameter2 == "Intercept") %>% 
               select(PedigreeNew, Value) %>% rename(Intercept = Value), 
             by = "PedigreeNew") %>% 
  inner_join(ec_ranges, by = "Parameter2") %>% 
  mutate(YMin = Intercept + Value*Min, 
         YMax = Intercept + Value*Max)

ylim <- range(c(rxn2$YMin, rxn2$YMax))
ylim <- c(ylim[1] - abs(0.04*ylim[1]), ylim[2] + 0.04*ylim[2])
clim <- c(0, 120)

rxn_plots <- list(
  ggplot(rxn2[1:10, ]) + theme_classic() + 
    geom_segment(aes(x = Min, y = YMin, xend = Max, yend = YMax, 
                     group = PedigreeNew), linewidth = 0.25) + 
    scale_y_continuous(limits = ylim) + 
    labs(y = expression("Grain yield (kg "*ha^-1*")"), tag = "(a)", 
         x = expression("Mid-season solar rad. (MJ "*m^-2*")")), 
  ggplot(rxn2[11:20, ]) + theme_classic() + 
    geom_segment(aes(x = Min, y = YMin, xend = Max, yend = YMax, 
                     group = PedigreeNew), linewidth = 0.25) + 
    scale_y_continuous(limits = ylim) + 
    labs(y = expression("Grain yield (kg "*ha^-1*")"), tag = " ", 
         x = expression(paste("Post-anthesis max. temp. (", degree, "C)"))), 
  ggplot(rxn2[21:30, ]) + theme_classic() + 
    geom_segment(aes(x = Min, y = YMin, xend = Max, yend = YMax, 
                     group = PedigreeNew), linewidth = 0.25) + 
    scale_y_continuous(limits = ylim) + 
    labs(y = expression("Grain yield (kg "*ha^-1*")"), tag = " ",
         x = expression(paste("Pre-anthesis max. temp. (", degree, "C)"))), 
  ggplot(rxn2[31:40, ]) + theme_classic() + 
    geom_segment(aes(x = Min, y = YMin, xend = Max, yend = YMax, 
                     group = PedigreeNew), linewidth = 0.25) + 
    scale_y_continuous(limits = ylim) + 
    labs(y = expression("Grain yield (kg "*ha^-1*")"), tag = " ",
         x = "Seasonal net evapotranspiration (mm)"), 
  ggplot(filter(rxn, Parameter2 == "Mid-season solar radiation")) + theme_classic() + 
    geom_histogram(aes(x = Value), colour = "black", fill = "brown", alpha = 0.75) + 
    scale_y_continuous(limits = clim) + 
    labs(x = expression("kg"~ha^-1~MJ^-1~m^-2), y = "# hybrids", tag = " "), 
  ggplot(filter(rxn, Parameter2 == "Post-anthesis max. temp.")) + theme_classic() + 
    geom_histogram(aes(x = Value), colour = "black", fill = "red", alpha = 0.75) + 
    scale_y_continuous(limits = clim) + 
    labs(x = expression("kg"~ha^-1~degree*C^-1), y = "# hybrids", tag = " "), 
  ggplot(filter(rxn, Parameter2 == "Pre-anthesis max. temp.")) + theme_classic() + 
    geom_histogram(aes(x = Value), colour = "black", fill = "red", alpha = 0.75) + 
    scale_y_continuous(limits = clim) + 
    labs(x = expression("kg"~ha^-1~degree*C^-1), y = "# hybrids", tag = " "), 
  ggplot(filter(rxn, Parameter2 == "Whole season net evapotranspiration")) + 
    theme_classic() + 
    geom_histogram(aes(x = Value), colour = "black", fill = "orange", alpha = 0.75) + 
    scale_y_continuous(limits = clim) + 
    labs(x = expression("kg"~ha^-1~mm^-1), y = "# hybrids", tag = " ")
)

rxn_plots <- lapply(rxn_plots, function(p) ggplotGrob(p + theme(axis.text = element_text(size = 4), axis.title = element_text(size = 5))))
# plot(arrangeGrob(grobs = rxn_plots, layout_matrix = matrix(1:8, ncol = 2)))


# B: genetic correlations -------------------------------------------------
phenos <- c("Intercept", "Whole season net evapotranspiration", "Mid-season solar radiation", 
            "Pre-anthesis max. temp.", "Post-anthesis max. temp.")
phenos2 <- c("Intercept", "Whole season net\nevapotranspiration", "Mid-season\nsolar radiation", 
            "Pre-anthesis\nmax. temp.", "Post-anthesis\nmax. temp.")
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
dimnames(cor_list)[[1]] <- phenos2
dimnames(cor_list)[[2]] <- phenos2

ql <- apply(cor_list, c(1, 2), quantile, probs = 0.025)
ql[upper.tri(ql, diag = TRUE)] <- NA
ql <- as_tibble(ql, rownames = "Phenotype1") %>%
  gather(Phenotype2, Lower, -Phenotype1) %>%
  filter(!is.na(Lower))

qu <- apply(cor_list, c(1, 2), quantile, probs = 0.975)
qu[upper.tri(qu, diag = TRUE)] <- NA
qu <- as_tibble(qu, rownames = "Phenotype1") %>%
  gather(Phenotype2, Upper, -Phenotype1) %>%
  filter(!is.na(Upper))

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
rownames(vcov) <- colnames(vcov) <- phenos2

corg <- cov2cor(vcov)
corg[upper.tri(corg, diag = TRUE)] <- NA
corg <- as_tibble(corg, rownames = "Phenotype1") %>%
  gather(Phenotype2, Mean, -Phenotype1) %>%
  filter(!is.na(Mean)) %>%
  arrange(Phenotype1, Phenotype2)

corg <- inner_join(corg, limits, by = c("Phenotype1", "Phenotype2")) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = phenos2, ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(phenos2), ordered = TRUE), 
         # Label = paste0(round(Mean, digits = 2), 
         #                if_else(sign(Lower) == sign(Upper), "*", ""), 
         #                "\n", Label))
         Label = paste0(round(Mean, digits = 2), 
                        if_else(sign(Lower) == sign(Upper), "*", "")))

pB <- ggplot(corg, aes(x = Phenotype1, y = Phenotype2)) + theme_classic() + 
  geom_tile(aes(fill = Mean)) + 
  geom_label(aes(label = Label), size = 2) + 
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) + 
  scale_x_discrete(position = "top") + 
  labs(x = "", y = "", fill = expression(r[G]^2), tag = "(b)") +
  theme(axis.text.x = element_text(angle = 90, size = 7), 
        axis.text.y = element_text(size = 7), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6), 
        legend.direction = "horizontal", 
        legend.position = "inside", 
        legend.position.inside = c(0.3, 0.075), 
        legend.key.width = unit(0.05, "npc"), 
        legend.key.height = unit(0.025, "npc"))


# C: heritabilities -------------------------------------------------------
rxn4 <- spread(rxn, Parameter, Value)
est <- map_df(3:7, function(p) {
  varE_A <- scan(paste0("data/bglr/h2/eigA_", names(rxn4)[p], "_varE.dat"))[-c(1:400)]
  varU_A <- scan(paste0("data/bglr/h2/eigA_", names(rxn4)[p], "_ETA_A_varU.dat"))[-c(1:400)]
  
  varE_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn4)[p], "_varE.dat"))[-c(1:400)]
  varUA_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn4)[p], "_ETA_A_varU.dat"))[-c(1:400)]
  varUD_D <- scan(paste0("data/bglr/h2/eigD_", names(rxn4)[p], "_ETA_D_varU.dat"))[-c(1:400)]
  
  tibble(Phenotype = names(rxn4)[p], 
         Model = rep(c("A", "A+D"), each = length(varE_A)), 
         h2 = c(varU_A/(varU_A + varE_A), 
                (varUA_D + varUD_D)/(varUA_D + varUD_D + varE_D)))
}) %>% 
  mutate(Phenotype = str_replace(Phenotype, "\\(Intercept\\)", "Intercept") %>%
           str_replace("TMAX_X0.825_X0.95", "Pre-anthesis\nmax. temp.") %>%
           str_replace("TMAX_X0.875_X1.425", "Post-anthesis\nmax. temp.") %>%
           str_replace("SR_X0.65_X1.425", "Mid-season\nsolar radiation") %>%
           str_replace("NET_X0.175_X1.4", "Whole season net\nevapotranspiration")) %>%
  group_by(Phenotype, Model) %>%
  summarise(Mean = mean(h2), 
            Lower = quantile(h2, probs = 0.025), 
            Upper = quantile(h2, probs = 0.975)) %>%
  ungroup() %>%
  mutate(Model = factor(Model), 
         Phenotype = factor(Phenotype, levels = phenos2, ordered = TRUE))

pC <- ggplot(est, aes(x = Phenotype, y = Mean, group = Model)) + theme_classic() +
  geom_pointrange(aes(ymin = Lower, ymax = Upper, colour = Model), 
                  position = position_dodge(0.5)) +
  scale_colour_brewer(type = "qual", palette = "Dark2") + 
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) + 
  labs(x = "", y = expression(h[SNP]^2), tag = "(c)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(size = 7), 
        axis.title = element_text(size = 10), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        legend.title = element_text(size = 9), 
        legend.text = element_text(size = 7), 
        legend.position = "inside", 
        legend.position.inside = c(0.75, 0.8), 
        legend.direction = "horizontal")


# Final figure ------------------------------------------------------------
lay <- matrix(c(1:8, rep(rep(9:10, each = 2), times = 2)), nrow = 4)
gp <- arrangeGrob(grobs = append(rxn_plots, list(ggplotGrob(pB), ggplotGrob(pC))), 
                  layout_matrix = lay)
ggsave("figures/final/Figure_4.pdf", gp, width = 180, height = 120, units = "mm", dpi = 600)
