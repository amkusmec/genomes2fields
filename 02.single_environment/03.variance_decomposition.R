### Re-run the stage one models with genotype forced to be a random effect.

library(tidyverse)
library(sommer)
library(parallel)

source("src/yield_var_comp.R")


# Identify variables present for each site --------------------------------
yield <- read_rds("data/phenotype/yield_agron0.rds")
variables <- yield %>%
  group_by(Environment, Year) %>%
  summarise(Replicate = !any(is.na(Replicate)), 
            Block = !any(is.na(Block)), 
            Row = !any(is.na(Row)), 
            Column = !any(is.na(Column)), 
            Stand = !any(is.na(Stand)), 
            StalkLodging = !any(is.na(StalkLodging)), 
            RootLodging = !any(is.na(RootLodging))) %>%
  ungroup()
variables$RootLodging[54] <- FALSE # OHH1_2015 has all 0's for RootLodging
yield <- filter(yield, Replicate != 0 | is.na(Replicate))


# Model selection for variance decomposition ------------------------------
# The main difference is the forced inclusion of `PedigreeNew` as a random--
# instead of fixed--effect.
cl <- makeCluster(32)
clusterEvalQ(cl, { library(tidyverse); library(sommer) })
clusterExport(cl, list("variables", "yield", "var_decomp"))

vcomp <- parLapply(cl, 1:nrow(variables), function(i) {
    r <- variables[i, ]
    temp <- filter(yield, Environment == r$Environment[1], Year == r$Year[1]) %>%
      mutate(Replicate = as.character(Replicate), 
             Block = paste(Replicate, Block, sep = "_"), 
             Rowf = as.character(Row), 
             Colf = as.character(Column))
    var_decomp(temp, r)
  })
# vcomp <- bind_rows(vcomp)
names(vcomp) <- with(variables, paste(Environment, Year, sep = "_"))
write_rds(vcomp, "data/phenotype/variance_components.rds")
stopCluster(cl)


# Get the unscaled variance components ------------------------------------
vcomp2 <- names(vcomp) %>%
  map_df(function(n) {
    x <- str_split(n, "_") %>% unlist()
    tibble(Year = x[2], 
           Environment = x[1], 
           Component = names(vcomp[[n]]$sigma) %>% str_remove(., "u:"), 
           Variance = unlist(vcomp[[n]]$sigma, use.names = FALSE), 
           SE = diag(vcomp[[n]]$sigmaSE))
  })


# Rename effects ----------------------------------------------------------
vcomp2 <- vcomp2 %>%
  mutate(Component = str_replace(Component, "Colf", "Column") %>%
           str_replace(., "PedigreeNew", "Genotype") %>%
           str_replace(., "Row$", "Smooth Spatial") %>%
           str_replace(., "Rowf", "Row") %>%
           str_replace(., "units", "Residual"), 
         Component = factor(Component, ordered = TRUE, 
                            levels = rev(c("Genotype", "Block", "Row", "Column", 
                                       "Smooth Spatial", "Residual"))), 
         Site = paste(Environment, Year, sep = "_")) %>%
  filter(Year != 2017)


# Plot of unscaled components ---------------------------------------------
ggplot(vcomp2, aes(x = Site, y = Variance)) + theme_classic() +
  geom_col(aes(fill = Component), colour = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(type = "qual", palette = "Set2")
ggsave("figures/single/var_unscaled.pdf", width = 11, height = 6, units = "in", dpi = 300)


# Plot of scale components ------------------------------------------------
vcomp2 %>%
  group_by(Site) %>%
  mutate(Variance = Variance/sum(Variance)) %>%
  ungroup() %>%
  ggplot(., aes(x = Site, y = Variance)) + theme_classic() +
    geom_col(aes(fill = Component), colour = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_brewer(type = "qual", palette = "Set2") +
    scale_y_continuous(labels = scales::percent)
ggsave("figures/single/var_scaled.pdf", width = 11, height = 6, units = "in", dpi = 300)


vcomp2 %>%
  filter(str_detect(Site, "NEH")) %>%
  group_by(Site) %>%
  mutate(Variance = Variance/sum(Variance)) %>%
  ungroup() %>%
  ggplot(., aes(x = Site, y = Variance)) + theme_classic() +
    geom_col(aes(fill = Component), colour = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_brewer(type = "qual", palette = "Set2") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "", y = expression(sigma[P]^2))
ggsave("figures/single/NEH_var_scaled.pdf", width = 6, height = 4, units = "in", dpi = 300)


gh2 <- map_dbl(vcomp, function(m) {
  gv <- drop(m$sigma$`u:PedigreeNew`)
  pev <- mean(diag(m$PevU$`u:PedigreeNew`$Yield))
  1 - pev/gv
})
