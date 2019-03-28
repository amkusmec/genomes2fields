library(tidyverse)
library(sommer)

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
vcomp <- 1:nrow(variables) %>%
  map_df(function(i) {
    r <- variables[i, ]
    temp <- filter(yield, Environment == r$Environment[1], Year == r$Year[1]) %>%
      mutate(Replicate = as.character(Replicate), 
             Block = paste(Replicate, Block, sep = "_"), 
             Rowf = as.character(Row), 
             Colf = as.character(Column))
    cat(temp$Site[1], "\n")
    var_decomp(temp, r)
  })


# Rename effects ----------------------------------------------------------
vcomp <- vcomp %>%
  mutate(Component = str_replace(Component, "Colf", "Column") %>%
           str_replace(., "PedigreeNew", "Genotype") %>%
           str_replace(., "Row$", "Smooth Spatial") %>%
           str_replace(., "Rowf", "Row") %>%
           str_replace(., "units", "Residual"), 
         Component = factor(Component, ordered = TRUE, 
                            levels = c("Genotype", "Block", "Row", "Column", 
                                       "Smooth Spatial", "Residual")), 
         Site = paste(Environment, Year, sep = "_"))


# Plot of unscaled components ---------------------------------------------
ggplot(vcomp, aes(x = Site, y = Variance)) + theme_classic() +
  geom_col(aes(fill = Component), colour = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(type = "qual", palette = "Set2")
ggsave("figures/single/var_unscaled.pdf", width = 11, height = 6, units = "in", dpi = 300)


# Plot of scale components ------------------------------------------------
vcomp %>%
  group_by(Site) %>%
  mutate(Variance = Variance/sum(Variance)) %>%
  ungroup() %>%
  ggplot(., aes(x = Site, y = Variance)) + theme_classic() +
    geom_col(aes(fill = Component), colour = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_brewer(type = "qual", palette = "Set2") +
    scale_y_continuous(labels = scales::percent)
ggsave("figures/single/var_scaled.pdf", width = 11, height = 6, units = "in", dpi = 300)
