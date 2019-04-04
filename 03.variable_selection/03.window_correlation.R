library(tidyverse)

windows <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  select(-BLUE) %>%
  gather(Variable, Value, -Site, -PedigreeNew) %>%
  separate(Variable, c("Variable", "Start", "End"), sep = "_", remove = TRUE) %>%
  mutate(Start = str_replace(Start, "X", "") %>% as.numeric(), 
         End = str_replace(End, "X", "") %>% as.numeric(), 
         Width = round(End - Start, 3)) %>%
  select(Site, PedigreeNew, Variable, Start, End, Width, everything())

variables <- unique(windows$Variable) %>% sort()
titles <- list("Net evapotranspiration", "Precipitation", "Solar radiation", 
            expression(T[max]), expression(T[min]))
for (i in seq_along(variables)) {
  windows %>%
    filter(Variable == variables[i], Width == 0.025) %>%
    select(-Variable, -End, -Width) %>%
    spread(Start, Value) %>%
    select(-Site, -PedigreeNew) %>%
    cor(., method = "pearson") %>%
    as_tibble(rownames = "Site1") %>%
    gather(Site2, R, -Site1) %>%
    ggplot(., aes(x = Site1, y = Site2, fill = R)) + theme_bw() +
      geom_raster() + labs(x = "", y = "", fill = "r") +
      scale_fill_distiller(type = "div", palette = "RdBu", 
                           direction = -1, limits = c(-1, 1)) +
      theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
      ggtitle(titles[[i]])
  ggsave(paste0("figures/select/", variables[i], "_0.025.pdf"), width = 10, 
         height = 8, units = "in", dpi = 300)
}
