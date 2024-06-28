library(tidyverse)
library(lubridate)


gbs <- read_rds("data/gbs/add_snps.rds")
taxa <- rownames(gbs$GD)

yield <- read_rds("data/phenotype/yield_agron0.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  filter(Site != "NEH3_2015") %>%
  filter(Replicate > 0 | is.na(Replicate)) %>%
  filter(PedigreeNew %in% taxa)
flowering <- yield %>%
  select(Site, Planted, Anthesis) %>%
  mutate(DTA = as.integer(Anthesis - Planted))

ggplot(flowering, aes(x = Site, y = DTA)) + theme_classic() + 
  geom_boxplot(outlier.shape = 1, outlier.colour = "red") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "", y = "Calendar Days to Anthesis") + 
  scale_y_continuous(breaks = seq(40, 100, 10))
ggsave("figures/single/phenological_variance.pdf", width = 8, height = 5, 
       units = "in", dpi = 300)
