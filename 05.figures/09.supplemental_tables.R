library(tidyverse)

data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  filter(Site != "NEH3_2015") %>%
  separate(Site, c("Environment", "Year"), sep = "_", remove = FALSE) %>%
  select(BLUE, PedigreeNew, Site, Environment, Year)


ped_count <- data %>%
  count(PedigreeNew) %>%
  arrange(PedigreeNew) %>%
  rename(Hybrid = PedigreeNew, Environments = n)
write_csv(ped_count, "data/tableS1_hybrids.csv")

meta <- read_csv("data/metadata/metadata_clean.csv")
env_count <- data %>%
  count(Year, Environment) %>%
  mutate(Year = as.integer(Year)) %>%
  inner_join(., meta, by = c("Year", "Environment")) %>%
  select(Year, Environment, Latitude, Longitude, n) %>%
  rename(Location = Environment, Hybrids = n)
write_csv(env_count, "data/tableS2_environments.csv")
