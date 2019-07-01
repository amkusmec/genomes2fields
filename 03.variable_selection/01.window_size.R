library(tidyverse)
library(lubridate)
library(purrrlyr)


# Identify hybrids with genomic information -------------------------------
# Only these hybrids will be included in GxE modeling
gbs <- read_rds("data/gbs/add_snps.rds")
taxa <- union(rownames(gbs$GD), rownames(gbs$GD17))


# Identify the minimum possible window size -------------------------------
# This is based on the percentage of PTT to anthesis contributed by each day
# to the earliest flowering hybrid across all sites. Counterintuitively, 
# the minimum window size is the maximum single-day contribution to PTT to 
# anthesis because this represents the minimum resolution across all sites.
yield <- read_rds("data/phenotype/yield_agron0.rds")
early <- yield %>%
  group_by(Site) %>%
  summarise(PTTA = min(PTTA)) %>%
  ungroup()

weather <- read_rds("data/weather/env_variables.rds") %>%
  mutate(Site = paste(Environment, Year, sep = "_"), 
         Date = paste(Year, Month, Day, sep = "-") %>% ymd())

min_window <- inner_join(weather, early, by = "Site") %>%
  mutate(Percent = PTT/PTTA) %>%
  filter(Percent > 0) %>%
  group_by(Site) %>%
  summarise(Percent = max(Percent)) %>%
  ungroup() %>%
  arrange(desc(Percent)) %>%
  slice(1L) %>%
  pull(Percent)


# Identify the range of developmental time --------------------------------
# This is relative to PTT to anthesis and accumulation of heat units beyond
# this threshold. The limit is set by the site where the fewest heat units
# were accumulated and the latest flowering hybrid in that site.

# Calculate PTT to harvest
yield <- yield %>%
  by_row(function(r) {
    a <- weather %>%
      filter(Environment == r$Environment, Date >= r$Planted, 
             Date <= r$Harvested) %>%
      summarise(A = sum(PTT, na.rm = TRUE)) %>%
      unlist()
    
    tibble(PTTH = a)
  }, .collate = "rows") %>%
  select(-.row)

max_window <- yield %>%
  mutate(Percent = PTTH/PTTA) %>%
  group_by(Site) %>%
  summarise(Percent = min(Percent)) %>%
  ungroup() %>%
  arrange(Percent) %>%
  slice(1L) %>%
  pull(Percent)

### Time will be scaled to 0-140% of PTT to anthesis broken up into windows
### of width 2.5%. This produces 150/2.5 = 56 windows. All possible windows from
### size 0-140% in 2.5% steps is 60*59/2 = 1540 windows per variable.
