library(tidyverse)
library(lubridate)
library(purrrlyr)


# Calculate median CHU to anthesis for each site --------------------------
chua <- read_rds("data/phenotype/yield_agron0.rds") %>%
  group_by(Site) %>%
  summarise(CHUA = median(CHUA)) %>%
  ungroup()


# Identify 2.5% CHUA windows from 0-150% ----------------------------------
weather <- read_rds("data/weather/env_variables.rds") %>%
  mutate(Site = paste(Environment, Year, sep = "_"), 
         Date = paste(Year, Month, Day, sep = "-") %>% ymd()) %>%
  inner_join(., chua, by = "Site") %>%
  group_by(Site) %>%
  mutate(CHU = cumsum(CHU)/CHUA) %>%
  ungroup() %>%
  filter(CHU <= 1.5)

thresholds <- seq(0.025, 1.5, 0.025)
weather <- weather %>%
  by_row(function(r) {
    min(which(r$CHU[1] - thresholds <= 0))
  }) %>%
  unnest(.out) %>%
  rename(Window = .out)


# Summarise variables by window -------------------------------------------
windowed <- weather %>%
  group_by(Site, Window) %>%
  summarise_at(c("TMIN", "TMAX", "PPT", "SR", "NET"), mean) %>%
  ungroup()


# Tensorial independent component analysis --------------------------------
# Create the data tensor
sites <- unique(windowed$Site)
etens <- array(0, c(ncol(windowed) - 2, length(sites), length(thresholds)))
dimnames(etens) <- list(names(windowed)[-c(1, 2)], sites, make.names(thresholds))
for (i in 1:(ncol(windowed) - 2)) {
  for (j in seq_along(sites)) {
    etens[i, j, ] <- windowed[[i + 2]][windowed$Site == sites[j]]
  }
}

# Create a list of data matrices
elist <- list()
for (i in 1:(ncol(windowed) - 2)) {
  temp <- names(windowed)[i + 2]
  names(windowed)[i + 2] <- "X"
  
  elist[[i]] <- windowed %>%
    select(Site, Window, X) %>%
    spread(Site, X) %>%
    select(-Window) %>%
    as.matrix()
  rownames(elist[[i]]) <- make.names(thresholds)
  
  names(windowed)[i + 2] <- temp
}
names(elist) <- names(windowed)[-c(1, 2)]



### Distance based (Euclidean) clustering
### Standardize data tensor for (window x variable) over site
### Compute Euclidean distance between sites as
### sqrt(sum(window x variable difference squared))
### hierarchical clustering

### do clustering for each variable separately as well

### plot all results using dendrograms and year/state color labels as for
### clusterings based on yield correlations
