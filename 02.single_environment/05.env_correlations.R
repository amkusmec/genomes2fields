library(tidyverse)
library(lubridate)
library(purrrlyr)
library(tensorBSS)


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
# write_rds(etens, "data/weather/env_tensor.rds")

### Can be run in the background with data/weather/env_tensor_decomp.R
# tica <- tgJADE(etens, maxiter = 500)
# write_rds(tica, "data/weather/tgJADE.rds")
tica <- read_rds("data/weather/tgJADE.rds")


# Hierarchical clustering using Euclidean distance ------------------------
# Center and standardize the data tensor for each combination of
# window and variable over sites. Using apply is more elegant, but this
# preserves the dimensions and dimension names.
for (i in 1:(dim(etens)[1])) {
  for (k in 1:(dim(etens)[3])) {
    etens[i, , k] <- scale(etens[i, , k], center = TRUE, scale = TRUE) %>%
      as.vector()
  }
}

# Euclidean distance matrix between sites using all variables
all_dist <- matrix(0, nrow = length(sites), ncol = length(sites))
dimnames(all_dist) <- list(sites, sites)
for (i in 1:(length(sites) - 1)) {
  for (j in (i + 1):length(sites)) {
    all_dist[j, i] <- all_dist[i, j] <- sqrt(sum((etens[, i, ] - etens[, j, ])^2))
  }
}

all_clust <- hclust(as.dist(all_dist), method = "ward.D")

all_yr_labs <- all_clust$labels %>% str_replace(., "[A-Z]{2}H[0-9]_", "") %>%
  factor() %>% as.integer()
all_yr_cols <- WGCNA::labels2colors(all_yr_labs)
all_st_labs <- all_clust$labels %>% str_replace(., "H[1-4]_201[4-7]", "") %>%
  factor() %>% as.integer()
all_st_cols <- WGCNA::labels2colors(all_st_labs)

pdf("figures/single/dendro_all.pdf", width = 11, height = 6)
WGCNA::plotDendroAndColors(all_clust, cbind(all_yr_cols, all_st_cols), 
                           c("Year", "State"), dendroLabels = NULL, 
                           addGuide = TRUE, main = "All variables")
dev.off()

# Euclidean distance matrices between sites for individual variables
variables <- names(windowed)[-c(1, 2)]
single <- list()
for (k in variables) {
  d <- matrix(0, nrow = length(sites), ncol = length(sites), 
              dimnames = list(sites, sites))
  for (i in 1:(length(sites) - 1)) {
    for (j in (i + 1):length(sites)) {
      d[j, i] <- d[i, j] <- sqrt(sum((etens[k, i, ] - etens[k, j, ])^2))
    }
  }
  
  clust <- hclust(as.dist(d), method = "ward.D")
  yr_labs <- clust$labels %>% str_replace(., "[A-Z]{2}H[0-9]_", "") %>%
    factor() %>% as.integer()
  yr_cols <- WGCNA::labels2colors(yr_labs)
  st_labs <- clust$labels %>% str_replace(., "H[1-4]_201[4-7]", "") %>%
    factor() %>% as.integer()
  st_cols <- WGCNA::labels2colors(st_labs)
  
  single[[k]] <- list(dist = d, clust = clust, yr_labs = yr_labs, 
                      yr_cols = yr_cols, st_labs = st_labs, st_cols = st_cols)
}

variable_names <- c("Minimum temperature", "Maximum temperature", 
                    "Precipitation", "Solar radiation", 
                    "Net evapotranspiration")
for (i in seq_along(variables)) {
  pdf(paste0("figures/single/dendro_", variables[i], ".pdf"), 
      width = 11, height = 6)
  WGCNA::plotDendroAndColors(single[[i]]$clust, cbind(single[[i]]$yr_cols, single[[i]]$st_cols), 
                             c("Year", "State"), dendroLabels = NULL, addGuide = TRUE, 
                             main = variable_names[i])
  dev.off()
}
