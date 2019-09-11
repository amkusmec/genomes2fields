library(tidyverse)
library(lubridate)
library(purrrlyr)
library(parallel)


# Prepare yield and weather data ------------------------------------------
weather <- read_rds("data/weather/env_variables.rds") %>%
  mutate(Site = paste(Environment, Year, sep = "_"), 
         Date = paste(Year, Month, Day, sep = "-") %>% ymd())

# Hybrids within sites are summarized by the mean CHU to anthesis and planting
# dates. Because there are never more than 2 replicates of a hybrid within a 
# site, the mean and median are identical.
yield <- read_rds("data/phenotype/yield_agron0.rds") %>%
  select(Site, PedigreeNew, CHUA, Planted, Harvested) %>%
  group_by(Site, PedigreeNew) %>%
  summarise(CHUA = mean(CHUA), 
            Planted = mean(Planted), 
            Harvested = mean(Harvested)) %>%
  ungroup()


# Find the dates to define windows ----------------------------------------
# Date on which a hybrid accumulates CHU equal to 2.5% increments of 150%
# CHU to anthesis
thresholds <- seq(0.025, 1.5, 0.025)
dates <- yield %>%
  by_row(function(r) {
    temp <- weather %>%
      filter(Site == r$Site, Date >= r$Planted, Date <= r$Harvested)
    chu <- cumsum(temp$CHU)/r$CHUA
    indices <- map_int(thresholds, function(x) {
      as.integer(min(which(chu - x >= 0)))
    })
    ymd(temp$Date[indices])
  }, .collate = "cols", .labels = FALSE) %>%
  mutate_all(funs(as_date))
names(dates) <- make.names(thresholds)

dates <- dates %>%
  mutate(X0 = as.character(yield$Planted)) %>%
  select(X0, everything())

date_idx <- matrix(0L, nrow = nrow(dates), ncol = ncol(dates))
for (i in 1:ncol(date_idx)) {
  for (j in 1:nrow(date_idx)) {
    date_idx[j, i] <- which(weather$Site == yield$Site[j] & 
                              weather$Date == dates[[i]][j])
  }
}


# Compute mean/summed environmental variables -----------------------------
# Temperatures are averaged over a window; all other variables are summed
cl <- makeCluster(5)
clusterExport(cl, list("dates", "date_idx", "weather"))

res <- parLapply(cl, c("TMIN", "TMAX", "PPT", "SR", "NET"), function(v) {
  temp <- matrix(0, nrow = nrow(dates), ncol = ncol(dates)*(ncol(dates) - 1)/2)
  colnames(temp) <- make.names(1:ncol(temp))
  counter <- 1
  for (i in 1:(ncol(dates) - 1)) {
    for (j in (i + 1):ncol(dates)) {
      temp[, counter] <- apply(date_idx[, c(i, j)], 1, function(m) {
        if (i == 1) {
          if (v %in% c("TMIN", "TMAX")) {
            mean(weather[[v]][m[1]:m[2]])
          } else {
            sum(weather[[v]][m[1]:m[2]])
          }
        } else {
          if (v %in% c("TMIN", "TMAX")) {
            mean(weather[[v]][(m[1] + 1):m[2]])
          } else {
            sum(weather[[v]][(m[1] + 1):m[2]])
          }
        }
      })
      
      colnames(temp)[counter] <- paste(v, colnames(dates)[i], colnames(dates)[j], sep = "_")
      counter <- counter + 1
    }
  }
  
  temp
})

stopCluster(cl)

res <- do.call("cbind", res)


# Combine environmental variables and yield data --------------------------
blue <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
sites <- rep(names(blue), times = sapply(blue, function(x) nrow(x$blue)))
blue <- blue %>%
  map_df(function(x) x$blue) %>%
  mutate(Site = sites)

yield <- cbind(yield, res) %>%
  inner_join(., blue, by = c("Site", "PedigreeNew")) %>%
  select(-CHUA, -Planted, -Harvested) %>%
  select(Site, PedigreeNew, BLUE, everything())


# Identify hybrids with genomic information -------------------------------
# Only these hybrids will be included in GxE modeling
gbs <- read_rds("data/gbs/add_snps.rds")
taxa <- union(rownames(gbs$GD), rownames(gbs$GD17))
yield <- filter(yield, PedigreeNew %in% taxa)

write_rds(yield, "data/phenotype/yield_blue_env.rds")
