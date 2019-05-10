library(tidyverse)
library(lubridate)
library(purrrlyr)


# Prepare yield and weather data ------------------------------------------
weather <- read_rds("data/weather/env_variables.rds") %>%
  mutate(Site = paste(Environment, Year, sep = "_"), 
         Date = paste(Year, Month, Day, sep = "-") %>% ymd())

yield <- read_rds("data/phenotype/yield_agron0_filtered.rds") %>%
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


# Prepare the fused windows -----------------------------------------------
fused <- read_rds("data/weather/fused_windows.rds") %>% 
  group_by(Beta) %>%
  summarise(Variable = unique(Variable), 
            Start = min(Start), 
            End = max(End)) %>%
  ungroup() %>%
  select(-Beta) %>%
  arrange(Variable, Start) %>%
  split(., .$Variable)


# Compute fused environmental variables -----------------------------------
res <- lapply(fused, function(df) {
  temp <- matrix(0, nrow = nrow(dates), ncol = nrow(df))
  colnames(temp) <- paste(df$Variable, make.names(df$Start), 
                          make.names(df$End), sep = "_")
  v <- unique(df$Variable)
  for (i in 1:nrow(df)) {
    idx <- as.integer(c(df$Start[i], df$End[i])/0.025 + 1)
    temp[, i] <- apply(date_idx[, idx], 1, function(m) mean(weather[[v]][m[1]:m[2]]))
  }
  
  temp
})

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

write_rds(yield, "data/phenotype/yield_blue_env_fused.rds")


# Prepare the data --------------------------------------------------------
yield <- filter(yield, !str_detect(Site, "2017$"))
ped_site <- yield %>%
  select(Site, PedigreeNew) %>%
  split(., .$Site)

# Response
y <- yield$BLUE

# Sites for LOO evaluation
sites <- yield$Site
site_levels <- unique(sites)

# Predictors
X <- yield[, -c(1:3)] %>% as.matrix()

# Function to construct the direct sum of a list of matrices
direct_sum <- function(matrices) {
  n <- sapply(matrices, nrow) %>% sum()
  m <- sapply(matrices, ncol) %>% sum()
  
  temp <- matrix(0, nrow = n, ncol = m)
  for (i in seq_along(matrices)) {
    if (i == 1) {
      start_r <- 0
      start_c <- 0
    }
    else {
      start_r <- sapply(matrices[1:(i - 1)], nrow) %>% sum()
      start_c <- sapply(matrices[1:(i - 1)], ncol) %>% sum()
    }
    
    temp[(start_r + 1):(start_r + nrow(matrices[[i]])),
         (start_c + 1):(start_c + ncol(matrices[[i]]))] <- matrices[[i]]
  }
  
  return(temp)
}

# Variance-covariance matrix of the BLUEs
d <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
d <- d[!str_detect(names(d), "2017$")]
d <- mapply(function(x, y) {
  idx <- which(rownames(x$vcov) %in% y$PedigreeNew)
  x$vcov[idx, idx]
}, x = d, y = ped_site)
R <- direct_sum(d)

# Clean-up
rm(d); gc()

# Compute V^(-0.5)
e <- eigen(R)
V <- e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)

z <- drop(V %*% y)
W <- V %*% X

# Center (and scale)
y_sc <- z - mean(z)
X_sc <- scale(W, center = TRUE, scale = TRUE)


# Exhaustive search of the model space ------------------------------------
combos <- combn(colnames(res), 5)

### This creates an enormous (~452 GB) object in memory. There are more elegant
### ways to do this, but none that are worth my time right now.
models <- lapply(1:ncol(combos), function(x) {
  cat("Model:", x, "/", ncol(combos), "\r")
  temp <- as_tibble(cbind(y_sc, X_sc[, combos[, x]]))
  full_form <- as.formula(paste0("y_sc ~ 0 + ", paste(combos[, x], collapse = "*")))
  step(lm(full_form, data = temp), direction = "backward", trace = 0, 
       k = log(nrow(temp)))
})
bic <- sapply(models, BIC)
write_rds(bic, "data/weather/fused_bic.rds")
rm(models); gc()

percent1 <- quantile(bic, probs = 0.01)
idx <- which(bic <= percent1)

temp <- as_tibble(cbind(y_sc, X_sc[, combos[, idx[2]]]))
full_form <- as.formula(paste0("y_sc ~ 0 + ", paste(combos[, idx[2]], collapse = "*")))
full_model <- lm(full_form, data = temp)
null_model <- lm(y_sc ~ 0, data = temp)
m <- step(full_model, scope = list(lower = null_model), direction = "backward", k = log(nrow(temp)))
View(combos[, idx])
