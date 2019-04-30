### Similar to SGL, SSR has desirable statistical properties. There are 
### computational time issues similar to SSR. Additionally, there is some odd
### internal parallelization that is going on. Individual processes take up to
### ~75% RAM (on a 384 GB machine!). There is also an error in returning results
### in a parallelized setting, but I haven't been able to tell if that is an 
### issue with the package or my code. The grouped version would also be off 
### interest if computationally tractable.

library(tidyverse)
library(spikeslab)

# Prepare the data --------------------------------------------------------
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$"))

# Response
y <- data$BLUE

# Sites for LOO evaluation
sites <- data$Site
site_levels <- unique(sites)

# Predictors
X <- data[, -c(1:3)] %>% as.matrix()


# SSVS with leave-one-site-out CV -----------------------------------------
res <- lapply(site_levels, function(s) {
  cat(which(site_levels == s), "\n")
  
  # Partition into training and testing sets
  idx <- which(sites == s)
  
  test_x <- X[idx, ]
  test_y <- y[idx]
  
  train_x <- X[-idx, ]
  train_y <- y[-idx]
  
  # Train the model
  ssvs <- cv.spikeslab(x = train_x, y = train_y, K = 20, plot.it = FALSE, 
                       n.iter1 = 2000, n.iter2 = 5000, mse = TRUE, 
                       bigp.smalln = TRUE, max.var = 5, parallel = 5)
  
  # Compute MSE on the left out site
  pred <- predict(ssvs, test_x)
  mse <- mean((pred$yhat.gnet - test_y)^2)
  
  return(list(ssvs = ssvs, mse = mse))
})

write_rds(res, "data/weather/ssvs_select.rds")
