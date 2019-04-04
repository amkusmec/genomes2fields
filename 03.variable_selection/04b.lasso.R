library(tidyverse)
library(glmnet)
library(parallel)


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


# LASSO with leave-one-site-out CV ----------------------------------------
cl <- makeCluster(length(site_levels))
clusterEvalQ(cl, library(glmnet))
clusterExport(cl, list("sites", "y", "X"))

res <- parLapply(cl, site_levels, function(s) {
  # Partition into training and testing sets
  idx <- which(sites == s)
  
  test_x <- X[idx, ]
  test_y <- y[idx]
  
  train_x <- X[-idx, ]
  train_y <- y[-idx]
  
  # Train the LASSO model
  lasso <- cv.glmnet(train_x, train_y, alpha = 1, nfolds = 20)
  
  # Compute MSE on the left out site
  pred <- predict(lasso, newx = test_x, type = "response", s = lasso$lambda.1se)
  mse <- mean((test_y - drop(pred))^2)
  
  return(list(lasso = lasso, mse = mse))
})

stopCluster(cl)
write_rds(res, "data/weather/lasso_select.rds")


# LASSO with all data -----------------------------------------------------
lasso <- cv.glmnet(X, y, alpha = 1, nfolds = 20)
write_rds(lasso, "data/weather/lasso_all.rds")
