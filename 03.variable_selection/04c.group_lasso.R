library(tidyverse)
library(SGL)


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

# Create an index of variable groups
group_idx <- colnames(X) %>%
  str_split(., "_") %>%
  sapply(., function(x) x[1]) %>%
  factor() %>% as.integer()

### For linear regression, SGL automatically centers the response, so no
### intercept needs to be added to the design matrix.



# SGL with leave-one-site-out CV ------------------------------------------
res <- lapply(site_levels, function(s) {
  cat(which(site_levels == s), "\n")
  
  # Partition into training and testing sets
  idx <- which(sites == s)
  
  test_x <- X[idx, ]
  test_y <- y[idx]
  
  train_x <- X[-idx, ]
  train_y <- y[-idx]
  
  # Estimate the sparse group LASSO model with 20-fold cross-validation
  sgl <- cvSGL(list(x = train_x, y = train_y), group_idx, nfold = 20)
  
  # `cvSGL` does not return a suitable SGL object for prediction, so create a 
  # fake object using information from the cv.SGL object and the procedure used
  # to calculate standardization variables: 
  # (https://github.com/cran/SGL/blob/master/R/SGLmain.R)
  sgl_fake <- list(beta = sgl$fit$beta, 
                   lambdas = sgl$fit$lambdas, 
                   type = "linear", 
                   intercept = mean(train_y))
  X.means <- apply(train_x, 2, mean)
  X.scale <- sapply(1:ncol(train_x), 
                    function(i) sqrt(sum((train_x[, i] - X.means[i])^2)))
  names(X.scale) <- names(X.means)
  sgl_fake$X.transform <- list(X.scale = X.scale, X.means = X.means)
  attr(sgl_fake, "class") <- "SGL"
  
  # Estimate the optimal lambda
  lambda_min <- which.min(sgl$lldiff)
  upper_min <- sgl$lldiff[lambda_min] + sgl$llSD[lambda_min]
  lambda_opt <- max(which((sgl$lldiff[1:lambda_min] - sgl$llSD[1:lambda_min]) > upper_min)) + 1
  
  # Compute MSE on the left out site
  pred <- predictSGL(sgl_fake, test_x, lambda_opt) %>% drop()
  mse <- mean((test_y - pred)^2)
  
  list(sgl = sgl, lambda_opt = lambda_opt, mse = mse)
})

write_rds(res, "data/weather/sgl_select.rds")
