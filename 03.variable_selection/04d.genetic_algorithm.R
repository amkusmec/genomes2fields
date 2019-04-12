library(tidyverse)
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


# Genetic algorithm function ----------------------------------------------
ga <- function(resp, pred, popsize, n = 5, maxiter = 50, run = 10, 
               pcrossover = 0.8, pmutation = 0.1) {
  # Initialize a population of solutions
  S0 <- matrix(FALSE, ncol = popsize, nrow = ncol(pred))
  S0 <- apply(S0, 2, function(x) {
    sample(c(rep(TRUE, n), rep(FALSE, ncol(pred) - n)), ncol(pred))
  })
  bic0 <- min(apply(S0, 2, function(x) {
    temp <- pred[, x]
    BIC(lm(resp ~ temp))
  }))
  
  no_change <- 0
  iter <- 0
  repeat {
    # Check for maximum iterations
    iter <- iter + 1
    if(iter == maxiter) break
    
    # Evaluate solutions
    bic <- apply(S0, 2, function(x) {
      temp <- pred[, x]
      BIC(lm(resp ~ temp))
    })
    
    # Has a better solution been produced?
    if (iter > 1) {
      if (min(bic) >= bic0) {
        no_change <- no_change + 1
      } else {
        bic0 <- min(bic)
        no_change <- 0
      }
      
      # Check for a run of no improvement of sufficient length to call convergence
      if (no_change == run) break
    }
    
    # Generate the next generation
    S1 <- matrix(FALSE, ncol = popsize, nrow = ncol(pred))
    pointer <- 2
    repeat {
      # Select two random parents
      parents <- sample(popsize, 2)
      
      # Crossover
      if (runif(1) <= pcrossover) {
        idx <- sample(which(S0[, parents[1]] | S0[, parents[2]]), n)
        S1[idx, pointer] <- TRUE
        pointer <- pointer + 1
      } else {
        if (pointer + 1 > popsize) {
          S1[, pointer] <- S0[, parents[sample(1:2, 1)]]
        } else {
          S1[, pointer] <- S0[, parents[1]]
          S1[, pointer + 1] <- S0[, parents[2]]
        }
        pointer <- pointer + 2
      }
      
      if (pointer > popsize) break
    }
    
    # Mutation
    S1 <- mapply(function(x, r) {
      if (r <= pmutation) {
        idx1 <- which(x)
        idx2 <- which(!x)
        idx1 <- idx1[sample(length(idx1), 1)]
        idx2 <- idx2[sample(length(idx2), 1)]
        x[idx1] <- !x[idx1]
        x[idx2] <- !x[idx2]
      }
      
      x
    }, x = split(t(S1[, -1]), 1:(popsize - 1)), r = runif(popsize - 1))
    
    # Save this generation plus the best solution from the previous generation
    S0 <- cbind(S0[, which.min(bic)], S1)
  }
  
  # Final evaluation
  bic <- apply(S0, 2, function(x) {
    temp <- pred[, x]
    BIC(lm(resp ~ temp))
  })
  
  colnames(pred)[S0[, which.min(bic)]]
}


# Leave-one-site-out cross validation -------------------------------------
cl <- makeCluster(length(site_levels))
clusterExport(cl, list("sites", "y", "X", "ga"))

res <- parLapply(cl, site_levels, function(s) {
  # Partition into training and testing sets
  idx <- which(sites == s)
  
  test_x <- X[idx, ]
  test_y <- y[idx]
  
  train_x <- X[-idx, ]
  train_y <- y[-idx]
  
  # Optimize a linear regression model using BIC
  g <- ga(train_y, train_x, popsize = 500, n = 5, maxiter = 50, run = 10, 
          pcrossover = 0.8, pmutation = 0.1)
  
  # Compute MSE on the left out site
  beta <- coef(lm(train_y ~ train_x[, g]))
  pred <- drop(beta[1] + test_x[, g] %*% beta[-1])
  mse <- mean((pred - test_y)^2)
  
  return(list(ga = g, mse = mse))
})

stopCluster(cl)
write_rds(res, "data/weather/ga_select.rds")


# Genetic algorithm with all data -----------------------------------------
g <- ga(y, X, popsize = 500, n = 5, maxiter = 50, run = 10, 
        pcrossover = 0.8, pmutation = 0.1)
write_rds(g, "data/weather/ga_all.rds")
