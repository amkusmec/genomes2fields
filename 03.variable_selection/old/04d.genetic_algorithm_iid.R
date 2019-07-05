library(tidyverse)


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
  minima <- means <- medians <- numeric(maxiter)
  
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
    cat("Iter:", iter, "/", maxiter, "\r")
    if(iter == maxiter) break
    
    # Evaluate solutions
    bic <- apply(S0, 2, function(x) {
      temp <- pred[, x]
      BIC(lm(resp ~ temp))
    })
    
    # Record some population statistics
    minima[iter] <- min(bic)
    means[iter] <- mean(bic)
    medians[iter] <- median(bic)
    
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
  
  # Final statistics
  minima[iter] <- min(bic)
  minima <- minima[1:iter]
  
  means[iter] <- mean(bic)
  means <- means[1:iter]
  
  medians[iter] <- median(bic)
  medians <- medians[1:iter]
  
  list(g = colnames(pred)[S0[, which.min(bic)]], 
       minima = minima, means = means, medians = medians)
}


# Leave-one-site-out cross-validation -------------------------------------
res <- lapply(site_levels, function(s) {
  cat("\n", which(site_levels == s), "\n")
  
  # Partition into training and testing sets
  idx <- which(sites == s)
  
  test_x <- X[idx, ]
  test_y <- y[idx]
  
  train_x <- X[-idx, ]
  train_y <- y[-idx]
  
  # Optimize a linear regression model using BIC
  g <- ga(train_y, train_x, popsize = 500, n = 5, maxiter = 500, run = 100, 
          pcrossover = 0.8, pmutation = 0.2)
  
  # Compute MSE on the left out site
  beta <- coef(lm(train_y ~ train_x[, g$g]))
  pred <- drop(beta[1] + test_x[, g$g] %*% beta[-1])
  mse <- mean((pred - test_y)^2)
  
  return(list(ga = g, mse = mse))
})

write_rds(res, "data/weather/ga_select.rds")


# Genetic algorithm with all data -----------------------------------------
g <- ga(y, X, popsize = 500, n = 5, maxiter = 500, run = 100, 
        pcrossover = 0.8, pmutation = 0.2)
write_rds(g, "data/weather/ga_all.rds")


# Variables selected by leave-one-out models
variables <- lapply(res, function(x) x$ga$g) %>%
  unlist(use.names = FALSE)
length(unique(variables))
sort(table(variables))

g$g[g$g %in% variables]

# Fit of leave-one-out models on the entire dataset
r2adj <- sapply(res, function(x) {
  lm(y ~ X[, x$ga$g]) %>% broom::glance() %>% pull(adj.r.squared)
})
r2adj_all <- lm(y ~ X[, g$g]) %>% broom::glance() %>% pull(adj.r.squared)

# Composite model
vmat <- matrix(0, nrow = ncol(X), ncol = length(res))
colnames(vmat) <- site_levels
rownames(vmat) <- colnames(X)
for (i in seq_along(res)) {
  vmat[rownames(vmat) %in% res[[i]]$ga$g, i] <- 1
}

wt <- sapply(res, function(x) x$mse)
wt <- (1/(wt/sum(wt)))/sum(1/(wt/sum(wt)))
comp <- rowSums(vmat*matrix(wt, nrow = nrow(vmat), ncol = ncol(vmat), byrow = TRUE))
comp <- names(sort(comp, decreasing = TRUE))[1:5]
r2adj_comp <- lm(y ~ X[, comp]) %>% broom::glance() %>% pull(adj.r.squared)

tibble(R2 = r2adj) %>%
  ggplot(., aes(x = R2)) + theme_classic() +
    geom_histogram(binwidth = 0.01, fill = "skyblue", colour = "black") +
    geom_vline(xintercept = round(r2adj_all, 2), 
               linetype = 2, colour = "orange") +
    geom_vline(xintercept = round(r2adj_comp, 2), 
               linetype = 2, colour = "red") +
    labs(x = expression(R[adj]^2), y = "Count")
ggsave("figures/select/ga_iid_r2.pdf", width = 6, height = 4, units = "in", dpi = 300)

sum(r2adj > r2adj_all)/length(res)

best_r2 <- which.max(r2adj)
site_levels[best_r2]
res[[best_r2]]$ga$g

best_mse <- which.min(sapply(res, function(x) x$mse))
site_levels[best_mse]
res[[best_mse]]$ga$g
