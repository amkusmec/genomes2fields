library(tidyverse)
library(parallel)


# Prepare the data --------------------------------------------------------
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$"))
ped_site <- data %>%
  select(Site, PedigreeNew) %>%
  split(., .$Site)

# Response
y <- data$BLUE

# Sites for LOO evaluation
sites <- data$Site
site_levels <- unique(sites)

# Predictors
X <- data[, -c(1:3)] %>% as.matrix()
X <- cbind(matrix(1, nrow = nrow(X), ncol = 1), X)
colnames(X)[1] <- "Intercept"

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
rm(data, d); gc()


# Transform variables -----------------------------------------------------
# Compute V^(-0.5)
e <- eigen(R)
V <- e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)

z <- drop(V %*% y)
W <- V %*% X


# Genetic algorithm function ----------------------------------------------
ga <- function(resp, pred, popsize, n = 5, maxiter = 50, run = 10, 
               pcrossover = 0.8, pmutation = 0.1) {
  minima <- means <- medians <- numeric(maxiter)
  
  # Initialize a population of solutions
  S0 <- matrix(FALSE, ncol = popsize, nrow = ncol(pred) - 1)
  S0 <- apply(S0, 2, function(x) {
    sample(c(rep(TRUE, n), rep(FALSE, ncol(pred) - n - 1)), ncol(pred) - 1)
  })
  
  # Force the intercept into the model
  S0 <- rbind(matrix(TRUE, nrow = 1, ncol = ncol(S0)), S0)
  bic0 <- min(apply(S0, 2, function(x) {
    temp <- pred[, x]
    BIC(lm(resp ~ 0 + temp))
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
      BIC(lm(resp ~ 0 + temp))
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
    S1 <- matrix(FALSE, ncol = popsize, nrow = ncol(pred) - 1)
    pointer <- 2
    repeat {
      # Select two random parents
      parents <- sample(popsize, 2)
      
      # Crossover
      if (runif(1) <= pcrossover) {
        idx <- sample(which(S0[-1, parents[1]] | S0[-1, parents[2]]), n)
        S1[idx, pointer] <- TRUE
        pointer <- pointer + 1
      } else {
        if (pointer + 1 > popsize) {
          S1[, pointer] <- S0[-1, parents[sample(1:2, 1)]]
        } else {
          S1[, pointer] <- S0[-1, parents[1]]
          S1[, pointer + 1] <- S0[-1, parents[2]]
        }
        pointer <- pointer + 2
      }
      
      if (pointer > popsize) break
    }
    S1 <- rbind(matrix(TRUE, nrow = 1, ncol = ncol(S1)), S1)
    
    # Mutation
    S1 <- mapply(function(x, r) {
      if (r <= pmutation) {
        idx1 <- which(x[-1]) + 1
        idx2 <- which(!x[-1]) + 1
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
    BIC(lm(resp ~ 0 + temp))
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


# Leave-one-site-out cross validation -------------------------------------
res <- lapply(site_levels, function(s) {
  cat("\n", which(site_levels == s), "\n")
  
  # Partition into training and testing sets
  idx <- which(sites == s)
  
  test_x <- W[idx, ]
  test_y <- z[idx]
  
  train_x <- W[-idx, ]
  train_y <- z[-idx]
  
  # Optimize a linear regression model using BIC
  g <- ga(train_y, train_x, popsize = 500, n = 5, maxiter = 500, run = 100, 
          pcrossover = 0.8, pmutation = 0.2)
  
  # Compute MSE on the left out site
  beta <- coef(lm(train_y ~ 0 + train_x[, g$g]))
  pred <- drop(test_x[, g$g] %*% beta)
  mse <- mean((pred - test_y)^2)
  
  return(list(ga = g, mse = mse))
})

write_rds(res, "data/weather/ga_select_R.rds")


# Genetic algorithm with all data -----------------------------------------
g <- ga(z, W, popsize = 500, n = 5, maxiter = 500, run = 100, 
        pcrossover = 0.8, pmutation = 0.2)
write_rds(g, "data/weather/ga_all_R.rds")



# Variables selected by leave-one-out models
variables <- lapply(res, function(x) x$ga$g) %>%
  unlist(use.names = FALSE)
variables <- variables[variables != "Intercept"]
length(unique(variables))
sort(table(variables))

# Fit of leave-one-out models on the entire dataset
r2adj <- sapply(res, function(x) {
  lm(z ~ W[, x$ga$g]) %>% broom::glance() %>% pull(adj.r.squared)
})
r2adj_all <- lm(z ~ W[, g$g]) %>% broom::glance() %>% pull(adj.r.squared)

# Composite model
vmat <- matrix(0, nrow = ncol(X) - 1, ncol = length(res))
colnames(vmat) <- site_levels
rownames(vmat) <- colnames(X)[-1]
for (i in seq_along(res)) {
  vmat[rownames(vmat) %in% res[[i]]$ga$g, i] <- 1
}

wt <- sapply(res, function(x) x$mse)
wt <- (1/(wt/sum(wt)))/sum(1/(wt/sum(wt)))
comp <- rowSums(vmat*matrix(wt, nrow = nrow(vmat), ncol = ncol(vmat), byrow = TRUE))
comp <- names(sort(comp, decreasing = TRUE))[1:5]
r2adj_comp <- lm(z ~ W[, c("Intercept", comp)]) %>% broom::glance() %>% pull(adj.r.squared)

tibble(R2 = r2adj) %>%
  ggplot(., aes(x = R2)) + theme_classic() +
    geom_histogram(binwidth = 0.001, fill = "skyblue", colour = "black") +
    geom_vline(xintercept = round(r2adj_all, 2), 
               linetype = 2, colour = "orange") +
    geom_vline(xintercept = round(r2adj_comp, 2), 
               linetype = 2, colour = "red") +
    labs(x = expression(R[adj]^2), y = "Count")
ggsave("figures/select/ga_het_r2.pdf", width = 6, height = 4, units = "in", dpi = 300)

sum(r2adj > r2adj_all)/length(res)

best_r2 <- which.max(r2adj)
site_levels[best_r2]
res[[best_r2]]$ga$g

best_mse <- which.min(sapply(res, function(x) x$mse))
site_levels[best_mse]
res[[best_mse]]$ga$g

# No relationship between MSE on the held out site and r2adj on the whole dataset
cor(r2adj, sapply(res, function(x) x$mse), method = "kendall")
plot(r2adj, sapply(res, function(x) x$mse), pch = 19)
