library(tidyverse)
library(parallel)


# Prepare the data --------------------------------------------------------
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  separate(Site, c("Environment", "Year"), sep = "_", remove = FALSE) %>%
  select(BLUE, PedigreeNew, Site, Environment, Year, everything())
ped_site <- data %>%
  select(Site, PedigreeNew) %>%
  split(., .$Site)

# Response
y <- data$BLUE
y <- y - mean(y)

# Predictors
X <- cbind(model.matrix(~ 0 + Environment + Year, data = data), 
           as.matrix(data[, -c(1:5)]))
X[, -c(1:30)] <- scale(X[, -c(1:30)], center = TRUE, scale = FALSE)

# Function to construct the direct sum of a list of matrices
source("src/direct_sum.R")

# Variance-covariance matrix of the BLUEs
d <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
d <- d[!str_detect(names(d), "2017$")]
d <- mapply(function(x, y) {
  idx <- which(rownames(x$vcov) %in% y$PedigreeNew)
  x$vcov[idx, idx]
}, x = d, y = ped_site)
R <- direct_sum(d)
wts <- 1/diag(R)


# Genetic algorithm function ----------------------------------------------
ga <- function(resp, pred, wts, popsize, n = 5, maxiter = 50, run = 10, 
               pcrossover = 0.8, pmutation = 0.1, gamma = 1, verbose = TRUE) {
  # Structures for saving information
  minima <- means <- medians <- numeric(maxiter)
  
  # Initialize a population of solutions
  nvar <- sum(stringr::str_detect(colnames(pred), "_"))
  offset <- ncol(pred) - nvar
  S0 <- matrix(FALSE, ncol = popsize, nrow = nvar)
  S0[, -1] <- apply(S0[, -1], 2, function(x) {
    m <- sample(1:n, 1)
    sample(c(rep(TRUE, m), rep(FALSE, nvar - m)), nvar)
  })
  
  # Identify the best starting model
  bic <- apply(S0, 2, function(x) {
    temp <- pred[, c(rep(TRUE, offset), x)]
    BIC(lm(resp ~ 0 + temp, weights = wts)) + 2*gamma*log(choose(nvar, sum(x)))
  })
  bic0 <- min(bic, na.rm = TRUE)
  
  no_change <- 0
  iter <- 0
  repeat {
    # Check for maximum iterations
    iter <- iter + 1
    if (verbose) cat("Iter:", iter, "/", maxiter, "; Dimension:", sum(S0[, which.min(bic)]), "\r")
    if(iter == maxiter) break
    
    # Evaluate solutions
    bic <- apply(S0, 2, function(x) {
      temp <- pred[, c(rep(TRUE, offset), x)]
      BIC(lm(resp ~ 0 + temp, weights = wts)) + 2*gamma*log(choose(nvar, sum(x)))
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
    S1 <- matrix(FALSE, ncol = popsize, nrow = nvar)
    pointer <- 2
    repeat {
      # Select two random parents
      parents <- sample(popsize, 2)
      
      # Crossover
      if (runif(1) <= pcrossover) {
        cross_point <- sample(2:(nvar - 1), 1)
        idx <- c(which(S0[1:cross_point, parents[1]]), 
                 which(S0[(cross_point + 1):nvar, parents[2]]))
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
        idx <- sample(nvar, 1)
        x[idx] <- !x[idx]
      }
      
      x
    }, x = split(t(S1[, -1]), 1:(popsize - 1)), r = runif(popsize - 1))
    
    # Save this generation plus the best solution from the previous generation
    S0 <- cbind(matrix(S0[, which.min(bic)], ncol = 1), S1)
  }
  
  # Final evaluation
  bic <- apply(S0, 2, function(x) {
    temp <- pred[, c(rep(TRUE, offset), x)]
    BIC(lm(resp ~ 0 + temp, weights = wts)) + 2*gamma*log(choose(nvar, sum(x)))
  })
  
  # Final statistics
  minima[iter] <- min(bic)
  minima <- minima[1:iter]
  
  means[iter] <- mean(bic)
  means <- means[1:iter]
  
  medians[iter] <- median(bic)
  medians <- medians[1:iter]
  
  list(g = colnames(pred)[c(rep(FALSE, offset), S0[, which.min(bic)])], 
       minima = minima, means = means, medians = medians)
}


# Genetic algorithm with all data -----------------------------------------
# We force fixed environment and year effects into the model and control for
# heteroscedasticity by using the inverse variance of the hybrid BLUEs.
# Model fitness is evaluated using the extended BIC to account for the 
# dimensionality of the model space.
seeds <- c(615151, 317119, 379563, 885645, 222345, 373540, 812586, 730468, 
           492480, 877399, 55276, 217964, 399946, 768559, 941469, 250222, 292795,
           467124, 780664, 183441, 963089, 355954, 825231, 130461, 554633, 
           682924, 561486, 331539, 946933, 487918, 740381, 99897, 922192,
           92956, 191295, 996894, 710191, 527575, 842501, 355097, 797741, 210842,
           410334, 444784, 769311, 217768, 794136, 586385, 197430, 440136, 360728, 
           740153, 754293, 958298, 189606, 900939, 741204, 884192, 701602, 937754,
           5838, 323438, 117217, 958705, 764388, 805319, 57676, 196027, 636488,
           72912, 820377, 744595, 298550, 221873, 660686, 301155, 711609, 746138,
           422521, 643464, 460942, 257990, 843141, 519712, 254114, 443789, 273042,
           136967, 747637, 139247, 857166, 10271, 120128, 963215, 605301, 93668, 
           692502, 84266, 513859, 764728)
cl <- makeCluster(length(seeds)/2)
clusterExport(cl, list("y", "X", "wts", "ga"))
g <- parLapply(cl, seeds, function(s) {
  set.seed(s)
  ga(y, X, wts, popsize = 200, n = 5, maxiter = 1000, run = 200, pcrossover = 0.8, 
     pmutation = 0.2, gamma = 1, verbose = FALSE)
})
write_rds(g, "data/weather/ga_het_resid.rds")
stopCluster(cl)


# Reduce the number of variables ------------------------------------------
# Get all the variables that were selected in any model
which.min(sapply(g, function(x) min(x$minima)))
vars <- lapply(g, function(x) x$g) %>% unlist(use.names = FALSE)
tt <- table(vars)

# For variables that were selected by at least 2 models...
frequent <- tibble(Variable = names(tt)[tt > 1]) %>%
  separate(Variable, c("Category", "Start", "End"), sep = "_", remove = FALSE) %>%
  mutate_at(c("Start", "End"), str_remove, pattern = "X") %>%
  mutate_at(c("Start", "End"), as.numeric) %>%
  split(., .$Category) %>%
  map_df(function(df) {
    # Construct a distance matrix
    X_sub <- X[, df$Variable]
    X_dist <- dist(t(X_sub))
    
    # Perform a fast complete linkage clustering
    X_h <- hclust(X_dist, method = "complete")
    
    # Seed a PAM clustering with the complete linkage results to choose an
    # optimal number of clusters
    X_kk <- dendextend::find_k(X_h)
    
    tibble(Variable = names(X_kk$pamobject$clustering), 
           Cluster = X_kk$pamobject$clustering)
  }) %>%
  separate(Variable, c("Category", "Start", "End"), sep = "_", remove = FALSE) %>%
  mutate_at(c("Start", "End"), str_remove, pattern = "X") %>%
  mutate_at(c("Start", "End"), as.numeric)

# Combine overlapping ranges within a cluster
# Code from: 
# https://stackoverflow.com/questions/53213418/r-collapse-and-merge-overlapping-time-intervals
ranges <- frequent %>%
  group_by(Category, Cluster) %>%
  arrange(Start) %>%
  mutate(indx = c(0, cumsum(as.numeric(lead(Start)) > 
                              cummax(as.numeric(End)))[-n()])) %>%
  group_by(Category, Cluster, indx) %>%
  summarise(Start = min(Start), 
            End = max(End)) %>%
  ungroup() %>%
  select(-indx) %>%
  distinct(Category, Start, End, .keep_all = TRUE)

write_csv(ranges, "data/weather/ga_windows.csv")
