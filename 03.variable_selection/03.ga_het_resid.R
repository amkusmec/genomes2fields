### Implements a genetic algorithm to select environmental variables using 
### weighted least squares and the extended BIC.

library(tidyverse)
library(parallel)


# Prepare the data --------------------------------------------------------
data <- read_rds("data/phenotype/yield_blue_env.rds") %>%
  filter(!str_detect(Site, "2017$")) %>%
  filter(Site != "NEH3_2015") %>%
  separate(Site, c("Environment", "Year"), sep = "_", remove = FALSE) %>%
  select(BLUE, PedigreeNew, Site, Environment, Year, everything())
net <- names(data)[which(str_detect(names(data), "NET"))]
data <- mutate_at(data, net, function(x) -1*x)
ped_site <- data %>%
  select(Site, PedigreeNew) %>%
  split(., .$Site)

# Response
y <- data$BLUE

# Predictors
X <- cbind(model.matrix(~ 1 + Site, data = data), 
           as.matrix(data[, -c(1:5)]))
X[, -c(1:45)] <- scale(X[, -c(1:45)], center = TRUE, scale = TRUE)

# Function to construct the direct sum of a list of matrices
# source("src/direct_sum.R")

# Variance-covariance matrix of the BLUEs
d <- read_rds("data/phenotype/yield_stage_one_all_agron0.rds")
d <- d[names(d) %in% names(ped_site)]
d <- mapply(function(x, y) {
  idx <- which(rownames(x$vcov) %in% y$PedigreeNew)
  temp <- x$vcov[idx, idx]
  c(temp[1, 1], temp[1, 1] + diag(temp)[-1] + 2*temp[-1, 1])
}, x = d, y = ped_site)
# R <- direct_sum(d)
wts <- 1/unlist(d, use.names = FALSE)

# # Transform the response and predictors to meet the assumption of iid residuals
# cm <- chol(R)
# Si <- solve(cm)
# z <- drop(Si %*% matrix(y, ncol = 1))
# W <- Si %*% X


# Genetic algorithm function ----------------------------------------------
ga <- function(resp, pred, wts, popsize, n = 5, maxiter = 50, run = 10, 
               pcrossover = 0.8, pmutation = 0.1, elitism = 0.05, gamma = 1, 
               verbose = TRUE) {
  # Structures for saving information
  minima <- means <- medians <- numeric(maxiter)
  
  # Initialize a population of solutions
  nvar <- sum(stringr::str_count(colnames(pred), "_") == 2)
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
    S0 <- S0[, order(bic, decreasing = FALSE)]
    S1 <- matrix(FALSE, ncol = popsize*(1 - elitism), nrow = nvar)
    pointer <- 1
    repeat {
      # Select two random parents
      parents <- sample(1:(popsize*elitism), 2)
      
      # Crossover
      if (runif(1) <= pcrossover) {
        cross_point <- sample(2:(nvar - 1), 1)
        idx <- c(which(S0[1:cross_point, parents[1]]), 
                 which(S0[(cross_point + 1):nvar, parents[2]]))
        S1[idx, pointer] <- TRUE
        pointer <- pointer + 1
      } else {
        if (pointer + 1 > ncol(S1)) {
          S1[, pointer] <- S0[, parents[sample(1:2, 1)]]
          pointer <- pointer + 1
        } else {
          S1[, pointer] <- S0[, parents[1]]
          S1[, pointer + 1] <- S0[, parents[2]]
          pointer <- pointer + 2
        }
      }
      
      if (pointer > ncol(S1)) break
    }
    
    # Mutation
    S1 <- mapply(function(x, r) {
      if (r <= pmutation) {
        idx <- sample(nvar, 1)
        x[idx] <- !x[idx]
      }
      
      x
    }, x = split(t(S1), 1:ncol(S1)), r = runif(ncol(S1)))
    
    # Save this generation plus the best solutions from the previous generation
    S0 <- cbind(S0[, 1:(200*0.05)], S1)
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
  ga(y, X, wts, popsize = 800, n = 5, maxiter = 1000, run = 200, pcrossover = 0.8, 
     pmutation = 0.5, elitism = 0.05, gamma = 1, verbose = FALSE)
})
write_rds(g, "data/weather/ga_het_resid_mean.rds")
stopCluster(cl)


# Reduce the number of variables ------------------------------------------
# Get all the variables that were selected in any model
nvars <- sapply(g, function(x) length(x$g))
vars <- lapply(g, function(x) x$g) %>% unlist(use.names = FALSE)
minima <- sapply(g, function(x) min(x$minima))
niters <- sapply(g, function(x) length(x$minima))

tibble(NVars = nvars, Minima = minima, NIters = niters) %>%
  gather(Term, Value, everything()) %>%
  ggplot(., aes(x = Value)) + theme_classic() +
    geom_histogram() + facet_wrap(~ Term, scales = "free")

tt <- table(vars)
enframe(tt) %>% mutate(value = as.integer(value)) %>% pull(value) %>% summary()
var_table <- tibble(Variable = names(tt)) %>%
  separate(Variable, c("EVar", "Start", "End"), sep = "_", remove = FALSE) %>%
  mutate(Start = str_remove(Start, "X") %>% as.numeric(),
         End = str_remove(End, "X") %>% as.numeric()) %>%
  group_by(EVar) %>%
  arrange(Start, End) %>%
  mutate(Index = 1:n()) %>%
  ungroup()
pB <- ggplot(var_table, aes(y = Index)) + theme_bw() +
  geom_segment(aes(yend = Index, x = Start, xend = End)) +
  facet_wrap(~ EVar, ncol = 2) + labs(x = "% CHU to Anthesis", y = "") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.5)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("figures/select/stage_one_variables.pdf", plot = pB, width = 5,
       height = 8, units = "in", dpi = 300)

X2 <- X[, c(1:45, which(colnames(X) %in% var_table$Variable))]
g2 <- ga(y, X2, wts, popsize = 800, n = 3, maxiter = 1000, run = 200, pcrossover = 0.8,
         pmutation = 0.5, elitism = 0.05, gamma = 1, verbose = TRUE)
write_rds(g2, "data/weather/ga_final_model.rds")


# Behavior of GA over replicates and final model
minima <- g %>%
  map(function(x) {
    c(x$minima, rep(NA, max(niters) - length(x$minima)))
  }) %>%
  unlist() %>%
  matrix(., nrow = max(niters), ncol = 100, byrow = FALSE)

perf <- as_tibble(minima) %>%
  gather(Replicate, Minimum, everything()) %>%
  filter(!is.na(Minimum)) %>%
  group_by(Replicate) %>%
  mutate(Iteration = 1:n()) %>%
  ungroup()
pA <- ggplot(perf, aes(x = Iteration, y = Minimum)) + theme_classic() +
  geom_line(aes(group = Replicate), colour = "grey80", alpha = 0.5) +
  geom_line(aes(y = M), data = perf %>%
              group_by(Iteration) %>%
              summarise(M = mean(Minimum))) +
  geom_line(data = tibble(Iteration = 1:length(g2$minima),
                          Minimum = g2$minima),
            colour = "red", size = 1) +
  labs(y = "eBIC")


# Evaluate all possible models --------------------------------------------
models <- lapply(1:5, function(i) {
  temp <- combn(g2$g, i)
  if (nrow(temp) < 5) {
    temp <- rbind(temp, matrix(NA, ncol = ncol(temp), nrow = 5 - i))
  }
  return(temp)
})
models <- do.call("cbind", models)

bic <- apply(models, 2, function(v) {
  v <- v[!is.na(v)]
  temp <- X[, c(1:45, which(colnames(X) %in% v))]
  BIC(lm(y ~ 0 + temp, weights = wts)) + 2*1*log(choose(ncol(X) - 45, length(v)))
})

# Add the model without environmental variables
# N.B.: The ith column of `models` corresponds to the (i+1)th element of `bic`.
bic <- c(BIC(lm(y ~ 0 + X[, 1:45], weights = wts)) + 
           2*1*log(choose(ncol(X) - 45, 0)), bic)

selected <- models[, which.min(bic) - 1]
pA <- pA + geom_hline(yintercept = min(bic), linetype = 2)
ggsave("figures/select/stage_one_GA.pdf", plot = pA, width = 5, height = 4,
       units = "in", dpi = 300)

# Plot of selected variables
pC <- tibble(V = selected[!is.na(selected)]) %>%
  separate(V, c("Variable", "Start", "End"), sep = "_", remove = TRUE) %>%
  mutate_at(c("Start", "End"), str_remove, pattern = "X") %>%
  mutate_at(c("Start", "End"), as.numeric) %>%
  mutate(Ypos = factor(Variable) %>% as.integer(),
         Ypos = Ypos + c(0.125, -0.125, 0, 0)) %>%
  ggplot(.) + theme_classic() +
    geom_segment(aes(x = Start, xend = End, y = Ypos, yend = Ypos,
                     colour = Variable), size = 2) +
    geom_vline(xintercept = 1, linetype = 2) +
    labs(x = "% CHU to anthesis", y = "") + guides(colour = "none") +
    scale_colour_manual(values = c("TMAX" = "red", "SR" = "goldenrod", "NET" = "brown")) +
    scale_y_continuous(breaks = 1:3, labels = c("NET", "SR" ,"TMAX")) +
    scale_x_continuous(limits = c(0, 1.5), labels = scales::percent)
ggsave("figures/select/selected_variables.pdf", plot = pC, width = 5, height = 4,
       units = "in", dpi = 300)

library(grid)
library(gridExtra)

lay <- matrix(c(1, 2, 3, 2), ncol = 2, byrow = TRUE)
grob1 <- grobTree(ggplotGrob(pA),
                  textGrob("A", x = unit(0.03, "npc"), y = unit(0.975, "npc"),
                           hjust = "left", vjust = "top",
                           gp = gpar(fontface = "bold", fontsize = 14)))
grob2 <- grobTree(ggplotGrob(pB),
                  textGrob("B", x = unit(0.03, "npc"), y = unit(0.975, "npc"),
                           hjust = "left", vjust = "top",
                           gp = gpar(fontface = "bold", fontsize = 14)))
grob3 <- grobTree(ggplotGrob(pC),
                  textGrob("C", x = unit(0.03, "npc"), y = unit(0.975, "npc"),
                           hjust = "left", vjust = "top",
                           gp = gpar(fontface = "bold", fontsize = 14)))
gp <- arrangeGrob(grob1, grob2, grob3, layout_matrix = lay)
ggsave("figures/select/ga_results.pdf", gp, width = 10, height = 8,
       units = "in", dpi = 300)

write_lines(selected[!is.na(selected)], "data/phenotype/selected_variables.txt")


# Replicate 30 is the best of the replicate models
selected <- selected[!is.na(selected)]
min_vars <- c(selected, g2$g, g[[30]]$g) %>% unique()
var_count <- tibble(Variable = min_vars, 
                    Replicate = Variable %in% g[[30]]$g, 
                    Second = Variable %in% g2$g, 
                    Final = Variable %in% selected)
