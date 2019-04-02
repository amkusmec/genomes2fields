library(tidyverse)
library(randomForest)
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


# Random forest with leave-one-site-out CV --------------------------------
cl <- makeCluster(length(site_levels))
clusterEvalQ(cl, library(randomForest))
clusterExport(cl, list("sites", "y", "X"))

res <- parLapply(cl, site_levels, function(s) {
  # Partition into training and testing sets
  idx <- which(sites == s)
  
  test_x <- X[idx, ]
  test_y <- y[idx]
  
  train_x <- X[-idx, ]
  train_y <- y[-idx]
  
  # Train the random forest
  rf <- randomForest(train_x, train_y, ntree = 1000, replace = TRUE, 
                     importance = TRUE, nPerm = 5, do.trace = FALSE, 
                     keep.forest = TRUE)
  
  # Compute MSE on the left out site
  pred <- predict(rf, test_x)
  mse <- mean((pred - test_y)^2)
  
  return(list(rf = rf, mse = mse))
})

stopCluster(cl)
write_rds(res, "data/weather/rf_select.rds")


# Random forest with all data ---------------------------------------------
rf <- randomForest(X, y, ntree = 1000, replace = TRUE, importance = TRUE, 
                   nPerm = 5, do.trace = FALSE, keep.forest = TRUE)
write_rds(rf, "data/weather/rf_all.rds")


# Concordance of variable ranks -------------------------------------------
imp <- sapply(res, function(x) importance(x[[1]])[, 1])

wt <- sapply(res, function(x) x[[2]])
wt <- (1/(wt/sum(wt)))/sum(1/(wt/sum(wt)))
comp <- rowSums(imp*matrix(wt, nrow = nrow(imp), ncol = ncol(imp), byrow = TRUE))

imp <- cbind(importance(rf)[, 1], comp, imp)
dimnames(imp) <- list(colnames(X), c("All", "Composite", site_levels))

imp_cor_tau <- cor(imp, method = "kendall")

as_tibble(imp_cor_tau, rownames = "Site1") %>%
  gather(Site2, Tau, -Site1) %>%
  ggplot(., aes(x = Site1, y = Site2, fill = Tau)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = expression(tau[b])) +
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) +
    theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/select/rf_importance_cor.pdf", width = 11, height = 8, 
       units = "in", dpi = 300)

imp_cor_r <- cor(imp, method = "pearson")

as_tibble(imp_cor_r, rownames = "Site1") %>%
  gather(Site2, R, -Site1) %>%
  ggplot(., aes(x = Site1, y = Site2, fill = R)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "r") +
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-1, 1)) +
    theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("figures/select/rf_importance_cor2.pdf", width = 11, height = 8, 
       units = "in", dpi = 300)
