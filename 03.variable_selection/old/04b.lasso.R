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


# Analysis ----------------------------------------------------------------
# Number of non-zero coefficients selected by cross-validation
nzero <- tibble(Site = site_levels, 
                NZero = map_int(res, function(r) {
                  idx <- which(r$lasso$lambda == r$lasso$lambda.1se)
                  r$lasso$nzero[idx]
                }))
ggplot(nzero, aes(x = Site, y = NZero)) + theme_classic() +
  geom_point(size = 3) + theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
  geom_hline(yintercept = lasso$nzero[which(lasso$lambda == lasso$lambda.1se)], 
             linetype = 2, colour = "red") +
  labs(x = "", y = "# Non-Zero Coefficients")
ggsave("figures/select/lasso_nzero.pdf", width = 10, height = 6, units = "in", dpi = 300)

mean(nzero$NZero); var(nzero$NZero)

# Variables selected by leave-one-out models
variables <- lapply(res, function(r) {
  idx <- max(which(r$lasso$nzero <= 5))
  r$lasso$glmnet.fit$beta@Dimnames[[1]][which(r$lasso$glmnet.fit$beta[, idx] != 0)]
})
names(variables) <- site_levels

length(unique(unlist(variables, use.names = FALSE)))
top5 <- names(sort(table(unlist(variables, use.names = FALSE)), decreasing = TRUE)[1:5])
r2adj_top5 <- lm(y ~ X[, top5]) %>% broom::glance() %>% pull(adj.r.squared)

# Quality of the models
r2adj <- sapply(variables, function(v) {
  lm(y ~ X[, v]) %>% broom::glance() %>% pull(adj.r.squared)
})
r2adj_all <- lm(y ~ X[, lasso$glmnet.fit$beta@Dimnames[[1]][which(lasso$glmnet.fit$beta[, max(which(lasso$nzero <= 5))] != 0)]]) %>%
  broom::glance() %>% pull(adj.r.squared)

tibble(R2 = r2adj) %>%
  ggplot(., aes(x = R2)) + theme_classic() +
    geom_histogram(binwidth = 0.01, fill = "skyblue", colour = "black") +
    geom_vline(xintercept = round(r2adj_all, 2), linetype = 2, colour = "orange") +
    labs(x = expression(R[adj]^2), y = "Count")
ggsave("figures/select/lasso_r2.pdf", width = 6, height = 4, units = "in", dpi = 300)

sum(r2adj >= r2adj_all)/length(site_levels)

best_r2 <- which.max(r2adj)
site_levels[best_r2]
res[[best_r2]]$lasso$glmnet.fit$beta@Dimnames[[1]][which(res[[best_r2]]$lasso$glmnet.fit$beta[, max(which(res[[best_r2]]$lasso$nzero <= 5))] != 0)]

best_mse <- which.min(sapply(res, function(x) x$mse))
site_levels[best_mse]
res[[best_mse]]$lasso$glmnet.fit$beta@Dimnames[[1]][which(res[[best_mse]]$lasso$glmnet.fit$beta[, max(which(res[[best_mse]]$lasso$nzero <= 5))] != 0)]
