library(tidyverse)
library(parallel)
library(genlasso)


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

cnames <- tibble(Name = colnames(X)) %>%
  separate(Name, c("Variable", "Start", "End"), sep = "_", remove = TRUE) %>%
  mutate(Start = str_replace(Start, "X", "") %>% as.double(), 
         End = str_replace(End, "X", "") %>% as.double(), 
         Width = End - Start)
X <- X[, near(cnames$Width, 0.025)]
cnames <- filter(cnames, near(Width, 0.025))

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

# Compute V^(-0.5)
e <- eigen(R)
V <- e$vectors %*% diag(1/sqrt(e$values)) %*% t(e$vectors)

z <- drop(V %*% y)
W <- V %*% X

# Center (and scale)
y_sc <- z - mean(z)
X_sc <- scale(W, center = TRUE, scale = TRUE)


# Set up the penalty matrix -----------------------------------------------
# Transitions, i.e., identify where the variables change
trans <- sapply(unique(cnames$Variable), function(v) max(which(cnames$Variable == v)))

# Basic fused LASSO penalty matrix
D <- -1*diag(1, ncol(X), ncol(X))
for (i in 1:(nrow(D) - 1)) { 
  if (!(i %in% trans)) {
    D[i, i + 1] <- 1
  }
}


# Fused LASSO -------------------------------------------------------------
# Grid of gamma values to explore for sparsity
gam <- c(1:5, 10, 15, 20)

# 1) Compute the regularization paths
# cl <- makeCluster(length(site_levels))
# clusterEvalQ(cl, library(genlasso))
# clusterExport(cl, list("sites", "X_sc", "y_sc", "D"))

paths <- list()
for (s in site_levels) {
  # cat("-", s, "\n")
  idx <- which(sites == s)
  
  # Get the training set
  train_x <- X_sc[-idx, ]
  train_y <- y_sc[-idx]
  
  paths[[s]] <- fusedlasso(train_y, train_x, D, gamma = 0, verbose = TRUE)
  # # Fit regularization paths under different levels of sparsity
  # sparse <- lapply(gam, function(g) {
  #   # cat("--", g, "\n")
  #   fusedlasso(train_y, train_x, D, gamma = g, verbose = FALSE)
  # })
  # names(sparse) <- paste0("G", seq_along(sparse))
  # 
  # return(sparse)
}
# stopCluster(cl)

names(paths) <- site_levels
# write_rds(paths, "data/weather/fused_lasso_paths.rds")
write_rds(paths, "data/weather/fused_lasso_paths0.rds")

plot_fused <- function(obj, df = 5) { 
  plot(coef(obj, df = df)$beta[, 1], type = "b", col = as.integer(factor(cnames$Variable)))
  abline(h = 0, lty = 2)
  abline(v = 40 + 60*(0:4), lty = 4) 
}

# 2) Extract the ranges of lambda and construct a grid on it
lambda_bounds <- lapply(paths, function(p) {
  # l <- lapply(p, function(g) range(g$lambda))
  # do.call("rbind", l)
  range(p$lambda)
})
lambda_bounds <- do.call("rbind", lambda_bounds)
# lambda_bounds <- c(quantile(lambda_bounds[, 1], probs = 0.25),
#                    quantile(lambda_bounds[, 2], probs = 0.75))
lambda_bounds <- c(max(lambda_bounds[, 1]), min(lambda_bounds[, 2]))

# lambda_grid <- seq(lambda_bounds[1], lambda_bounds[2], length.out = 40)
lambda_grid <- seq(lambda_bounds[1], lambda_bounds[2], length.out = 50)

# 3) Calculate the prediction error for each fold and combination of parameters
pe <- lapply(site_levels, function(s) {
  # Get the test set data
  idx <- which(sites == s)
  test_x <- X_sc[idx, ]
  test_y <- y_sc[idx]
  
  # Cross-validation error for each combination of gamma and lambda
  # cv <- lapply(paths[[s]], function(g) {
  #   # Get the coefficients for the given lambdas
  #   coefs <- coef(g, lambda = lambda_grid)
  # 
  #   # Calculate the errors
  #   yhat <- test_x %*% coefs$beta
  # 
  #   # Total error
  #   cvt <- apply(yhat, 2, function(x) sum((test_y - x)^2))
  #   names(cvt) <- paste0("L", seq_along(cvt))
  # 
  #   # Mean error
  #   cvk <- cvt/length(test_y)
  # 
  #   return(list(cvt = cvt, cvk = cvk))
  # })
  
  coefs <- coef(paths[[s]], lambda = lambda_grid)
  yhat <- test_x %*% coefs$beta
  cvt <- apply(yhat, 2, function(x) sum((test_y - x)^2))
  names(cvt) <- paste0("L", seq_along(cvt))
  cvk <- cvt/length(test_y)
  
  # return(cv)
  return(list(cvt = cvt, cvk = cvk))
})
names(pe) <- site_levels

# 4) Compute the CV means and standard errors
# pe <- lapply(pe, function(x) {
#   list(cvt = do.call("rbind", lapply(x, function(y) y$cvt)), 
#        cvk = do.call("rbind", lapply(x, function(y) y$cvk)))
# })
# 
# cv <- matrix(0, nrow = length(gam), ncol = length(lambda_grid))
# cv_se <- matrix(0, nrow = length(gam), ncol = length(lambda_grid))
# for (i in 1:nrow(cv)) {
#   for (j in 1:ncol(cv)) {
#     cv[i, j] <- sapply(pe, function(x) x$cvt[i, j]) %>% sum()
#     cv_se[i, j] <- sapply(pe, function(x) x$cvk[i, j]) %>% sd()
#   }
# }

pe <- list(cvt = do.call("rbind", lapply(pe, function(y) y$cvt)), 
           cvk = do.call("rbind", lapply(pe, function(y) y$cvk)))

cv <- apply(pe$cvt, 2, sum)/length(y)
cv_se <- apply(pe$cvk, 2, sd)/sqrt(length(site_levels))

# 5) Select the 1 SE CV parameters
cv_1 <- cv + cv_se
# rownames(cv) <- rownames(cv_1) <- paste0("G", seq_along(gam))
# colnames(cv) <- colnames(cv_1) <- paste0("L", seq_along(lambda_grid))

tibble(MSE = cv, SE = cv_se, Lambda = lambda_grid) %>% 
  ggplot(., aes(x = Lambda, y = MSE)) + theme_classic() + 
    geom_pointrange(aes(ymin = MSE - SE, ymax = MSE + SE)) + 
    geom_line() + labs(x = expression(lambda), y = "CV-MSE")
ggsave("figures/select/fused_lasso_cv0.pdf", width = 6, height = 4, 
       units = "in", dpi = 300)

cv_min <- which.min(cv)
cv_1se <- which(cv <= cv_1[cv_min])
idx <- which.max(cv[cv_1se])
# idx2 <- arrayInd(cv_min, .dim = dim(cv))
# 
# rbind(cv, cv_1) %>%
#   as_tibble(rownames = "Gamma") %>%
#   mutate(Table = rep(c("CV", "CV+1SE"), each = 8)) %>%
#   gather(Lambda, Value, L1:L40) %>%
#   mutate(Gamma = str_replace(Gamma, "G", "") %>% as.integer(), 
#          Gamma = gam[Gamma],
#          Gamma = factor(Gamma, levels = gam, ordered = TRUE), 
#          Lambda = str_replace(Lambda, "L", "") %>% as.integer(), 
#          Lambda = lambda_grid[Lambda], 
#          Region = if_else(Value <= cv_1[cv_min] & Table == "CV", 
#                           "Yes", as.character(NA)), 
#          Region = if_else(Value == cv[cv_min] & Table == "CV", 
#                           "Min", Region), 
#          Region = if_else(Value == cv[idx] & Table == "CV", 
#                           "CV", Region)) %>%
#   ggplot(., aes(x = Gamma, y = Lambda, fill = Value)) + theme_bw() +
#     geom_tile() + geom_point(aes(shape = Region)) + facet_wrap(~ Table) + 
#     labs(fill = "CV-MSE", x = expression(gamma), y = expression(lambda)) +
#     scale_fill_distiller(type = "seq", palette = "Reds", direction = 1) +
#     scale_shape_manual(values = c("Yes" = 3, "Min" = 16, "CV" = 8)) +
#     guides(shape = FALSE)
# ggsave("figures/select/fused_lasso_cv.pdf", width = 10, height = 8, 
#        units = "in", dpi = 300)


# Fit of the optimal regularization ratio to all the data -----------------
# fused_all <- fusedlasso(y_sc, X_sc, D, gamma = gam[idx2[1, 1]], verbose = TRUE)
# coefs <- coef(fused_all, lambda = lambda_grid[idx2[1, 2]])

fused_all <- fusedlasso(y_sc, X_sc, D, gamma = 0, verbose = TRUE)
coefs <- coef(fused_all, lambda = lambda_grid[idx])

plot(coefs$beta[, 1], type = "b", col = as.integer(factor(cnames$Variable)))
abline(h = 0, lty = 2); abline(v = 40 + 60*(0:4), lty = 4)
