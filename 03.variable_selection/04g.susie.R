library(tidyverse)
library(susieR)


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


# Estimate the scaled prior variance --------------------------------------
pve <- apply(X, 2, function(x) lm(y ~ x) %>% broom::glance() %>% pull(r.squared))
prior_effect_var <- mean(pve)



# Leave-one-site-out cross-validation -------------------------------------
res <- lapply(site_levels, function(s) {
  cat("\n", which(site_levels == s), "\n")
  
  # Partition into training and testing sets
  idx <- which(sites == s)
  
  test_x <- X[idx, ]
  test_y <- y[idx]
  
  train_x <- X[-idx, ]
  train_y <- y[-idx]
  
  # SuSiE
  susie(train_x, train_y, L = 5, scaled_prior_variance = prior_effect_var, 
        estimate_residual_variance = TRUE, estimate_prior_variance = TRUE)
})

write_rds(res, "data/weather/susie_select.rds")


# SuSiE with all data -----------------------------------------------------
sus <- susie(X, y, L = 5, scaled_prior_variance = prior_effect_var,
             estimate_residual_variance = TRUE, estimate_prior_variance = TRUE)
write_rds(sus, "data/weather/susie_all.rds")

lm(y ~ X[, sort(unlist(sus$sets$cs, use.names = FALSE))[-1]]) %>%
  broom::glance() %>% pull(adj.r.squared)
