library(tidyverse)
library(SGL)
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

# Create an index of variable groups
group_idx <- colnames(X) %>%
  str_split(., "_") %>%
  sapply(., function(x) x[1]) %>%
  factor() %>% as.integer()

# For linear regression, SGL automatically centers the response, so no intercept
# needs to be added to the design matrix.

s <- site_levels[1]
idx <- which(sites == s)

test_x <- X[idx, ]
test_y <- y[idx]

train_x <- X[-idx, ]
train_y <- y[-idx]

var <- apply(train_x, 2, var)
sum(near(var, 0))

### SGL is introducing NA's into X during centering/standardization
sgl <- SGL(list(x = train_x, y = train_y), group_idx)
