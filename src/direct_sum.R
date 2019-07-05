require(tidyverse)

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
