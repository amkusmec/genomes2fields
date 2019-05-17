# Function to construct the direct sum of a list of matrices
direct_sum <- function(matrices) {
  n <- sum(sapply(matrices, nrow))
  m <- sum(sapply(matrices, ncol))
  
  temp <- matrix(0, nrow = n, ncol = m)
  for (i in seq_along(matrices)) {
    if (i == 1) {
      start_r <- 0
      start_c <- 0
    }
    else {
      start_r <- sum(sapply(matrices[1:(i - 1)], nrow))
      start_c <- sum(sapply(matrices[1:(i - 1)], ncol))
    }
    
    temp[(start_r + 1):(start_r + nrow(matrices[[i]])),
         (start_c + 1):(start_c + ncol(matrices[[i]]))] <- matrices[[i]]
  }
  
  return(temp)
}
