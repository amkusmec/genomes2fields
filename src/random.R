### Variable definitions:
###  - taxa = character vector of sample IDs; used as row/col names for matrices
###  - n = integer size of desired training set

# Select a random subset of size n from taxa
random <- function(taxa, n) sample(taxa, n, replace = FALSE)
