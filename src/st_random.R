### Variable definitions:
###  - taxa = character vector of sample IDs; used as row/col names for matrices
###  - n = integer size of desired training set
###  - groups = named character vector of group assignments; names match taxa

source("~/genomes2fields/src/sub_n.R")
source("~/genomes2fields/src/random.R")

# Select a random subset of size n stratified by a grouping variable
st_random <- function(taxa, groups, n) {
  ng <- sub_n(groups, n)
  unlist(sapply(sort(unique(groups)), 
                function(g) random(taxa[groups == g], ng[g])),
         use.names = FALSE)
}
