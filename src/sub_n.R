### Variable definitions:
###  - n = integer size of desired training set
###  - groups = named character vector of group assignments; names match taxa

# Calculate the number of individuals to select from each group
# Used for stratified methods (except uniform)
sub_n <- function(groups, n) {
  tt <- table(groups)
  ng <- n*(log(tt)/sum(log(tt)))
  ng <- round(ng, 0)
  ng[ng > tt] <- tt[ng > tt]
  ng
}
