suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
library(parallel)

options(readr.num_columns = 0)

parser <- ArgumentParser()
parser$add_argument("-s", "--stat", type = "character")
parser$add_argument("-r", "--reference", type = "character")
parser$add_argument("-x", "--index", type = "character")
parser$add_argument("-i", "--imputed", type = "character")
parser$add_argument("-o", "--output", type = "character")
parser$add_argument("-c", "--count", type = "integer")
args <- parser$parse_args()


# Load the reference genotypes --------------------------------------------
ref <- read_tsv(args$reference, comment = "##")


# Load SNP statistics -----------------------------------------------------
s <- read_tsv(args$stat) %>%
  filter(snpid %in% ref$ID)


# Load imputed SNPs -------------------------------------------------------
impute <- read_tsv(args$imputed, comment = "##")


# Load the index ----------------------------------------------------------
idx <- read_csv(args$index, col_names = FALSE) %>%
  mutate_if(is.integer, function(x) x + 1)
idx2 <- tibble(ID = idx[[1]], 
               idx = split(t(idx[, -1]), rep(1:nrow(idx), each = ncol(idx) - 1)))
rm(idx)


# Combine the index and statistics ----------------------------------------
s <- inner_join(s, idx2, by = c("snpid" = "ID"))


# Calculate accuracy ------------------------------------------------------
cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(tidyverse))

split_idx <- clusterSplit(cl, 1:nrow(impute))
split_idx <- rep(1:detectCores(), times = sapply(split_idx, length))

ref <- split(ref, split_idx)
impute <- split(impute, split_idx)
s <- split(s, split_idx)

res <- clusterMap(cl, function(ref, impute, s) {
  map_dbl(1:nrow(s), function(i) {
    known <- ref[i, s$idx[[i]]] %>% t() %>% as.vector() %>%
      str_split(., ":") %>% sapply(., function(x) x[1]) %>%
      str_replace(., "1/0", "0/1")
    unknown <- impute[i, s$idx[[i]]] %>% t() %>% as.vector() %>%
      str_split(., ":") %>% sapply(., function(x) x[1]) %>%
      str_replace(., "\\|", "/") %>% str_replace(., "1/0", "0/1")
    sum(known == unknown)/length(known)
  })
}, ref = ref, impute = impute, s = s)
stopCluster(cl)

s <- bind_rows(s) %>%
  mutate(acc = unlist(res, use.names = FALSE)) %>%
  select(-idx)

write_csv(s, paste0(args$output, "_", args$count, ".csv"))
