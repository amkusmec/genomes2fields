library(tidyverse)
library(purrrlyr)
library(igraph)


# Create the graph --------------------------------------------------------
ppi_txt <- read_tsv("data/highConfidentPPIs.txt", col_names = c("Node1", "Node2"))
proteins <- c(unique(ppi_txt$Node1), unique(ppi_txt$Node2)) %>%
  unique() %>% sort()
ppi_adj <- matrix(0, nrow = length(proteins), ncol = length(proteins), 
                  dimnames = list(proteins, proteins))
for (i in 1:nrow(ppi_txt)) {
  if (i %% 500 == 0) cat(i, "\r")
  ppi_adj[ppi_txt$Node1[i], ppi_txt$Node2[i]] <- 1
  ppi_adj[ppi_txt$Node2[i], ppi_txt$Node1[i]] <- 1
}
ppi <- graph_from_adjacency_matrix(ppi_adj, mode = "undirected", diag = FALSE)

# pc <- cluster_edge_betweenness(ppi, weights = NULL, directed = FALSE)
# write_rds(pc, "data/ppi_clustering.rds")
pc <- read_rds("data/ppi_clustering.rds")
pc <- cluster_leading_eigen(ppi)


cand <- read_rds("data/gemma/nearest_genes.rds")
gene_table <- read_rds("data/gemma/candidate_genes.rds") %>%
  by_row(function(r) {
    sapply(unique(cand$Gene), function(g) str_detect(r$Genes, g)) %>%
      sum()
  }, .to = "NGenes") %>%
  unnest(NGenes) %>%
  mutate(Genes = str_split(Genes, ", "))
bgenes <- unlist(gene_table$Genes, use.names = FALSE) %>% unique()

modules <- enframe(membership(pc)) %>%
  rename(Protein = name, Module = value) %>%
  mutate(Module = as.integer(Module), 
         Gene = if_else(str_detect(Protein, "AC"), str_remove(Protein, "P"), 
                        str_remove(Protein, "_P[0-9]{2}")))

cand_modules <- filter(modules, Gene %in% cand$Gene) %>% 
  pull(Module) %>% unique()
ld_modules <- filter(modules, Gene %in% bgenes) %>%
  pull(Module) %>% unique()

source("src/go_enrich.R")
module_enrich <- cand_modules %>%
  map_dbl(function(m) {
    k <- sum(unique(cand$Gene) %in% modules$Gene[modules$Module == m])
    K <- length(unique(cand$Gene))
    n <- sum(modules$Module == m)
    N <- length(unique(modules$Protein))
    1 - phyper(k, K, N - K, n)
  }) %>%
  p.adjust(method = "fdr")
enriched_modules <- cand_modules[module_enrich <= 0.05]

# Both enriched modules have 1 gene in them
# cand_enrich <- enriched_modules %>%
#   map(function(m) {
#     filter(modules, Module == m) %>% pull(Gene) %>% go_enrich()
#   })
# names(cand_enrich) <- make.names(enriched_modules)

ld_module_enrich <- ld_modules %>%
  map_dbl(function(m) {
    k <- sum(bgenes %in% modules$Gene[modules$Module == m])
    K <- length(bgenes)
    n <- sum(modules$Module == m)
    N <- length(unique(modules$Protein))
    1 - phyper(k, K, N - K, n)
  }) %>%
  p.adjust(method = "fdr")
ld_enriched_modules <- ld_modules[ld_module_enrich <= 0.05]

# All modules have <= 3 genes
# ld_enrich <- setdiff(ld_enriched_modules, enriched_modules) %>%
#   map(function(m) {
#     if (m == 1244) {
#       NULL
#     } else {
#       filter(modules, Module == m) %>% pull(Gene) %>% go_enrich()
#     }
#   })
# names(ld_enrich) <- make.names(ld_enriched_modules)
# ld_enrich <- ld_enrich[!sapply(ld_enrich, is.null)]

terms <- cand_enrich %>%
  map(function(df) filter(df, q_val <= 0.05) %>% pull(Name))

# Terms unique to modules
u <- map(seq_along(terms), function(i) setdiff(terms[[i]], reduce(terms[-i], union)))
names(u) <- names(terms)

# Terms shared by at least 2 modules
u2 <- setdiff(reduce(terms, union), unlist(u, use.names = FALSE))

write_rds(list(Unique = u, Shared = u2), "data/gemma/ppi_enrich.rds")
