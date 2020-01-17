library(tidyverse)
library(WGCNA)

allowWGCNAThreads(nThreads = 8)

# xref <- read_tsv("~/anno/gene_model_xref_v2.txt", skip = 4) %>%
#   select(v2_gene_model:v2_chr, v3_gene_model) %>%
#   mutate(v2_chr = str_remove(v2_chr, "chr") %>% as.integer())
# files <- list.files("~/expr_plasticity/data/expression/TASSEL_formatted", 
#                     "*\\.txt", full.names = TRUE)
# expr <- read_tsv(files[1], col_names = TRUE, skip = 2)
# expr0 <- expr[, names(expr) %in% xref$v3_gene_model]
# rownames(expr0) <- expr$Taxa
# idx <- match(names(expr0), xref$v3_gene_model)
# names(expr0) <- xref$v2_gene_model[idx]

expr0 <- read_rds("~/SAM/networks/sam_rpkm.rds")

gsg <- goodSamplesGenes(expr0, verbose = 3)
if (!gsg$allOK) {
  expr0 <- expr0[gsg$goodSamples, gsg$goodGenes]
}

# powers <- c(c(1:10), seq(12, 20, 2))
# sft <- pickSoftThreshold(expr0, powerVector = powers, verbose = 5)
# adj <- adjacency(expr0, power = 5)
# dissTOM <- 1 - TOMsimilarity(adj)
# write_rds(dissTOM, "data/SAM_dissTOM.rds")
dissTOM <- read_rds("data/SAM_dissTOM.rds")

geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05)

MEList <- moduleEigengenes(expr0, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss = 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
# mergeList <- seq(0.1, 0.9, 0.05) %>%
#   map(function(i) {
#     mergeCloseModules(expr0, dynamicColors, cutHeight = i, verbose = 3)
#   })
# merge_num <- sapply(mergeList, function(m) ncol(m$newMEs))
# plot(seq(0.1, 0.9, 0.05), merge_num)
merge <- mergeCloseModules(expr0, dynamicColors, cutHeight = 0.6, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05)
moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

source("src/go_enrich_at.R")
cand <- read_rds("data/gemma/nearest_genes.rds")
gene_table <- read_rds("data/gemma/candidate_genes.rds") %>%
  by_row(function(r) {
    sapply(unique(cand$Gene), function(g) str_detect(r$Genes, g)) %>%
      sum()
  }, .to = "NGenes") %>%
  unnest(NGenes) %>%
  mutate(Genes = str_split(Genes, ", "))
bgenes <- unlist(gene_table$Genes, use.names = FALSE) %>% unique()

cand_modules <- mergedColors[names(expr0) %in% unique(cand$Gene)] %>%
  unique()
ld_modules <- mergedColors[names(expr0) %in% bgenes] %>% unique()

module_enrich <- cand_modules %>%
  map_dbl(function(m) {
    k <- sum(unique(cand$Gene) %in% names(expr0)[mergedColors == m])
    K <- length(unique(cand$Gene))
    n <- sum(mergedColors == m)
    N <- ncol(expr0)
    1 - phyper(k, K, N - K,  n)
  }) %>%
  p.adjust(method = "fdr")
enriched_modules <- cand_modules[module_enrich <= 0.05]

cand_enrich <- enriched_modules %>%
  map(function(m) {
    names(expr0)[mergedColors == m] %>% go_enrich()
  })
names(cand_enrich) <- enriched_modules

ld_module_enrich <- ld_modules %>%
  map_dbl(function(m) {
    k <- sum(bgenes %in% names(expr0)[mergedColors == m])
    K <- length(bgenes)
    n <- sum(mergedColors == m)
    N <- ncol(expr0)
    1 - phyper(k, K, N - K, n)
  }) %>%
  p.adjust(method = "fdr")
ld_enriched_modules <- ld_modules[ld_module_enrich <= 0.05]

ld_enrich <- setdiff(ld_enriched_modules, enriched_modules) %>%
  map(function(m) {
    names(expr0)[mergedColors == m] %>% go_enrich()
  })
names(ld_enrich) <- setdiff(ld_enriched_modules, enriched_modules)

terms <- c(cand_enrich, ld_enrich) %>%
# terms <- cand_enrich %>%
  map(function(df) filter(df, q_val <= 0.05) %>% pull(Name))

# Terms unique to modules
u <- map(seq_along(terms), function(i) setdiff(terms[[i]], reduce(terms[-i], union)))
names(u) <- names(cand_enrich)

# Terms shared by at least 2 modules
# u2 <- setdiff(reduce(terms, union), unlist(u, use.names = FALSE))

write_rds(list(cand_enrich = cand_enrich, ld_enrich = ld_enrich, 
               Unique = u), "data/gemma/coexpress_enrich.rds")
