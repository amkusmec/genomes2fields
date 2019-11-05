library(tidyverse)
library(sommer)


pcs <- read_rds("data/gbs/pca_covariates.rds") %>%
  as_tibble(rownames = "PedigreeNew")
rxn <- read_rds("data/phenotype/rxn_norm_parameters.rds") %>%
  spread(Variable, Value) %>%
  inner_join(., pcs, by = "PedigreeNew")

snps <- read_rds("data/gbs/add_snps.rds")
K <- A.mat(snps$GD - 1)

fix1 <- paste0("cbind(", paste(names(rxn)[2:8], collapse = ", "), ")")
fix2 <- paste(paste0("PC", 1:9), collapse = " + ")
fix <- paste0(fix1, " ~ 1 + ", fix2) %>% as.formula()

ans <- mmer(fix, random = ~ vs(PedigreeNew, Gu = K, Gtc = unsm(7)), 
            rcov = ~ vs(units, Gtc = unsm(7)), 
            data = rxn, date.warning = FALSE)
write_rds(ans, "data/phenotype/genetic_corr.rds")

bounds <- read_rds("data/phenotype/bootstrap_intervals.rds")

x <- ans$sigma$`u:PedigreeNew`
xx <- cov2cor(x)
xx[upper.tri(xx, diag = TRUE)] <- NA
xx <- xx %>%
  as_tibble(rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R))
xx %>%
  inner_join(., bounds, by = c("Phenotype1", "Phenotype2")) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R.x)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "") + 
    geom_text(aes(label = G)) +
    scale_fill_gradient2(low = "blue", high = "red") + 
    scale_x_discrete(position = "top") + 
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
ggsave("figures/select/rxn_cor_g.pdf", width = 9, height = 7, units = "in", dpi = 300)

y <- ans$sigma$`u:units`
yy <- cov2cor(y)
yy[upper.tri(yy, diag = TRUE)] <- NA
yy <- yy %>%
  as_tibble(rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R))
yy %>%
  inner_join(., bounds, by = c("Phenotype1", "Phenotype2")) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R.x)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "") +
    geom_text(aes(label = R.y)) +
    scale_fill_gradient2(low = "blue", high = "red") + 
    scale_x_discrete(position = "top") + 
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
ggsave("figures/select/rxn_cor_r.pdf", width = 9, height = 7, units = "in", dpi = 300)

z <- cov(rxn[, 2:8])
zz <- cov2cor(z)
zz[upper.tri(zz, diag = TRUE)] <- NA
zz <- zz %>%
  as_tibble(rownames = "Phenotype1") %>%
  gather(Phenotype2, R, -Phenotype1) %>%
  filter(!is.na(R))
zz %>% inner_join(., bounds, by = c("Phenotype1", "Phenotype2")) %>%
  mutate(Phenotype1 = factor(Phenotype1, levels = sort(unique(Phenotype1)), ordered = TRUE), 
         Phenotype2 = factor(Phenotype2, levels = rev(sort(unique(Phenotype2))), ordered = TRUE)) %>%
  ggplot(., aes(x = Phenotype1, y = Phenotype2, fill = R.x)) + theme_classic() +
    geom_tile() + labs(x = "", y = "", fill = "") +
    geom_text(aes(label = R.y)) +
    scale_fill_gradient2(low = "blue", high = "red") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.0625))
ggsave("figures/select/rxn_cor_p.pdf", width = 9, height = 7, units = "in", dpi = 300)


# Correlations as graphs --------------------------------------------------
# Threshold by bootstrapped confidence intervals
P <- cov2cor(z); G <- cov2cor(x); R <- cov2cor(y)
for (i in rownames(G)) {
  for (j in colnames(G)) {
    if (i == j) {
      P[i, j] <- G[i, j] <- R[i, j] <- 0
    } else {
      r <- which((bounds$Phenotype1 == i & bounds$Phenotype2 == j) |
                   (bounds$Phenotype2 == i & bounds$Phenotype1 == j))
      P[i, j] <- if_else(bounds$P_lower[r] > 0 | bounds$P_upper[r] < 0, 
                         P[i, j], 0)
      G[i, j] <- if_else(bounds$G_lower[r] > 0 | bounds$G_upper[r] < 0, 
                         G[i, j], 0)
      R[i, j] <- if_else(bounds$R_lower[r] > 0 | bounds$R_upper[r] < 0, 
                         R[i, j], 0)
    }
  }
}

library(igraph)
vertex_colors <- c("lightgreen", "lightsteelblue2", "orange", "wheat1")

rownames(G) <- colnames(G) <- 
  c("Hybrid", "lnMSE", "NET1", "NET2", "TMAX", "TMIN1", "TMIN2")
Gnet <- graph_from_adjacency_matrix(G, mode = "undirected", weighted = TRUE)
E(Gnet)$sign <- sign(E(Gnet)$weight)
E(Gnet)$color <- if_else(E(Gnet)$sign < 0, "red", "blue")
E(Gnet)$weight <- abs(E(Gnet)$weight)
E(Gnet)$width <- E(Gnet)$weight*10

Gc <- cluster_fast_greedy(Gnet)
V(Gnet)$color <- vertex_colors[membership(Gc)]
E(Gnet)$lty <- crossing(Gc, Gnet) + 1

pdf("figures/select/network_g.pdf", width = 7, height = 7)
plot(Gnet, layout = layout_with_fr, edge.curved = 0.1, vertex.size = 25, 
     main = "Genetic Correlations")
dev.off()


rownames(R) <- colnames(R) <- 
  c("Hybrid", "lnMSE", "NET1", "NET2", "TMAX", "TMIN1", "TMIN2")
Rnet <- graph_from_adjacency_matrix(R, mode = "undirected", weighted = TRUE)
E(Rnet)$sign <- sign(E(Rnet)$weight)
E(Rnet)$color <- if_else(E(Rnet)$sign < 0, "red", "blue")
E(Rnet)$weight <- abs(E(Rnet)$weight)
E(Rnet)$width <- E(Rnet)$weight*10

Rc <- cluster_fast_greedy(Rnet)
V(Rnet)$color <- vertex_colors[membership(Rc)]
E(Rnet)$lty <- crossing(Rc, Rnet) + 1

pdf("figures/select/network_r.pdf", width = 7, height = 7)
plot(Rnet, layout = layout_with_fr, edge.curved = 0.1, vertex.size = 25, 
     main = "Residual Correlations")
dev.off()


rownames(P) <- colnames(P) <- 
  c("Hybrid", "lnMSE", "NET1", "NET2", "TMAX", "TMIN1", "TMIN2")
Pnet <- graph_from_adjacency_matrix(P, mode = "undirected", weighted = TRUE)
E(Pnet)$sign <- sign(E(Pnet)$weight)
E(Pnet)$color <- if_else(E(Pnet)$sign < 0, "red", "blue")
E(Pnet)$weight <- abs(E(Pnet)$weight)
E(Pnet)$width <- E(Pnet)$weight*10

Pc <- cluster_fast_greedy(Pnet)
V(Pnet)$color <- vertex_colors[membership(Pc)]
E(Pnet)$lty <- crossing(Pc, Pnet) + 1

pdf("figures/select/network_p.pdf", width = 7, height = 7)
plot(Pnet, layout = layout_with_fr, edge.curved = 0.1, vertex.size = 25, 
     main = "Phenotypic Correlations")
dev.off()


# Structure of genetic correlation matrix ---------------------------------
G_svd <- svd(x)

# G is ill-conditioned and most of the genetic variances (99.9%) lies along
# the first eigenvector
max(G_svd$d)/min(G_svd$d)
G_svd$d^2/sum(G_svd$d^2)

tibble(Phenotype = rownames(x), 
       E = G_svd$v[, 1]/max(abs(G_svd$v[, 1]))) %>%
  ggplot(., aes(x = Phenotype, y = E)) + theme_classic() + 
    geom_col(colour = "black") + labs(x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# "autonomy" (Hansen and Houle 2008)
auto <- map_dbl(1:7, function(i) {
  1 - (1/x[i, i])*(x[i, -i, drop = FALSE] %*% solve(x[-i, -i]) %*% t(x[i, -i, drop = FALSE]))
})
names(auto) <- rownames(x)

# modular autonomy
m1 <- c(4, 6, 7) # NET2, TMIN1, TMIN2
m2 <- c(1, 3, 5) # Hybrid, NET1, TMAX
a1 <- x[m1, m1] - x[m1, m2] %*% solve(x[m2, m2]) %*% t(x[m1, m2])
a2 <- x[m2, m2] - x[m2, m1] %*% solve(x[m1, m1]) %*% t(x[m2, m1])
a1/x[m1, m1]
a2/x[m2, m2]

# (integrated) 0 <= a <= 1 (autonomous)
# proportion of genetic variance independent from other traits
auto
tibble(Phenotype = rownames(x), 
       Autonomy = auto, 
       Fill = as.character(membership(Gc))) %>% 
  inner_join(., tibble(Phenotype = c(colnames(a1), colnames(a2), "lnMSE"), 
                       CA = c(diag(a1/x[m1, m1]), diag(a2/x[m2, m2]), NA)),
             by = "Phenotype") %>%
  ggplot(., aes(x = Phenotype, y = Autonomy)) + theme_classic() +
    geom_col(aes(fill = Fill), colour = "black") + 
    geom_point(aes(y = CA)) + guides(fill = FALSE) +
    scale_fill_manual(values = c("1" = "lightgreen", "2" = "lightsteelblue2", 
                                 "3" = "orange")) +
    labs(x = "", y = expression(paste("% ", sigma[g]^2, "(y|X)"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/select/genetic_integration.pdf", width = 6, height = 4, units = "in", dpi = 300)
