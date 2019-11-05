source("src/sim_setup.R")
source("src/plot_sim.R")
ngen <- 15

# Simulation-specific random seed
set.seed(906565)

res <- list()
for (i in 1:ncol(populations)) {
  cat("Replicate:", i, "\r")
  
  pop <- newPop(founderPop)
  pop <- pop[populations[, i]]
  genMean <- list(meanG(pop))
  genVar <- list(varA(pop))
  
  for (j in 1:ngen) {
    pop2 <- selectInd(pop, nInd = 40, trait = 1, use = "bv")
    pop <- randCross(pop2, nCrosses = 400, nProgeny = 1)
    genMean[[j + 1]] <- meanG(pop)
    genVar[[j + 1]] <- varA(pop)
  }
  
  res[[i]] <- list(genMean = genMean, genVar = genVar)
}
write_rds(res, "data/simulation/sim_yield_nE.rds")

p <- plot_sim(res, scaled = TRUE, pheno = names(rxn)[2:8])
ggsave("figures/simulation/sim_yield_nE.pdf", p, 
       width = 10, height = 6, units = "in", dpi = 300)
