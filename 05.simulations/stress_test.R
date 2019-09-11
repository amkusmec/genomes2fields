library(AlphaSimR)


haplotypes <- readr::read_rds("data/gbs/synthetic_hybrids_haplotypes.rds")
genMap <- readr::read_rds("data/gbs/sim_map.rds")
for (i in seq_along(genMap)) {
  genMap[[i]] <- genMap[[i]] - min(genMap[[i]])
}

G <- 1.5*diag(2) - 0.5
E <- diag(2)

founderPop <- newMapPop(genMap, haplotypes, inbred = FALSE)
SP <- SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = 1000, mean = c(0, 0), var = c(1, 2), corA = G)
SP$setVarE(h2 = c(0.5, 0.5))

pop <- newPop(founderPop)
pop <- pop[sample(701, 400, replace = FALSE)]
genMean <- list(meanG(pop))
genVar <- list(varA(pop))

start_time <- proc.time()[3]
for (i in 1:50) {
  b <- smithHazel(c(1, 2), varG(pop), varP(pop))
  pop2 <- selectInd(pop, nInd = 40, trait = selIndex, b = b)
  pop <- randCross(pop2, nCrosses = 400)
  genMean[[i + 1]] <- meanG(pop)
  genVar[[i + 1]] <- varA(pop)
}
end_time <- proc.time()[3]

end_time - start_time

plot(0:50, sapply(genMean, function(x) x[1]), type = "l", ylab = "Gen. Mean", 
     ylim = c(-10, 50))
lines(0:50, sapply(genMean, function(x) x[2]), col = "red")

plot(0:50, sapply(genVar, function(x) x[1, 1]), type = "l", ylab = "Gen. Var", 
     ylim = c(0, 3))
lines(0:50, sapply(genVar, function(x) x[2, 2]), col = "red")
