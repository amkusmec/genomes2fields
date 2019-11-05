# Add a public function to `SimParam-class` to allow the addition of arbitrary
# traits with additive and dominance components

require(AlphaSimR)

SimParam$set(
  "public", 
  "addTraitAD_manual", 
  function(qtlPos, addEff, domEff, mean = 0, var = 1, useVarA = TRUE, force = FALSE) {
    # Check for a running simulation
    if (!force) {
      private$.isRunning()
    }
    
    # Basic error checking
    stopifnot(length(qtlPos) == private$.nChr)
    stopifnot(ncol(addEff) == ncol(domEff))
    
    # Construct a custom `LociMap` object
    qtlLoci <- new("LociMap", 
                   nLoci = as.integer(sum(sapply(qtlPos, length))), 
                   lociPerChr = as.integer(sapply(qtlPos, length)), 
                   lociLoc = as.integer(do.call("c", qtlPos)))
    
    nTraits <- length(mean)
    for (i in 1:nTraits) {
      trait = new("TraitAD", 
                  qtlLoci, 
                  addEff = addEff[, i], 
                  domEff = domEff[, i], 
                  intercept = 0)
      tmp <- calcGenParam(trait, private$.founderPop, self$nThreads)
      
      if (useVarA) {
        scale <- sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
      } else {
        scale <- sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
      }
      
      trait@addEff <- trait@addEff*scale
      trait@domEff <- trait@domEff*scale
      trait@intercept <- mean[i] - mean(tmp$gv*scale)
      
      if (useVarA) {
        private$.addTrait(trait, var[i], popVar(tmp$gv*scale)[1])
      } else {
        private$.addTrait(trait, popVar(tmp$bv*scale)[1], var[i])
      }
    }
    
    invisible(self)
  }
)