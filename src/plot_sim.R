require(tidyverse)
require(abind)

plot_sim <- function(res, scaled = FALSE, pheno = NULL) {
  vm <- lapply(res, function(x) do.call("rbind", x$genMean))
  vg <- lapply(res, function(x) {
    do.call("rbind", lapply(x$genVar, function(y) diag(y))) 
  })
  
  ngen <- nrow(vm[[1]]); npheno <- ncol(vm[[1]])
  if (scaled) {
    for (i in seq_along(vm)) {
      vm[[i]] <- sapply(1:npheno, function(j) (vm[[i]][, j] - vm[[i]][1, j])/sqrt(vg[[i]][1, j]))
      vg[[i]] <- sapply(1:npheno, function(j) vg[[i]][, j]/vg[[i]][1, j])
    }
  }
  
  vm <- abind(vm, along = 3)
  vg <- abind(vg, along = 3)
  
  if (!is.null(pheno)) {
    dimnames(vm)[[2]] <- pheno
    dimnames(vg)[[2]] <- pheno
  }
  
  vm_mean <- apply(vm, c(1, 2), mean)
  vm_sd <- apply(vm, c(1, 2), sd)
  vg_mean <- apply(vg, c(1, 2), mean)
  vg_sd <- apply(vg, c(1, 2), sd)
  
  d <- rbind(vm_mean, vm_sd, vg_mean, vg_sd) %>%
    as_tibble(.name_repair = "minimal") %>%
    mutate(Generation = rep(0:(ngen - 1), times = 4), 
           Measure = rep(c("Genetic Gain", "Genetic Variance"), each = ngen*2),
           Statistic = rep(rep(c("Mean", "SD"), each = ngen), times = 2)) %>%
    gather(Phenotype, Value, -Generation, -Measure, -Statistic) %>%
    spread(Statistic, Value)
  
  if (scaled) {
    d <- d %>% 
      mutate(Measure = paste("Scaled", Measure))
  }
  
  p <- ggplot(d, aes(x = Generation, y = Mean, colour = Phenotype)) +
    theme_bw() + geom_line(aes(group = Phenotype)) + 
    geom_pointrange(aes(ymin = Mean - SD, ymax = Mean + SD)) + 
    geom_hline(yintercept = 0, linetype = 2) + 
    scale_colour_brewer(type = "qual", palette = "Paired") + 
    facet_wrap(~ Measure, scales = "free_y") + 
    labs(x = "Generation", y = "")
  return(p)
}
