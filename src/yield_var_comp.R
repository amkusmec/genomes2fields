### Utility function to conduct model selection for variance decomposition

require(tidyverse)
require(sommer)

var_decomp <- function(temp, r) {
  # Construct the basic fixed formula
  if (r$Replicate[1]) {
    base_fix <- "Yield ~ 1 + Replicate"
  } else {
    base_fix <- "Yield ~ 1"
  }
  
  # Construct the basic random formula
  if (r$Block[1]) {
    base_ran <- "~ vs(PedigreeNew) + vs(Block)"
  } else {
    base_ran <- "~ vs(PedigreeNew)"
  }
  
  # Identify possible fixed effects
  fterms <- c()
  if (r$Stand[1]) fterms <- c(fterms, "Stand")
  if (r$StalkLodging[1]) fterms <- c(fterms, "StalkLodging")
  if (r$RootLodging[1]) fterms <- c(fterms, "RootLodging")
  
  # Identify possible random effects
  rterms <- c()
  if (r$Row[1]) rterms <- c(rterms, "vs(Rowf)")
  if (r$Column[1]) rterms <- c(rterms, "vs(Colf)")
  if (r$Row[1] & r$Column[1]) rterms <- c(rterms, "vs(spl2D(Row, Column))")
  
  # Calculate the base model
  if (is.null(fterms) & is.null(rterms)) {
    fix_form <- base_fix
    ran_form <- base_ran
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  } else if (is.null(fterms) & !is.null(rterms)) {
    fix_form <- base_fix
    ran_form <- paste(c(base_ran, rterms), collapse = " + ")
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  } else if (!is.null(fterms) & is.null(rterms)) {
    fix_form <- paste(c(base_fix, fterms), collapse = " + ")
    ran_form <- base_ran
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  } else if (!is.null(fterms) & !is.null(rterms)) {
    fix_form <- paste(c(base_fix, fterms), collapse = " + ")
    ran_form <- paste(c(base_ran, rterms), collapse = " + ")
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  }
  
  # Random effects selection
  if (!is.null(rterms)) {
    repeat {
      if (length(rterms) == 1) {
        rmodels <- list()
        if (base_ran == "~ ") {
          rmodels[[1]] <- lm(as.formula(fix_form), data = temp)
        } else {
          ran_form <- base_ran
          rmodels[[1]] <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                               rcov = ~ vs(units), data = temp, verbose = FALSE)
        }
      } else {
        rmodels <- lapply(seq_along(rterms), function(i) {
          if (base_ran == "~ ") {
            ran_form <- paste0("~ ", paste(rterms[-i], collapse = " + "))
          } else {
            ran_form <- paste(base_ran, paste(rterms[-i], collapse = " + "), sep = " + ")
          }
          mmer(as.formula(fix_form), random = as.formula(ran_form), 
               rcov = ~ vs(units), data = temp, verbose = FALSE)
        })
      }
      
      bic <- c(base_model$BIC, sapply(rmodels, function(m) {
        if (class(m) == "mmer") {
          m$BIC
        } else {
          BIC(m)
        }
      }))
      idx <- which.min(bic)
      # idx <- max(which(bic == min(bic)))
      
      if (idx == 1) break
      
      cat("Drop", rterms[idx - 1], "\n")
      rterms <- rterms[-(idx - 1)]
      base_model <- rmodels[[idx - 1]]
      
      if (length(rterms) == 0) break
    }
    
    if (length(rterms) == 0) {
      ran_form <- base_ran
      base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                         rcov = ~ vs(units), data = temp, verbose = FALSE)
    } else {
      ran_form <- paste(base_ran, paste(rterms, collapse = " + "), sep = " + ")
      base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                         rcov = ~ vs(units), data = temp, verbose = FALSE)
    }
  }
  
  if (!is.null(fterms)) {
    repeat {
      if (length(fterms) == 1) {
        fmodels <- list()
        fmodels[[1]] <- mmer(as.formula(base_fix), random = as.formula(ran_form), 
                             rcov = ~ vs(units), data = temp, verbose = FALSE)
      } else {
        fmodels <- lapply(seq_along(fterms), function(i) {
          fix_form <- paste(base_fix, paste(fterms[-i], collapse = " + "), sep = " + ")
          mmer(as.formula(fix_form), random = as.formula(ran_form), 
               rcov = ~ vs(units), data = temp, verbose = FALSE)
        })
      }
      
      bic <- c(base_model$BIC, sapply(fmodels, function(m) m$BIC))
      idx <- which.min(bic)
      
      if (idx == 1) break
      
      cat("Drop", fterms[idx - 1], "\n")
      fterms <- fterms[-(idx - 1)]
      base_model <- fmodels[[idx - 1]]
      
      if (length(fterms) == 0) break
    }
    
    if (length(fterms) == 0) {
      fix_form <- base_fix
    } else {
      fix_form <- paste(base_fix, paste(fterms, collapse = " + "), sep = " + ")
    }
    
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  }
  
  # # Get the unscaled variance components and return the results
  # tibble(Year = r$Year[1], 
  #        Environment = r$Environment[1], 
  #        Component = names(base_model$sigma) %>% str_replace(., "u:", ""), 
  #        Variance = unlist(base_model$sigma, use.names = FALSE))
  return(base_model)
}
