### Utility function to conduct stage one analysis for grain yield

require(tidyverse)
require(sommer)

stage_one <- function(temp, r) {
  # Construct the basic fixed formula
  if (r$Replicate[1]) {
    # base_fix <- "Yield ~ 0 + PedigreeNew + Replicate"
    base_fix <- "Yield ~ 1 + PedigreeNew + Replicate"
  } else {
    # base_fix <- "Yield ~ 0 + PedigreeNew"
    base_fix <- "Yield ~ 1 + PedigreeNew"
  }
  
  # Construct the basic random formula
  if (r$Block[1]) {
    base_ran <- "~ vs(Block)"
  } else {
    base_ran <- "~ "
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
  if (is.null(fterms) & base_ran == "~ " & is.null(rterms)) {
    fix_form <- base_fix
    base_model <- lm(as.formula(fix_form), data = temp)
  } else if (!is.null(fterms) & base_ran == "~ " & is.null(rterms)) {
    fix_form <- paste(base_fix, paste(fterms, collapse = " + "), sep = " + ")
    base_model <- lm(as.formula(fix_form), data = temp)
  } else if (is.null(fterms) & base_ran != "~ " & is.null(rterms)) {
    fix_form <- base_fix
    ran_form <- base_ran
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  } else if (is.null(fterms) & base_ran == "~ " & !is.null(rterms)) {
    fix_form <- base_fix
    ran_form <- paste0("~ ", paste(rterms, collapse = " + "))
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  } else if (!is.null(fterms) & base_ran != "~ " & is.null(rterms)) {
    fix_form <- paste(base_fix, paste(fterms, collapse = " + "), sep = " + ")
    ran_form <- base_ran
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  } else if (!is.null(fterms) & base_ran == "~ " & !is.null(rterms)) {
    fix_form <- paste(base_fix, paste(fterms, collapse = " + "), sep = " + ")
    ran_form <- paste0("~ ", paste(rterms, collapse = " + "))
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  } else if (is.null(fterms) & base_ran != "~ " & !is.null(rterms)) {
    fix_form <- base_fix
    ran_form <- paste(base_ran, paste(rterms, collapse = " + "), sep = " + ")
    base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                       rcov = ~ vs(units), data = temp, verbose = FALSE)
  } else {
    fix_form <- paste(base_fix, paste(fterms, collapse = " + "), sep = " + ")
    ran_form <- paste(base_ran, paste(rterms, collapse = " + "), sep = " + ")
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
      
      if (idx == 1) break
      
      cat("Drop", rterms[idx - 1], "\n")
      rterms <- rterms[-(idx - 1)]
      base_model <- rmodels[[idx - 1]]
      
      if (length(rterms) == 0) break
    }
    
    if (base_ran == "~ " & length(rterms) == 0) {
      base_model <- lm(as.formula(fix_form), data = temp)
    } else if (base_ran != "~ " & length(rterms) == 0) {
      ran_form <- base_ran
      base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                         rcov = ~ vs(units), data = temp, verbose = FALSE)
    } else if (base_ran == "~ " & length(rterms) > 0) {
      ran_form <- paste0("~ ", paste(rterms, collapse = " + "))
      base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                         rcov = ~ vs(units), data = temp, verbose = FALSE)
    } else {
      ran_form <- paste(base_ran, paste(rterms, collapse = " + "), sep = " + ")
      base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                         rcov = ~ vs(units), data = temp, verbose = FALSE)
    }
  }
  
  if (!is.null(fterms)) {
    if (class(base_model) == "mmer") { # Fixed effects selection with random terms
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
    } else { # Fixed effects selection without random terms
      repeat {
        if (length(fterms) == 1) {
          fmodels <- list()
          fmodels[[1]] <- lm(as.formula(base_fix), data = temp)
        } else {
          fmodels <- lapply(seq_along(fterms), function(i) {
            fix_form <- paste(base_fix, paste(fterms[-i], collapse = " + "), sep = " + ")
            lm(as.formula(fix_form), data = temp)
          })
        }
        
        bic <- c(BIC(base_model), sapply(fmodels, BIC))
        idx <- which.min(bic)
        
        if (idx == 1) break
        
        cat("Drop", fterms[idx - 1], "\n")
        fterms <- fterms[-(idx - 1)]
        base_model <- fmodels[[idx - 1]]
        
        if (length(fterms) == 0) break
      }
    }
    
    if (length(fterms) == 0) {
      fix_form <- base_fix
    } else {
      fix_form <- paste(base_fix, paste(fterms, collapse = " + "), sep = " + ")
    }
    
    if (class(base_model) == "mmer") {
      base_model <- mmer(as.formula(fix_form), random = as.formula(ran_form), 
                         rcov = ~ vs(units), data = temp, verbose = FALSE)
    } else {
      base_model <- lm(as.formula(fix_form), data = temp)
    }
  }
  
  # Get the fixed genotype effects and variance-covariance matrices
  if (class(base_model) == "mmer") {
    # co <- coef(base_model) %>%
    #   as_tibble() %>%
    #   dplyr::select(-Trait)
    
    co <- coef(base_model) %>%
      as_tibble() %>%
      dplyr::select(-Trait) %>%
      mutate(Effect = as.character(Effect), 
             Estimate = c(Estimate[1], Estimate[2:n()] + Estimate[1]))
    
    vc <- base_model$VarBeta
    # idx <- which(str_detect(co$Effect, "PedigreeNew"))
    idx <- which(str_detect(co$Effect, "PedigreeNew") | co$Effect == "(Intercept)")
    vc <- vc[idx, idx]
    
    # co <- filter(co, str_detect(Effect, "PedigreeNew")) %>%
    #   mutate(Effect = str_replace(Effect, "PedigreeNew", "")) %>%
    #   rename(PedigreeNew = Effect, BLUE = Estimate)
    
    co <- filter(co, str_detect(Effect, "PedigreeNew") | Effect == "(Intercept)") %>%
      mutate(Effect = str_remove(Effect, "PedigreeNew")) %>%
      rename(PedigreeNew = Effect, BLUE = Estimate)
    co$PedigreeNew[1] <- setdiff(unique(temp$PedigreeNew), co$PedigreeNew)
    
    rownames(vc) <- co$PedigreeNew
    colnames(vc) <- co$PedigreeNew
  } else {
    co <- coef(base_model)
    co <- tibble(PedigreeNew = names(co), 
                 BLUE = co) %>%
      mutate(PedigreeNew = str_replace(PedigreeNew, "PedigreeNew", ""), 
             BLUE = c(BLUE[1], BLUE[2:n()] + BLUE[1]))
    
    vc <- vcov(base_model)
    idx <- which(str_detect(co$PedigreeNew, "/") | co$PedigreeNew == "(Intercept)")
    
    co <- filter(co, str_detect(PedigreeNew, "/") | PedigreeNew == "(Intercept)")
    co$PedigreeNew[1] <- setdiff(unique(temp$PedigreeNew), co$PedigreeNew)
    
    vc <- vc[idx, idx]
    rownames(vc) <- colnames(vc) <- co$PedigreeNew
  }
  
  # Return the results
  return(list(blue = co, vcov = vc, model = base_model$call))
}
