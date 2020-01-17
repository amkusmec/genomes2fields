require(tidyverse)
require(sommer)

get_results <- function(base_model, temp) {
  co <- coef(base_model) %>%
    as_tibble() %>%
    dplyr::select(-Trait) %>%
    mutate(Effect = as.character(Effect), 
           Estimate = c(Estimate[1], Estimate[2:n()] + Estimate[1]))
  
  vc <- base_model$VarBeta
  idx <- which(str_detect(co$Effect, "PedigreeNew") | co$Effect == "(Intercept)")
  vc <- vc[idx, idx]
  
  co <- filter(co, str_detect(Effect, "PedigreeNew") | Effect == "(Intercept)") %>%
    mutate(Effect = str_remove(Effect, "PedigreeNew")) %>%
    rename(PedigreeNew = Effect, BLUE = Estimate)
  co$PedigreeNew[1] <- setdiff(unique(temp$PedigreeNew), co$PedigreeNew)
  
  rownames(vc) <- co$PedigreeNew
  colnames(vc) <- co$PedigreeNew
  
  return(list(blue = co, vcov = vc, model = base_model$call))
}
