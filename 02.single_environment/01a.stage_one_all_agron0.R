library(tidyverse)
library(sommer)

source("src/yield_stage_one.R")
source("src/get_mmer_results.R")


# Identify variables present for each site --------------------------------
yield <- read_rds("data/phenotype/yield_agron0.rds") %>%
  filter(Replicate > 0 | is.na(Replicate))
# TXH2_2015 has some entries marked as replicate "0"
variables <- yield %>%
  group_by(Environment, Year) %>%
  summarise(Replicate = !any(is.na(Replicate)), 
            Block = !any(is.na(Block)), 
            Row = !any(is.na(Row)), 
            Column = !any(is.na(Column)), 
            Stand = !any(is.na(Stand)), 
            StalkLodging = !any(is.na(StalkLodging)), 
            RootLodging = !any(is.na(RootLodging))) %>%
  ungroup()
variables$RootLodging[54] <- FALSE # OHH1_2015 remove RootLodging from consideration


# Compute stage one estimates of hybrid effects ---------------------------
res <- list()
for (i in c(46:nrow(variables))) {
  r <- variables[i, ]
  temp <- filter(yield, Environment == r$Environment[1], Year == r$Year[1]) %>%
    mutate(Replicate = as.character(Replicate), 
           Block = paste(Replicate, Block, sep = "_"), 
           Rowf = as.character(Row), 
           Colf = as.character(Column))
  cat(temp$Site[1], "\n")
  res[[temp$Site[1]]] <- stage_one(temp, r)
}

### Some manual modifications
###  - (33) MNH1_2015
###  - (42) NEH1_2015
###  - (45) NEH3_2015

### (33) MNH1_2015 = keep only Colf
i <- 33; r <- variables[i, ]
temp <- filter(yield, Environment == r$Environment[1], Year == r$Year[1]) %>%
  mutate(Replicate = as.character(Replicate), 
         Block = paste(Replicate, Block, sep = "_"), 
         Rowf = as.character(Row), 
         Colf = as.character(Column))
# Model with Rowf and Colf only is singular
base_model1 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Rowf) + vs(spl2D(Row, Column)), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model2 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Colf) + vs(spl2D(Row, Column)), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model1$BIC; base_model2$BIC

base_model3 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Colf), rcov = ~ units, 
                    data = temp, verbose = FALSE)
base_model4 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(spl2D(Row, Column)), rcov = ~ units, 
                    data = temp, verbose = FALSE)
base_model2$BIC; base_model3$BIC; base_model4$BIC

base_model5 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block), rcov = ~ units, 
                    data = temp, verbose = FALSE)
base_model3$BIC; base_model5$BIC

res[[temp$Site[1]]] <- get_results(base_model3, temp)

### (42) NEH1_2015 = keep Rowf and Colf
i <- 42; r <- variables[i, ]
temp <- filter(yield, Environment == r$Environment[1], Year == r$Year[1]) %>%
  mutate(Replicate = as.character(Replicate), 
         Block = paste(Replicate, Block, sep = "_"), 
         Rowf = as.character(Row), 
         Colf = as.character(Column))
base_model1 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Rowf) + vs(Colf), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model2 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Rowf) + vs(spl2D(Row, Column)), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model3 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Colf) + vs(spl2D(Row, Column)), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model1$BIC; base_model3$BIC

base_model4 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Rowf), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model5 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Colf), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model1$BIC; base_model5$BIC

base_model6 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block), rcov = ~ units, 
                    data = temp, verbose = FALSE)
base_model1$BIC; base_model6$BIC

res[[temp$Site[1]]] <- get_results(base_model1, temp)

### (45) NEH3_2015
i <- 45; r <- variables[i, ]
temp <- filter(yield, Environment == r$Environment[1], Year == r$Year[1]) %>%
  mutate(Replicate = as.character(Replicate), 
         Block = paste(Replicate, Block, sep = "_"), 
         Rowf = as.character(Row), 
         Colf = as.character(Column))
base_model1 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Rowf) + vs(Colf), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model2 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Rowf) + vs(spl2D(Row, Column)), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model3 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Colf) + vs(spl2D(Row, Column)), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model1$BIC; base_model2$BIC

base_model4 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(Rowf), rcov = ~ units, 
                    data = temp, verbose = FALSE)
base_model5 <- mmer(Yield ~ 1 + PedigreeNew + Replicate, 
                    random = ~ vs(Block) + vs(spl2D(Row, Column)), 
                    rcov = ~ units, data = temp, verbose = FALSE)
base_model2$BIC; base_model4$BIC; base_model5$BIC

res[[temp$Site[1]]] <- get_results(base_model2, temp)


# Save the results --------------------------------------------------------
res <- res[order(names(res))]
write_rds(res, "data/phenotype/yield_stage_one_all_agron0.rds")
