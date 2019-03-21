library(tidyverse)
library(sommer)

source("src/yield_stage_one.R")


# Identify variables present for each site --------------------------------
yield <- read_rds("data/phenotype/yield_agron40.rds")
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


# Compute stage one estimates of hybrid effects ---------------------------
res <- list()
for (i in c(1:41, 43:44, 46:53, 55:58, 60:nrow(variables))) {
  r <- variables[i, ]
  temp <- filter(yield, Environment == r$Environment[1], Year == r$Year[1]) %>%
    mutate(Replicate = as.character(Replicate), 
           Block = paste(Replicate, Block, sep = "_"), 
           Rowf = as.character(Row), 
           Colf = as.character(Column))
  cat(temp$Site[1], "\n")
  res[[temp$Site[1]]] <- stage_one(temp, r)
}

write_rds(res, "data/phenotype/yield_stage_one_all_agron40.rds")

### Some manual modifications
###  - (42) NEH1_2015 = remove Rowf
###  - (45) NEH3_2015 = remove Colf and StalkLodging
###  - (54) OHH1_2015 = remove RootLodging from consideration
###  - (59) TXH2_2017 = remove Replicate == "0"
