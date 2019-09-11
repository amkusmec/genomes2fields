library(tidyverse)
library(readxl)


# Prepare the NAM markers -------------------------------------------------
# Removal of markers was determined interactively and then hard-coded into
# this script.
map <- read_xls("data/gbs/NAM/NAM_map_and_genos-121025/NAM_map_20080419.xls") %>%
  select(marker, ch, cumulative) %>%
  mutate(marker = as.character(marker))
positions <- read_xlsx("data/gbs/NAM/NAM_map_and_genos-121025/NAM_1144SNPs_AGPv2_positions.xlsx") %>%
  select(SNP_NAME:AGPv2_pos) %>%
  mutate(SNP_NAME = as.character(SNP_NAME),
         AGPv2_Chr = as.numeric(as.character(AGPv2_Chr)),
         AGPv2_pos = as.numeric(as.character(AGPv2_pos)))
map <- inner_join(map, positions, by = c("marker" = "SNP_NAME")) %>%
  select(marker, AGPv2_Chr, AGPv2_pos, cumulative) %>%
  rename(Marker = marker, Chr = AGPv2_Chr, Physical = AGPv2_pos, 
         Genetic = cumulative) %>%
  filter(!is.na(Physical)) %>%
  arrange(Chr, Physical)

map <- map[-c(243, 279, 280, 319, 450, 460, 503, 518, 533, 698, 702, 746, 
              748, 838, 887, 961), ]
map <- map[-873, ]


# Function to interpolate genetic positions -------------------------------
# Peforms linear interpolation of genetic positions between two flanking markers
p2g <- function(map, markers) {
  # Use the first interval to estimate the equation for markers prior to the
  # first mapped marker
  m <- (map$Genetic[2] - map$Genetic[1])/(map$Physical[2] - map$Physical[1])
  offset <- map$Genetic[2] - m*map$Physical[2]
  markers <- markers %>%
    mutate(Genetic = if_else(Physical <= map$Physical[1], 
                             m*Physical + offset + abs(offset), Genetic))
  
  # Iterate through each interval
  offset <- abs(offset)
  for (j in 1:(nrow(map) - 1)) {
    m <- (map$Genetic[j + 1] - map$Genetic[j])/(map$Physical[j + 1] - map$Genetic[j])
    b <- map$Genetic[j + 1] - m*map$Physical[j + 1]
    markers <- markers %>%
      mutate(Genetic = if_else(Physical >= map$Physical[j] & Physical <= map$Physical[j + 1],
                               m*Physical + b + offset, Genetic))
  }
  
  # Use the final interval to estimate the equation for markers beyond the
  # last mapped marker.
  m <- (map$Genetic[nrow(map)] - map$Genetic[nrow(map) - 1])/
    (map$Physical[nrow(map)] - map$Physical[nrow(map) - 1])
  b <- map$Genetic[nrow(map)] - m*map$Physical[nrow(map)]
  markers <- markers %>%
    mutate(Genetic = if_else(Physical >= map$Physical[nrow(map)],
                             m*Physical + b + offset, Genetic))
  
  return(markers)
}


# Compute the new map -----------------------------------------------------
markers <- read_rds("data/gbs/add_snps.rds")$GM %>%
  select(SNP:Position) %>%
  rename(Marker = SNP, Chr = Chromosome, Physical = Position) %>%
  mutate(Genetic = 0)
markers <- map2(.x = split(map, map$Chr), 
                .y = split(markers, markers$Chr), p2g) %>%
  bind_rows() %>%
  mutate(Genetic = Genetic/100)

write_csv(markers, "data/gbs/sim_map.csv")

genMap <- map %>% split(., .$Chr) %>%
  map(function(df) pull(df, Genetic))
names(genMap) <- NULL
write_rds(genMap, "data/gbs/sim_map.rds")
