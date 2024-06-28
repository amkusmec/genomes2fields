library(tidyverse)

cl <- read_csv("data/gemma/closest_vs_ld_perm.csv") %>%
  filter(perm <= 0.05) %>%
  mutate(Background = "LD") %>%
  select(Background, GO, Name, perm, white_draws, white, black, draws) %>%
  rename(`GO ID` = GO, `GO Terms` = Name, `Permutation p-value` = perm, 
         Successes = white_draws, `Max. Successes` = white, Failures = black, 
         Trials = draws)
cg <- read_csv("data/gemma/closest_vs_genome_perm.csv") %>%
  filter(perm <= 0.05) %>%
  mutate(Background = "Genome") %>%
  select(Background, GO, Name, perm, white_draws, white, black, draws) %>%
  rename(`GO ID` = GO, `GO Terms` = Name, `Permutation p-value` = perm, 
         Successes = white_draws, `Max. Successes` = white, Failures = black, 
         Trials = draws)

write_csv(cl, "data/tableS6_GO_terms_LD.csv")
write_csv(cg, "data/tableS7_GO_terms_genome.csv")
