library(tidyverse)


cl_terms <- read_csv("data/gemma/closest_vs_ld_go_terms.csv")
cg_terms <- read_csv("data/gemma/closest_vs_genome_go_terms.csv")
lg_terms <- read_csv("data/gemma/ld_vs_genome_go_terms.csv")

cl_permute <- read_rds("data/gemma/closest_vs_ld_permute.rds")
cg_permute <- read_rds("data/gemma/closest_vs_genome_permute.rds")
lg_permute <- read_rds("data/gemma/ld_vs_genome_permute.rds")


cl_enrich <- cl_permute %>%
  map(function(df) filter(df, p_val <= 0.05) %>% pull(GO)) %>%
  unlist(use.names = FALSE) %>%
  enframe() %>%
  rename(GO = value) %>%
  count(GO) %>%
  mutate(p = n/1e4)
cg_enrich <- cg_permute %>%
  map(function(df) filter(df, p_val <= 0.05) %>% pull(GO)) %>%
  unlist(use.names = FALSE) %>%
  enframe() %>%
  rename(GO = value) %>%
  count(GO) %>%
  mutate(p = n/1e4)
lg_enrich <- lg_permute %>%
  map(function(df) filter(df, p_val <= 0.05) %>% pull(GO)) %>%
  unlist(use.names = FALSE) %>%
  enframe() %>%
  rename(GO = value) %>%
  count(GO) %>%
  mutate(p = n/1e4)


cl_terms <- left_join(cl_terms, cl_enrich, by = "GO") %>%
  rename(perm = p)
cg_terms <- left_join(cg_terms, cg_enrich, by = "GO") %>%
  rename(perm = p)
lg_terms <- left_join(lg_terms, lg_enrich, by = "GO") %>%
  rename(perm = p)

cl <- cl_terms %>% filter(perm <= 0.05) %>% pull(Name)
cg <- cg_terms %>% filter(perm <= 0.05) %>% pull(Name)
lg <- lg_terms %>% filter(perm <= 0.05) %>% pull(Name)

write_csv(cl_terms, "data/gemma/closest_vs_ld_perm.csv")
write_csv(cg_terms, "data/gemma/closest_vs_genome_perm.csv")
write_csv(lg_terms, "data/gemma/ld_vs_genome_perm.csv")


cl_terms %>%
  select(GO, Name, perm, white_draws, white, black, draws) %>%
  rename(`GO ID` = GO, `GO Term` = Name, `Permutation p-value` = perm, 
         Successes = white_draws, `Max. Successes` = white, Failures = black, 
         Trials = draws) %>%
  write_csv("data/gemma/table_s4.csv")

cg_terms %>%
  select(GO, Name, perm, white_draws, white, black, draws) %>%
  rename(`GO ID` = GO, `GO Terms` = Name, `Permutation p-value` = perm, 
         Successes = white_draws, `Max. Successes` = white, Failures = black, 
         Trials = draws) %>%
  write_csv("data/gemma/table_s5.csv")
