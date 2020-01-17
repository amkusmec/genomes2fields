require(tidyverse)
require(purrrlyr)


# Load data files ---------------------------------------------------------
ortho <- read_csv("~/anno/at_orthologs.csv") %>%
  rename(GID = `Gene stable ID`, 
         ATGID = `Arabidopsis thaliana gene stable ID`) %>%
  filter(!is.na(ATGID))

go <- read_delim("~/anno/ATH_GO_SLIM2.txt", delim = "\t", progress = FALSE) %>%
  select(Locus.name, GO.ID, GO.term, Ontology) %>%
  distinct(Locus.name, GO.ID, .keep_all = TRUE) %>%
  filter(Ontology == "P") %>%
  select(-Ontology) %>%
  rename(ATGID = Locus.name, GO = GO.ID, Name = GO.term)

anno <- inner_join(ortho, go, by = "ATGID")


# GO term enrichment ------------------------------------------------------
go_enrich <- function(cand, background = anno) {
  tibble(Gene = cand) %>%
    inner_join(., anno, by = c("Gene" = "GID")) %>%
    count(GO) %>%
    dplyr::rename(white_draws = n) %>%
    by_row(function(r) {
      m <- length(unique(background$GID[background$GO == r$GO[1]]))
      n <- length(unique(background$GID)) - m
      tibble(white = m, black = n, draws = length(cand), 
             p_val = phyper(r$white_draws[1] - 1, m, n, draws, lower.tail = FALSE))
    }, .collate = "rows") %>%
    mutate(q_val = p.adjust(p_val, method = "fdr")) %>%
    left_join(., distinct(anno, GO, Name), by = "GO")
}
