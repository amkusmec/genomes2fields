require(tidyverse)
require(purrrlyr)


# Load data files ---------------------------------------------------------
xref <- read_tsv("~/anno/gene_model_xref_v2.txt", skip = 4) %>%
  select(v2_gene_model:v2_chr, v3_gene_model) %>%
  mutate(v2_chr = str_remove(v2_chr, "chr") %>% as.integer())
gamer <- read_tsv("~/anno/maize_v3.agg.nr.gaf", skip = 1) %>%
  filter(!is.na(db_object_type), aspect == "P") %>%
  select(db_object_symbol, term_accession) %>%
  dplyr::rename(Gene = db_object_symbol, GO = term_accession) %>%
  inner_join(., xref, by = c("Gene" = "v3_gene_model"))
obo <- read_lines("~/anno/go-basic.obo")
obo <- split(obo, cumsum(str_detect(obo, "\\[Term\\]")))
obo <- obo[sapply(obo, function(o) any(str_detect(o, "^namespace\\: biological_process")))]
obo <- obo %>%
  map_df(function(o) {
    tibble(GO = c(str_remove(o[2], "^id\\: "), 
                  str_remove(o[str_detect(o, "^alt-id\\: ")], "^alt-id\\: ")), 
           Name = str_remove(o[str_detect(o, "^name\\: ")], "^name\\: "))
  })


# GO term enrichment ------------------------------------------------------
go_enrich <- function(cand, background = gamer) {
  tibble(Gene = cand) %>%
    inner_join(., gamer, by = c("Gene" = "v2_gene_model")) %>%
    count(GO) %>%
    dplyr::rename(white_draws = n) %>%
    by_row(function(r) {
      m <- length(unique(background$v2_gene_model[background$GO == r$GO[1]]))
      n <- length(unique(background$v2_gene_model)) - m
      tibble(white = m, 
             black = n, 
             draws = length(cand), 
             p_val = phyper(r$white_draws[1] - 1, m, n, draws, lower.tail = FALSE))
    }, .collate = "rows") %>%
    mutate(q_val = p.adjust(p_val, method = "fdr")) %>%
    left_join(., obo, by = "GO")
}
