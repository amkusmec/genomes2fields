### RUN ONCE

### Creates a custom annotation package for calculating the information content
### of GO terms.

library(tidyverse)


ortho <- read_csv("~/anno/at_orthologs.csv") %>%
  rename(GID = `Gene stable ID`, 
         ATGID = `Arabidopsis thaliana gene stable ID`) %>%
  filter(!is.na(ATGID))

go <- read_delim("~/anno/ATH_GO_SLIM2.txt", delim = "\t", progress = FALSE) %>%
  select(Locus.name, GO.ID, GO.term, Ontology) %>%
  distinct(Locus.name, GO.ID, .keep_all = TRUE) %>%
  filter(Ontology == "P") %>%
  select(-Ontology) %>%
  rename(ATGID = Locus.name, GO = GO.ID, NAME = GO.term)

anno <- inner_join(ortho, go, by = "ATGID") %>%
  select(-ATGID, -NAME) %>%
  mutate(EVIDENCE = "IEA")

gff <- read_delim("~/anno/ZmB73_5b_FGS.gff", comment = "#", delim = "\t", 
                  progress = FALSE, col_names = FALSE) %>%
  filter(X3 == "gene", !is.na(X1)) %>%
  select(X1, X4, X5, X9) %>%
  rename(CHROMOSOME = X1, START = X4, END = X5, INFO = X9) %>%
  arrange(CHROMOSOME, START) %>%
  separate(INFO, c("GID", "NAME", "BIOTYPE"), sep = ";", remove = TRUE) %>%
  filter(str_detect(BIOTYPE, "protein")) %>%
  mutate(GID = str_remove(GID, "ID=")) %>%
  select(GID, CHROMOSOME, START, END)

all_genes <- intersect(gff$GID, anno$GID)

anno <- filter(anno, GID %in% all_genes) %>%
  distinct(GID, GO, NAME, .keep_all = TRUE)
gff <- gff %>%
  filter(GID %in% all_genes) %>%
  distinct(GID, CHROMOSOME, START, END)

library(AnnotationForge)
makeOrgPackage(gene_info = as.data.frame(gff), 
               go = as.data.frame(anno), 
               version = "0.1", maintainer = "Aaron Kusmec <amkusmec@iastate.edu>", 
               author = "Aaron Kusmec <amkusmec@iastate.edu>", 
               outputDir = ".", tax_id = "381124", 
               genus = "Zea", species = "mays", 
               goTable = "go")

install.packages("./org.Zmays.eg.db", repos = NULL)
