library(tidyverse)


# Collect potentially interesting structure results -----------------------
snps <- read_rds("data/gbs/add_snps.rds")
taxa <- rownames(snps$GD)
testers <- c("LH198", "LH82", "LH185", "PHZ51", "PHB47", "LH195", "3IIH6")

st2 <- read.table("data/structure/g2f_hyb_rep1.2.meanQ", header = FALSE) %>%
  as_tibble() %>%
  mutate(Taxa = taxa, 
         NAM = str_detect(Taxa, "Z[0-9]{3}E[0-9]{4}")) %>%
  separate(Taxa, c("Parent1", "Parent2"), sep = "/", remove = FALSE) %>%
  mutate(Tester = Parent2 %in% testers) %>%
  gather(Population, Proportion, V1:V2) %>%
  arrange(Taxa, desc(Proportion)) %>%
  distinct(Taxa, .keep_all = TRUE)
st5 <- read.table("data/structure/g2f_hyb_rep1.5.meanQ", header = FALSE) %>%
  as_tibble() %>%
  mutate(Taxa = taxa, 
         NAM = str_detect(Taxa, "Z[0-9]{3}E[0-9]{4}")) %>%
  separate(Taxa, c("Parent1", "Parent2"), sep = "/", remove = FALSE) %>%
  mutate(Tester = Parent2 %in% testers) %>%
  gather(Population, Proportion, V1:V5) %>%
  arrange(Taxa, desc(Proportion)) %>%
  distinct(Taxa, .keep_all = TRUE)
st13 <- read.table("data/structure/g2f_hyb_rep1.13.meanQ", header = FALSE) %>%
  as_tibble() %>%
  mutate(Taxa = taxa, 
         NAM = str_detect(Taxa, "Z[0-9]{3}E[0-9]{4}")) %>%
  separate(Taxa, c("Parent1", "Parent2"), sep = "/", remove = FALSE) %>%
  mutate(Tester = Parent2 %in% testers) %>%
  gather(Population, Proportion, V1:V13) %>%
  arrange(Taxa, desc(Proportion)) %>%
  distinct(Taxa, .keep_all = TRUE)


with(st2, table(NAM, Population))
with(st5, table(NAM, Population))
with(st13, table(NAM, Population))

filter(st5, Tester) %>%
  with(., table(Parent2, Population))

filter(st13, Tester) %>%
  with(., table(Parent2, Population))
