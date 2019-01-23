library(tidyverse)
library(lubridate)


# Download files ----------------------------------------------------------
### 2017 data is currently downloaded from a private cooperator link

### See 01.munge_metadata.R for information about configuring icommands to
### download data from CyVerse

### This line will need to change based on where you installed icommands.
Sys.setenv("PATH" = paste("/mnt/01/amkusmec/bin/icommands", Sys.getenv("PATH"), sep = ":"))

setwd("data/phenotype")
system("iinit")
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Nov_2016_V.3/", 
              "a._2014_hybrid_phenotypic_data/_g2f_2014_hybrid_data_description.txt"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Nov_2016_V.3/", 
              "a._2014_hybrid_phenotypic_data/g2f_2014_hybrid_no_outliers.csv"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Mar_2017/", 
              "a._2015_hybrid_phenotypic_data/_g2f_2015_hybrid_data_description.txt"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "Carolyn_Lawrence_Dill_G2F_Mar_2017/", 
              "a._2015_hybrid_phenotypic_data/g2f_2015_hybrid_data_no_outliers.csv"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "GenomesToFields_G2F_2016_Data_Mar_2018/", 
              "a._2016_hybrid_phenotypic_data/_g2f_2016_hybrid_data_description.txt"))
system(paste0("iget /iplant/home/shared/commons_repo/curated/", 
              "GenomesToFields_G2F_2016_Data_Mar_2018/", 
              "a._2016_hybrid_phenotypic_data/g2f_2016_hybrid_data_no_outliers.csv"))
system("iexit")
setwd("../..")


# Process the yield data --------------------------------------------------
yield_2014 <- read_csv("data/phenotype/g2f_2014_hybrid_no_outliers.csv", 
                       col_types = "icccicciiiiiiidcccciiiiiiiddddcccidddi") %>%
  filter(is.na(`Plot Discarded`),
         is.na(Filler),
         !str_detect(Pedigree, "\\?\\?\\?TX205"),
         str_detect(Pedigree, "/")) %>%
  mutate(`Date Planted` = mdy(`Date Planted`),
         `Date Harvested` = mdy(`Date Harvested`),
         `Anthesis [date]` = mdy(`Anthesis [date]`)) %>%
  select(`Field-Location`, Pedigree, Experiment, Rep:Pass,
         `Date Planted`:`Stand Count [plants]`, 
         `Plant height [cm]`:`Grain yield [bu/A]`) %>%
  rename(Environment = `Field-Location`, Set = Experiment, Replicate = Rep, 
         Row = Range, Column = Pass, Planted = `Date Planted`,
         Harvested = `Date Harvested`, Anthesis = `Anthesis [date]`,
         Stand = `Stand Count [plants]`, PlantHeight = `Plant height [cm]`, 
         EarHeight = `Ear height [cm]`, RootLodging = `Root lodging [plants]`, 
         StalkLodging = `Stalk lodging [plants]`, GrainMoisture = `Grain Moisture [percent]`, 
         TestWeight = `Test Weight [lbs/bu]`, PlotWeight = `Plot Weight [lbs]`, 
         Yield = `Grain yield [bu/A]`) %>%
  mutate(Environment = str_replace(Environment, "[a-z]$", "")) %>%
  filter(!(Environment %in% c("GAH1", "ONH1", "ONH2") & PlotWeight < 3)) %>%
  filter(!is.na(Anthesis)) %>%
  filter(Yield < 340) %>%
  mutate(Year = 2014) %>%
  select(Year, Environment:Block, Row:Anthesis, Stand, StalkLodging, RootLodging, 
         Yield:GrainMoisture, PlantHeight, EarHeight)

yield_2015 <- read_csv("data/phenotype/g2f_2015_hybrid_data_no_outliers.csv", 
                       col_types = "icccicciiiiidcccciiiiiiiddddccciiiiiiiiiiiiiiiiii") %>%
  filter(is.na(`Plot Discarded`),
         is.na(Filler),
         str_detect(Pedigree, "/")) %>%
  mutate(`Date Planted` = mdy(`Date Planted`),
         `Date Harvested` = mdy(`Date Harvested`),
         `Anthesis [date]` = mdy(`Anthesis [date]`)) %>%
  select(`Field-Location`, Pedigree:Pass, `Date Planted`:`Stand Count [plants]`,
         `Plant height [cm]`:`Grain yield [bu/acre]`) %>%
  rename(Environment = `Field-Location`, Row = Range, Column = Pass, 
         Planted = `Date Planted`, Harvested = `Date Harvested`,
         Anthesis = `Anthesis [date]`, Stand = `Stand Count [plants]`, 
         PlantHeight = `Plant height [cm]`, EarHeight = `Ear height [cm]`, 
         RootLodging = `Root lodging [plants]`, StalkLodging = `Stalk lodging [plants]`,
         GrainMoisture = `Grain Moisture [percent]`, TestWeight = `Test weight [lbs]`,
         PlotWeight = `Plot Weight [lbs]`, Yield = `Grain yield [bu/acre]`) %>%
  filter(!(Environment %in% c("MOH1", "NYH1", "NYH2", "KSH1", "ILH2"))) %>%
  filter(!is.na(Anthesis)) %>%
  mutate(Year = 2015, Set = as.integer(NA)) %>%
  select(Year, Environment, Pedigree, Set, Replicate:Block, Row:Anthesis, Stand, 
         StalkLodging, RootLodging, Yield:GrainMoisture, PlantHeight, EarHeight)

### N.B. There are no cooperator comments in the README with respect to site-
### specific issues with particular trials or exclusion criteria.
yield_2016 <- read_csv("data/phenotype/g2f_2016_hybrid_data_no_outliers.csv") %>%
  filter(is.na(`Plot Discarded`), is.na(Filler), 
         str_detect(Pedigree, "/")) %>%
  mutate(`Date Plot Planted` = mdy(`Date Plot Planted`),
         `Date Plot Harvested` = mdy(`Date Plot Harvested`),
         `Anthesis [date]` = mdy(`Anthesis [date]`),
         Year = 2016, Block = as.integer(NA), Set = as.integer(NA)) %>%
  select(`Field-Location`, Pedigree, Replicate:`Silking [date]`, 
         `Plant Height [cm]`:`Grain Yield [bu/acre]`, Year:Set) %>%
  rename(Environment = `Field-Location`, Row = Range, Column = Pass,
         Planted = `Date Plot Planted`, Harvested = `Date Plot Harvested`,
         Anthesis = `Anthesis [date]`, PlantHeight = `Plant Height [cm]`, 
         EarHeight = `Ear Height [cm]`, Stand = `Stand Count [plants]`, 
         RootLodging = `Root Lodging [plants]`, StalkLodging = `Stalk Lodging [plants]`, 
         GrainMoisture = `Grain Moisture [percent]`, 
         TestWeight = `Test Weight [lbs/bu]`, PlotWeight = `Plot Weight [lbs]`, 
         Yield = `Grain Yield [bu/acre]`) %>%
  filter(!is.na(Anthesis), !is.na(Harvested)) %>%
  select(Year, Environment:Pedigree, Set, Replicate, Block, Row:Anthesis, 
         Stand, StalkLodging:RootLodging, Yield:GrainMoisture,
         PlantHeight:EarHeight)

### N.B. There are no cooperator comments in the README with respect to site-
### specific issues with particular trials or exclusion criteria.
yield_2017 <- read_csv("data/phenotype/g2f_2017_hyrbid_data_no_outliers.csv") %>%
  filter(is.na(`Plot Discarded`), is.na(Comments), 
         !(Filler %in% c("No seed", "Planting error"))) %>%
  select(`Field-Location`, Pedigree, Replicate:`Silking [date]`, 
         `Plant Height [cm]`:`Grain yield bu/A`) %>%
  mutate(Year = 2017, Block = as.integer(NA), Set = as.integer(NA),
         `Date Planted` = mdy(`Date Planted`),
         `Date Harvested` = mdy(`Date Harvested`),
         `Anthesis [date]` = mdy(`Anthesis [date]`),
         Range = as.integer(Range), Pass = as.integer(Pass),
         `Stalk Lodging [plants]` = as.integer(`Stalk Lodging [plants]`),
         `Root Lodging [plants]` = as.integer(`Root Lodging [plants]`)) %>%
  rename(Environment = `Field-Location`, Row = Range, Column = Pass, 
         Planted = `Date Planted`, Harvested = `Date Harvested`, 
         Anthesis = `Anthesis [date]`, PlantHeight = `Plant Height [cm]`, 
         EarHeight = `Ear Height [cm]`, Stand = `Stand Count [plants]`, 
         RootLodging = `Root Lodging [plants]`, StalkLodging = `Stalk Lodging [plants]`, 
         GrainMoisture = `Grain Moisture [percent]`, TestWeight = `Test Weight [lbs/bu]`, 
         PlotWeight = `Plot Weight [lbs]`, Yield = `Grain yield bu/A`) %>%
  filter(!is.na(Anthesis)) %>%
  filter(!str_detect(Pedigree, "ex pvp"), !str_detect(Pedigree, "Maizex"), 
         !str_detect(Pedigree, "stax")) %>%
  mutate(Pedigree = str_replace_all(Pedigree, " ", "") %>%
           str_replace(., "DK3IIH6", "3IIH6"), 
         Pedigree = if_else(str_detect(Pedigree, "x") & !str_detect(Pedigree, "/"), 
                            str_replace(Pedigree, "x", "/"), Pedigree), 
         Pedigree = if_else(str_detect(Pedigree, "\\*"), 
                            str_replace(Pedigree, "\\*", "/"), Pedigree)) %>%
  filter(str_detect(Pedigree, "/")) %>%
  select(Year, Environment:Pedigree, Set, Replicate, Block, Row:Anthesis, 
         Stand, StalkLodging:RootLodging, Yield:GrainMoisture, 
         PlantHeight:EarHeight)


# Combine all data sets and save ------------------------------------------
yield <- bind_rows(yield_2014, yield_2015, yield_2016, yield_2017) %>%
  mutate(Pedigree = if_else(str_detect(Pedigree, "ARG"), 
                            "ARGENTINE_FLINTY_COMPOSITE-C(1)-37-B-B-B2-1-B25/LH195", 
                            Pedigree))

write_rds(yield, "data/phenotype/yield_munged.rds")
