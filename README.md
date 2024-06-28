## Data-driven identification of environmental variables influencing phenotypic plasticity to facilitate breeding for future climates

### Introduction

This project uses hybrid maize grain yield trials from the [Genomes to Fields Initiative](https://www.genomes2fields.org/home/) in 2014-2016 to develop and apply a genetic algorithm for the identification of environmental covariates driving grain yield phenotypic plasticity. A paper describing the results is forthcoming in **New Phytologist**.

### Description

Prior to running any of the scripts in the repository, please run `00.install_packages.R` to install all R packages required to reproduce these analyses. The necessary directory structure has been preserved through the use of `.gitignore` files.

All data necessary to reproduce the results are publicly available from the sources in the "Data Availabilitly Statement" of the paper linked above.

01. **curation**: These scripts curate and assemble the yield, weather, and genotypic data.
02. **single_environment**: These scripts perform single-environment phenotypic analyses of the yield data. They produce best linear unbiased estimates of grain yield and associated variances suitable for a two-stage analysis.
03. **variable_selection**: These scripts run a genetic algorithm to identify a subset of environmental variables that are most informative of the grain yield plasticity patterns across 45 trials. Reaction norm parameters for 727 maize hybrids are then estimated using Bayesian random regressions.
04. **genetic_analysis**: These scripts perform GWAS on the reaction norm parameters, analyze pleiotropy, identify candidate genes, and perform Gene Ontology enrichment analysis.
05. **figures**: These scripts construct the main text and more complicated supplementary figures and tables.

### License

This repository is free and open source for use in research and is licensed under the terms of [GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/#).
