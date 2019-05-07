### Run once:

install.packages(c("devtools", "tidyverse", "lubridate", "lutz", "prism", 
                   "purrrlyr", "jsonlite", "readxl", "mgcv", "sommer", 
                   "evolqg", "tensorBSS", "WGCNA", "randomForest", 
                   "glmnet", "SGL", "spikeslab"))

### Ubuntu 14.04 requires the installation of an old version of `rgdal`
devtools::install_version("rgdal", version = "1.2-20", 
                          repos = "http://cran.us.r-project.org")
install.packages("elevatr")

### Work-around due to using R 3.4.4; if using a newer version of R, you
### may be able to use the normal installation procedure for Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
install.packages("isva")

devtools::install_github("stephenslab/susieR@0.8.0")


### Need to figure out how to install `ashr` directly and not as part of `mashr`
### unless we end up using the latter after figuring out how to deal with the
### missing data issue.

