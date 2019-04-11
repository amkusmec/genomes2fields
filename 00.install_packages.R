### Run once:

install.packages(c("devtools", "tidyverse", "lubridate", "lutz", "prism", 
                   "purrrlyr", "jsonlite", "readxl", "mgcv", "sommer", 
                   "evolqg", "tensorBSS", "WGCNA", "randomForest", 
                   "glmnet", "SGL"))

### Ubuntu 14.04 requires the installation of an old version of `rgdal`
devtools::install_version("rgdal", version = "1.2-20", 
                          repos = "http://cran.us.r-project.org")
install.packages("elevatr")

### Work-around due to using R 3.4.4; if using a newer version of R, you
### may be able to use the normal installation procedure for Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
install.packages("isva")
