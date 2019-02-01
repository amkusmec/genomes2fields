library(tidyverse)


# Download data for imputation --------------------------------------------
### See 01.munge_metadata.R for information about configuring icommands.

### This line will need to change based on where you installed icommands.
Sys.setenv("PATH" = paste("/mnt/01/amkusmec/bin/icommands", Sys.getenv("PATH"), sep = ":"))

# Get HapMap1 data from CyVerse
setwd("data/gbs/HapMap1")
system("iinit")
system(paste0("iget /iplant/home/shared/panzea/genotypes/HapMap/v1/", 
              "maizeHapMapV1_B73RefGenV2_20110309.hmp.txt.gz"))
system("gunzip maizeHapMapV1_B73RefGenV2_20110309.hmp.txt.gz")

# Get HapMap2 data from CyVerse
setwd("../HapMap2")
for (i in 1:10) {
  system(paste0("iget /iplant/home/shared/panzea/genotypes/HapMap/v2/", 
                "maizeHapMapV2_B73RefGenV2_201203028_chr", i, ".hmp.txt.gz"))
  system(paste0("gunzip maizeHapMapV2_B73RefGenV2_201203028_chr", i, ".hmp.txt.gz"))
}

# Get NAM marker data from CyVerse
setwd("../NAM")
system("iget /iplant/home/shared/panzea/genotypes/SNPs/NAM_map_and_genos-120731.zip")
system("unzip NAM_map_and_genos-120731.zip")
system("iexit")


# Get the RNA-seq SNP matrix ----------------------------------------------
setwd("..")
system(paste0("cp /media/analyses/NSF-CNV.NAM_Ear+Root+Shoot+Tassel+Apex.", 
              "AGPv2-20181003/SNPs/NAM_Ear+Root+Shoot+Tassel+Apex.AGPv2-", 
              "20181003.SNP_Matrix.Recalled.genotypes.vcf.gz NAM_rna.vcf.gz"))
system("gunzip NAM_rna.vcf.gz")
