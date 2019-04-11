#!/bin/bash

cd data/gbs

# Impute two families of NAM RILs with the phased and imputed HapMap1, HapMap2,
# and RNA-seq SNPs for the NAM parents
for i in `seq 1 10`; do
  java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom${i}.vcf ref=NAM_reference_chrom${i}.vcf.gz out=NAM_imputed_chrom${i} nthreads=8
done
