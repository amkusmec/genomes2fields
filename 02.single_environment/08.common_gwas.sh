#!/bin/bash

cd data/single_site_gwas

# Single-environment GWAS with environment-specific kinship matrices
for i in `seq 1 46`; do
  ../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p yield_blue.txt -n ${i} -a snp_info.txt -snps common_snps.txt -gk 2 -maf 0 -o common_g2f_related${i}
  ../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p yield_blue.txt -n ${i} -a snp_info.txt -snps common_snps.txt -c covariates.txt -k output/common_g2f_related${i}.sXX.txt -lmm 4 -maf 0 -o common_single${i}
done
