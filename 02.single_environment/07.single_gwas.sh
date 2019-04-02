#!/bin/bash

cd data/gemma

# Single-environment GWAS with environment-specific kinship matrices
for i in `seq 1 46`; do
  ../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p yield_blue.txt -n ${i} -gk 2 -o g2f_related${i}
  ../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p yield_blue.txt -n ${i} -a snp_info.txt -c covariates.txt -k output/g2f_related${i}.sXX.txt -lmm 4 -maf 0.025 -o single${i}
done
