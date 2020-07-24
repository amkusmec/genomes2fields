#!/bin/bash

### Run GWAS with the MLM in GEMMA.

cd data/gemma

../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 6 -gk 2 -o g2f_related

../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 1 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_Hybrid
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 2 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_NET_X0.175_X1.4
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 3 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_SR_X0.65_X1.425
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 4 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_TMAX_X0.825_X0.95
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 5 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_TMAX_X0.875_X1.425
