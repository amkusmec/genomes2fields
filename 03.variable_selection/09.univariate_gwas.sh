#!/bin/bash

cd data/gemma

../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 8 -gk 2 -o g2f_related

../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 1 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_Hybrid
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 2 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_lnMSE
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 3 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_NET_X0.025_X0.45
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 4 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_NET_X0.65_X1.275
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 5 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_TMAX_X0.775_X0.875
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 6 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_TMIN_X0.05_X1.5
../../src/gemma-0.98.1-linux-static -g g2f_hybrids.txt -p rxn_norm.txt -n 7 -a snp_info.txt -c covariates.txt -k output/g2f_related.sXX.txt -lmm 4 -maf 0.025 -o norm_TMIN_X0.15_X0.65
