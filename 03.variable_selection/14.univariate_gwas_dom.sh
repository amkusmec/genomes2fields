#!/bin/bash

cd data/gemma/dominance

../../../src/gemma-0.98.1-linux-static -g g2f_hybrids_dom.txt -p ../rxn_norm.txt -n 1 -a ../snp_info.txt -c covariates_dom.txt -k dom_mat.txt -lmm 4 -maf 0 -o norm_Hybrid_dom
../../../src/gemma-0.98.1-linux-static -g g2f_hybrids_dom.txt -p ../rxn_norm.txt -n 2 -a ../snp_info.txt -c covariates_dom.txt -k dom_mat.txt -lmm 4 -maf 0 -o norm_lnMSE_dom
../../../src/gemma-0.98.1-linux-static -g g2f_hybrids_dom.txt -p ../rxn_norm.txt -n 3 -a ../snp_info.txt -c covariates_dom.txt -k dom_mat.txt -lmm 4 -maf 0 -o norm_NET_X0.025_X0.45_dom
../../../src/gemma-0.98.1-linux-static -g g2f_hybrids_dom.txt -p ../rxn_norm.txt -n 4 -a ../snp_info.txt -c covariates_dom.txt -k dom_mat.txt -lmm 4 -maf 0 -o norm_NET_X0.65_X1.275_dom
../../../src/gemma-0.98.1-linux-static -g g2f_hybrids_dom.txt -p ../rxn_norm.txt -n 5 -a ../snp_info.txt -c covariates_dom.txt -k dom_mat.txt -lmm 4 -maf 0 -o norm_TMAX_X0.775_X0.875_dom
../../../src/gemma-0.98.1-linux-static -g g2f_hybrids_dom.txt -p ../rxn_norm.txt -n 6 -a ../snp_info.txt -c covariates_dom.txt -k dom_mat.txt -lmm 4 -maf 0 -o norm_TMIN_X0.05_X1.5_dom
../../../src/gemma-0.98.1-linux-static -g g2f_hybrids_dom.txt -p ../rxn_norm.txt -n 7 -a ../snp_info.txt -c covariates_dom.txt -k dom_mat.txt -lmm 4 -maf 0 -o norm_TMIN_X0.15_X0.65_dom
