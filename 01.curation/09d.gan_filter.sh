#!/bin/bash

cd data/gbs
for i in `seq 1 10`; do
  ~/htslib/bgzip -d GAN_reference_chrom${i}.vcf.gz
  python3 ~/snptools/src/snpstat.py -i GAN_reference_chrom${i}.vcf \ 
    -o GAN_reference_chrom${i}.stat -mi 4
  python3 ~/snptools/src/filter.py -s GAN_reference_chrom${i}.stat \
    -i GAN_reference_chrom${i}.vcf -o GAN_reference_chrom${i}_filt -mi 4 \
    -ht 0.05
done
