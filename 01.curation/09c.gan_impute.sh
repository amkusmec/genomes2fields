#!/bin/bash

cd data/gbs
for i in `seq 1 10`; do
  java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom${i}.vcf.gz gp=true out=GAN_reference_chrom${i} nthreads=8
  bgzip -d GAN_reference_chrom${i}.vcf.gz
  python3 ~/snptools/src/snpstat.py -i GAN_reference_chrom${i}.vcf -o GAN_reference_chrom${i}.stat -mi 4
done
