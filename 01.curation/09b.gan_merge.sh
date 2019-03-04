#!/bin/bash

for i in `seq 1 10`; do
  ~/htslib/bgzip NAM_imputed_chrom${i}.vcf
  ~/htslib/tabix -p vcf NAM_imputed_chrom${i}.vcf.gz
  bcftools view -S rils_to_keep.txt --force-samples -o NAM_imputed_chrom${i}_sub.vcf.gz \
    -O z NAM_imputed_chrom${i}.vcf.gz
  ~/htslib/tabix -p vcf NAM_imputed_chrom${i}_sub.vcf.gz
  
  ~/htslib/bgzip G2F_Ames_combined_chrom${i}.vcf
  ~/htslib/tabix -p vcf G2F_Ames_combined${i}.vcf.gz
  bcftools merge -m snps -O z -o GAN_combined_chrom${i}.vcf.gz \
    NAM_imputed_chrom${i}_sub.vcf.gz G2F_Ames_combined_chrom${i}.vcf.gz
  
  ~/htslib/bgzip -d GAN_combined_chrom${i}.vcf.gz
  python3 ~/snptools/src/snpstat.py -i GAN_combined_chrom${i}.vcf \
    -o GAN_combined_chrom${i}.stat -mi 4
  python3 ~/snptools/src/filter.py -s GAN_combined_chrom${i}.stat \
    -i GAN_combined_chrom${i}.vcf -o GAN_combined_chrom${i}_filt -mi 4 \
    -n 0.9 -f 0.05 -ht 0.05
  ~/htslib/bgzip GAN_combined_chrom${i}_filt.vcf
done

rm *filtered*
