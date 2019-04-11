#!/bin/bash

for i in `seq 1 10`; do
  # Filter NAM RIL SNPs
  bgzip NAM_imputed_chrom${i}.vcf
  tabix -p vcf NAM_imputed_chrom${i}.vcf.gz
  bcftools view -S rils_to_keep.txt --force-samples -O z -o NAM_imputed_chrom${i}_sub.vcf.gz NAM_imputed_chrom${i}.vcf.gz
  tabix -p vcf NAM_imputed_chrom${i}_sub.vcf.gz
  bcftools view -O z -o NAM_imputed_chrom${i}_filt.vcf.gz --include 'F_MISSING<=0.8 & MAF>=0.025 & F_PASS(GT="het")<=0.05' NAM_imputed_chrom${i}_sub.vcf.gz
  tabix -p vcf NAM_imputed_chrom${i}_filt.vcf.gz
  
  # Filter and impute combined G2F/Ames SNPs
  bgzip G2F_Ames_combined_chrom${i}.vcf
  tabix -p vcf G2F_Ames_combined_chrom${i}.vcf.gz
  bcftools view -O z -o G2F_Ames_combined_chrom${i}_filt.vcf.gz --include 'F_MISSING<=0.8 & MAF>=0.025 & F_PASS(GT="het")<=0.05' G2F_Ames_combined_chrom${i}.vcf.gz
  tabix -p vcf G2F_Ames_combined_chrom${i}_filt.vcf.gz
  java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=G2F_Ames_combined_chrom${i}_filt.vcf.gz out=G2F_Ames_combined_chrom${i}_imputed nthreads=8
  
  # Filter imputed G2F/Ames SNPs
  bgzip -d G2F_Ames_combined_chrom${i}_imputed.vcf.gz
  bgzip G2F_Ames_combined_chrom${i}_imputed.vcf
  tabix -p vcf G2F_Ames_combined_chrom${i}_imputed.vcf.gz
  bcftools view -O z -o G2F_Ames_combined_chrom${i}_imp_filt.vcf.gz --include 'F_MISSING<=0.8 & MAF>=0.025 & F_PASS(GT="het")<=0.05' G2F_Ames_combined_chrom${i}_imputed.vcf.gz
  tabix -p vcf G2F_Ames_combined_chrom${i}_imp_filt.vcf.gz
  
  # Combine both SNP sets
  bcftools merge -m snps -O z -o GAN_combined_chrom${i}.vcf.gz NAM_imputed_chrom${i}_filt.vcf.gz G2F_Ames_combined_chrom${i}_imp_filt.vcf.gz
done
