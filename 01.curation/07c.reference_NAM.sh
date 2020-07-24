#!/bin/bash

### Convert the remaining SNP files, combine, and filter. Estimate imputation 
### accuracy through masking. Create reference NAM parent genotypes through 
### self-imputation.

# Convert HapMap to VCF and recode indels
cd data/gbs/HapMap1
~/tassel-5-standalone/run_pipeline.pl -Xmx50g -h maizeHapMapV1.hmp.txt \
  -export maizeHapMapV1.vcf -exportType VCF
python3 ../../../src/recode_indels.py -i maizeHapMapV1.vcf -o maizeHapMapV1_recode
rm maizeHapMapV1.vcf

cd ../HapMap2
~/tassel-5-standalone/run_pipeline.pl -Xmx50g -h maizeHapMapV2.hmp.txt \
  -export maizeHapMapV2.vcf -exportType VCF
python3 ../../../src/recode_indels.py -i maizeHapMapV2.vcf -o maizeHapMapV2_recode
rm maizeHapMapV2.hmp.txt maizeHapMapV2.vcf

cd ..
python3 ../../src/recode_indels.py -i NAM_rna.vcf -o NAM_rna_recode

cd NAM
~/tassel-5-standalone/run_pipeline.pl -Xmx10g -h NAM_markers.hmp.txt -separate \
  -export NAM_markers_chrom -exportType VCF

# Combine all the parent VCF files
cd ..
~/tassel-5-standalone/run_pipeline.pl -Xmx200g -fork1 -vcf HapMap1/maizeHapMapV1_recode.vcf \
  -fork2 -vcf HapMap2/maizeHapMapV2_recode.vcf -fork3 -vcf NAM_rna_recode.vcf \
  -combine4 -input1 -input2 -input3 -mergeGenotypeTables -export NAM_merged.vcf \
  -exportType VCF
~/tassel-5-standalone/run_pipeline.pl -Xmx100g -vcf NAM_merged.vcf -separate \
  -export NAM_merged_chrom -exportType VCF
rm NAM_merged.vcf

# QC on NAM parents
for i in `seq 1 10`; do
  python3 ~/snptools/src/snpstat.py -i NAM_merged_chrom${i}.vcf -o NAM_merged_chrom${i}.stat -mi 4
  python3 ~/snptools/src/filter.py -s NAM_merged_chrom${i}.stat -i NAM_merged_chrom${i}.vcf -o NAM_merged_chrom${i}_filt -mi 4 -n 0.8 -f 0.05 -ht 0.05
done
rm *_filt_filtered.vcf

# Estimate imputation accuracy
for i in `seq 1 10`; do
  for j in `seq 1 10`; do
    echo $i $j
    python3 ../../src/mask_vcf.py -i NAM_merged_chrom${i}_filt.vcf -o NAM_merged_chrom${i}_mask
    java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom${i}_mask.vcf out=NAM_merged_chrom${i}_acc nthreads=4
    Rscript ../../src/imputation_accuracy.R -s NAM_merged_chrom${i}.stat -r NAM_merged_chrom${i}_filt.vcf -x NAM_merged_chrom${i}_mask_index.txt -i NAM_merged_chrom${i}_acc.vcf.gz -o NAM_merged_chrom${i}_acc -c $j
  done
done

rm *mask*

# Create reference genotypes for the NAM parents
for i in `seq 1 10`; do
  java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom${i}_filt.vcf out=NAM_reference_chrom${i} nthreads=4
done
