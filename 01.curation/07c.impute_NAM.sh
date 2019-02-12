#!/bin/bash

# Convert HapMap to VCF and recode indels
cd data/gbs/HapMap1
~/tassel-5-standalone/run_pipeline.pl -Xmx50g -h maizeHapMapV1.hmp.txt \
  -export maizeHapMapV1.vcf -exportType VCF
python3 ../../../src/recode_indels.py -i maizeHapMapV1.vcf -o maizeHapMapV1_recode

cd ../HapMap2
~/tassel-5-standalone/run_pipeline.pl -Xmx50g -h maizeHapMapV2.hmp.txt \
  -export maizeHapMapV2.vcf -exportType VCF
python3 ../../../src/recode_indels.py -i maizeHapMapV2.vcf -o maizeHapMapV2_recode

cd ..
python3 ../../src/recode_indels.py -i NAM_rna.vcf -o NAM_rna_recode

cd NAM
# ~/tassel-5-standalone/run_pipeline.pl -Xmx10g -h NAM_markers.hmp.txt \
#   -export NAM_markers -exportType VCF
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

# Create reference genotypes for the NAM parents
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom1.vcf out=NAM_reference_chrom1 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom2.vcf out=NAM_reference_chrom2 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom3.vcf out=NAM_reference_chrom3 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom4.vcf out=NAM_reference_chrom4 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom5.vcf out=NAM_reference_chrom5 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom6.vcf out=NAM_reference_chrom6 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom7.vcf out=NAM_reference_chrom7 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom8.vcf out=NAM_reference_chrom8 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom9.vcf out=NAM_reference_chrom9 nthreads=4
java -Xmx50g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged_chrom10.vcf out=NAM_reference_chrom10 nthreads=4

# Impute the NAM RILs
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom1.vcf ref=NAM_reference_chrom1.vcf.gz window=45 out=NAM_imputed_chrom1 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom2.vcf ref=NAM_reference_chrom2.vcf.gz window=45 out=NAM_imputed_chrom2 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom3.vcf ref=NAM_reference_chrom3.vcf.gz window=45 out=NAM_imputed_chrom3 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom4.vcf ref=NAM_reference_chrom4.vcf.gz window=45 out=NAM_imputed_chrom4 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom5.vcf ref=NAM_reference_chrom5.vcf.gz window=45 out=NAM_imputed_chrom5 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom6.vcf ref=NAM_reference_chrom6.vcf.gz window=45 out=NAM_imputed_chrom6 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom7.vcf ref=NAM_reference_chrom7.vcf.gz window=45 out=NAM_imputed_chrom7 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom8.vcf ref=NAM_reference_chrom8.vcf.gz window=45 out=NAM_imputed_chrom8 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom9.vcf ref=NAM_reference_chrom9.vcf.gz window=45 out=NAM_imputed_chrom9 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom10.vcf ref=NAM_reference_chrom10.vcf.gz window=45 out=NAM_imputed_chrom10 nthreads=8
