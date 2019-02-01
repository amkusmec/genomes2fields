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
~/tassel-5-standalone/run_pipeline.pl -Xmx10g -h NAM_markers.hmp.txt \
  -export NAM_markers.vcf -exportType VCF

# Combine all the parent VCF files
cd ..
~/tassel-5-standalone/run_pipeline.pl -Xmx200g -fork1 -vcf HapMap1/maizeHapMapV1_recode.vcf \
  -fork2 -vcf HapMap2/maizeHapMapV2_recode.vcf -fork3 -vcf NAM_rna_recode.vcf \
  -combine4 -input1 -input2 -input3 -mergeGenotypeTables -export NAM_merged.vcf \
  -exportType VCF

# Create reference genotypes for the NAM parents
java -Xmx200g -jar ../../src/beagle.28Sep18.793.jar gt=NAM_merged.vcf \
  out=NAM_reference map=NAM/NAM_markers.map nthreads=24

# Impute the NAM RILs
java -Xmx300g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers.vcf \
  ref=NAM_reference.vcf.gz out=NAM_imputed map=NAM/NAM_markers.map nthreads=24
