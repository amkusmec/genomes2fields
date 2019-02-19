#!/bin/bash

cd data/gbs

# Impute two families of NAM RILs with the phased and imputed HapMap1, HapMap2,
# and RNA-seq SNPs for the NAM parents
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom1.vcf ref=NAM_reference_chrom1.vcf.gz out=NAM_imputed_chrom1 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom2.vcf ref=NAM_reference_chrom2.vcf.gz out=NAM_imputed_chrom2 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom3.vcf ref=NAM_reference_chrom3.vcf.gz out=NAM_imputed_chrom3 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom4.vcf ref=NAM_reference_chrom4.vcf.gz out=NAM_imputed_chrom4 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom5.vcf ref=NAM_reference_chrom5.vcf.gz out=NAM_imputed_chrom5 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom6.vcf ref=NAM_reference_chrom6.vcf.gz out=NAM_imputed_chrom6 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom7.vcf ref=NAM_reference_chrom7.vcf.gz out=NAM_imputed_chrom7 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom8.vcf ref=NAM_reference_chrom8.vcf.gz out=NAM_imputed_chrom8 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom9.vcf ref=NAM_reference_chrom9.vcf.gz out=NAM_imputed_chrom9 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=NAM/NAM_markers_chrom10.vcf ref=NAM_reference_chrom10.vcf.gz out=NAM_imputed_chrom10 nthreads=8
