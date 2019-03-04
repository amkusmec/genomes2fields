#!/bin/bash

cd data/gbs
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom1_filt.vcf.gz out=GAN_reference_chrom1 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom2_filt.vcf.gz out=GAN_reference_chrom2 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom3_filt.vcf.gz out=GAN_reference_chrom3 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom4_filt.vcf.gz out=GAN_reference_chrom4 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom5_filt.vcf.gz out=GAN_reference_chrom5 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom6_filt.vcf.gz out=GAN_reference_chrom6 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom7_filt.vcf.gz out=GAN_reference_chrom7 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom8_filt.vcf.gz out=GAN_reference_chrom8 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom9_filt.vcf.gz out=GAN_reference_chrom9 nthreads=8
java -Xmx100g -jar ../../src/beagle.28Sep18.793.jar gt=GAN_combined_chrom10_filt.vcf.gz out=GAN_reference_chrom10 nthreads=8
