#!/bin/bash

###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 07a RADseq filtering - Stylophora pistillata #####

# Infer population structure using the model based approach implemented in ADMIXTURE

cd ~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_ADMIXTURE_20210222
#create PED file from LE filtered vcf file
~/RADseq/Genotyping/tools/plink --vcf ../spis_vcftools_filtered_vcf/spis.LE.filtered.recode.indnames.vcf --make-bed --allow-extra-chr --out spis.367ind.LE.plink2recode
~/RADseq/Genotyping/tools/plink --bfile spis.367ind.LE.plink2recode --recode 12 --allow-extra-chr --out spis.367ind.LE.plink

# Run ADMIXTURE using cv=10. This was done in SYMBIOMICS /home/buitracn/RADseq/RAD_Big-dataset/spis-bam-files/stacks2.0/gstacks.448samples/test.Thierry/p1.mac4.r0.8.448sample/NEW/03_vcftools_filtered_vcf/For_ADMIXTURE
export PATH="/home/buitracn/Software/anaconda2/bin:$PATH"

cd run_cve_10fold_1
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_2
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_3
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_4
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_5
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_6
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

cd  ../run_cve_10fold_7
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

cd  ../run_cve_10fold_8
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

cd  ../run_cve_10fold_9
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

cd  ../run_cve_10fold_10
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 spis.367ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4
