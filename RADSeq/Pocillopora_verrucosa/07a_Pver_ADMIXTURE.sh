#!/bin/bash

###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 07a RADseq filtering - Pocillopora verrucosa #####

# Infer population structure using the model based approach implemented in ADMIXTURE
#!/bin/bash

# Infer population structure using the model based approach implemented in ADMIXTURE

cd ~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_ADMIXTURE_20210223
#create PED file from LE filtered vcf file
~/RADseq/Genotyping/tools/plink --vcf ../03_vcftools_filtered_vcf/pver.LE.filtered.recode.vcf --make-bed --allow-extra-chr --out pver.296ind.LE.plink2recode #Ped file couldn't be created directly from the vcf file because Chromosomes were not sorted
~/RADseq/Genotyping/tools/plink --bfile pver.296ind.LE.plink2recode --recode 12 --allow-extra-chr --out pver.296ind.LE.plink

# Run ADMIXTURE using cv=10
# for i in 1 2  3 4 5 6 7 8 9 10; do mkdir run_cve_10fold_"$i"; done
export PATH="/home/buitracn/Software/anaconda2/bin:$PATH"

cd run_cve_10fold_1
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_2
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_3
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_4
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_5
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4


cd  ../run_cve_10fold_6
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

cd  ../run_cve_10fold_7
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

cd  ../run_cve_10fold_8
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

cd  ../run_cve_10fold_9
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

cd  ../run_cve_10fold_10
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 ; \
do admixture -s $(date +%s%N | cut -b10-19) --cv=10 pver.296ind.LE.plink.ped $K -j40| tee log${K}.out ; done
#grep -h CV log*.out | sort -k4

