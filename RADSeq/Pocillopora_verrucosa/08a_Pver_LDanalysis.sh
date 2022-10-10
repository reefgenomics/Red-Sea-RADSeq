#!/bin/bash

###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 08a Linkage Disequilibrium decay - Pocillopora verrucosa #####
cd ~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_LD_20210307

# Create symbolic links to the vcf files to recode that contains only the whitelist of markers that pass filter_rad filters (except LD) and was outputed by STACKS with the --ordered-output flag
ln -s /home/buitracn/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_filter_rad_20210222@1117/genetic_diversity_in_STACKS/For_Pi_STACKS/ByREEF/populations.snps.vcf .
sed 's/_filtered-sorted//g' populations.snps.vcf > populations.snps.indnames.vcf

## Recode vcf file to contain only the individuals of each genetic cluster (K=7)
for i in PCL1 PCL2; do vcftools --vcf populations.snps.indnames.vcf  --keep "$i".txt --maf 0.05 --out pver.majorclust."$i" --recode --recode-INFO-all;done 

# estimate LD (r2) for individual genetic cluster
# assess ld from vcf file directly
for i in PCL1 PCL2 ; do ~/RADseq-Big-project/tools/plink --vcf pver.majorclust."$i".recode.vcf --make-bed --allow-extra-chr --out pver.majorclust."$i".sorted.plink; done # to sort markers in chromosomes
for i in PCL1 PCL2 ; do ~/RADseq-Big-project/tools/plink --bfile pver.majorclust."$i".sorted.plink --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 999999 --ld-window-kb 2095 --out pver."$i".LD.r2.plink; done

# Select only SNPs in the 10 largest Scaffolds
for i in PCL1 PCL2 ; do grep -wFf Pver.10largest.Scaffolds pver."$i".LD.r2.plink.ld > pver."$i".LD.r2.plink.10largestScaff.ld; done

# To count how many SNPs per population MAF >5% in the 10 largest Scaffolds
for i in PCL1 PCL2 ; do grep -wFf Pver.10largest.Scaffolds pver.majorclust."$i".recode.vcf | grep -v "#" | wc -l ; done
# 13730
# 14037

