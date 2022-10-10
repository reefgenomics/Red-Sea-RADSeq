#!/bin/bash

###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 08a Linkage Disequilibrium decay - Stylophora pistillata #####

cd ~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_LD_20210303

# Create symbolic links to the vcf files to recode that contains only the whitelist of markers that pass filter_rad filters (except LD) and was outputed by STACKS with the --ordered-output flag
ln -s /home/buitracn/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_filter_rad_20210221@1444/genetic_diversity_in_STACKS/For_Pi_STACKS/ByREEF/populations.snps.vcf .
sed 's/_filtered-sorted//g' populations.snps.vcf > populations.snps.indnames.vcf

## Recode vcf file to contain only the individuals of each genetic..SCLuster (K=6) and a minor allele frequcy greater than 5%
cd K6
for i in SCL1 SCL2 SCL3 SCL4 SCL5 SCL6; do vcftools --vcf ../populations.snps.indnames.vcf  --keep "$i".txt --maf 0.05 --out spis.majorcLust."$i" --recode --recode-INFO-all;done 

# estimate LD (r2) for individual genetic..SCLuster
# assess ld from vcf file directly
for i in SCL1 SCL2 SCL3 SCL4 SCL5 SCL6 ; do ~/RADseq-Big-project/tools/plink --vcf spis.majorclust."$i".recode.vcf --make-bed --allow-extra-chr --out spis.majorclust."$i".sorted.plink; done # to sort markers in chromosomes
for i in SCL1 SCL2 SCL3 SCL4 SCL5 SCL6 ; do ~/RADseq-Big-project/tools/plink --bfile spis.majorclust."$i".sorted.plink --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 999999 --ld-window-kb 2095 --out spis."$i".LD.r2.plink; done

# Select only SNPs in the 10 largest Scaffolds
for i in SCL1 SCL2 SCL3 SCL4 SCL5 SCL6 ; do grep -wFf ../Spis.10largest.Scaffolds spis."$i".LD.r2.plink.ld > spis."$i".LD.r2.plink.10largestScaff.ld; done

# To count how many SNPs per population MAF >5% in the 10 largest Scaffolds
for i in SCL1 SCL2 SCL3 SCL4 SCL5 SCL6  ; do grep -wFf ../Spis.10largest.Scaffolds spis.majorclust."$i".recode.vcf | grep -v "#" | wc -l ; done
# 6637
# 7062
# 6980
# 7582
# 6054
# 7889

