###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 09a Candidate SNPs for Positive Selection  - Stylophora pistillata #####

#### 09a.01. Formating SNP file for BAYESCAN

# libraries
library(hierfstat)
library(vcfR)
library(adegenet)
library(poppr)

# set working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_positive_selection_analyses_20210308")

# reformat the vcffile to contain only the individuals that had over 90% assignment to a single genetic cluster
# ln -s ~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_individuals_split_by_cluster_20210302/K6/spis.genclust.strata.K6.tsv .
# awk '{print $1}' spis.genclust.strata.K6.tsv | sed 1,1d > spis.genclust.ind.noadmix #312 individuals
# vcftools --vcf spis.LE.filtered.recode.indnames.vcf --keep spis.genclust.ind.noadmix --recode --recode-INFO-all --out spis.noadmix.clust

# generate a strata file
spis.strata <- read.delim("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_individuals_split_by_cluster_20210302/K6/spis.genclust.strata.K6.tsv", h=T, sep = "")
spis.strata$STRATA <- factor(spis.strata$STRATA, levels = c("SCL1", "SCL2", "SCL3", "SCL4", "SCL5", "SCL6"), ordered = T)


# read vcf file and generate a genind object
spis.vcfR <- read.vcfR(file = "spis.noadmix.clust.recode.vcf")
spis.genind <- vcfR2genind(spis.vcfR,  ind.names= spis.strata$INDIVIDUALS, pop = spis.strata$STRATA)

# convert genind to hierfstat and write a BAYESCAN format
spis.hierfstat <- genind2hierfstat(spis.genind)
write.bayescan(spis.hierfstat, fn = "./BAYESCAN/spis.popgenclust.bayescan") #write a bayescan file

# create an SNPs id dictionary BAYESCAN
system("grep -v '#' spis.noadmix.clust.recode.vcf | awk 'BEGIN{OFS="\t";}{sub(/\|size.*$/,"",$1); print $1,$2,$3}' | sed 's/Spis.scaffold//' | awk -v OFS='\t' '{print$0, NR}'  > snps.id.dictionary.BAYESCAN.txt") #CHROM, POS, ID, index

# There are some SNPs that turn to be monomorphic in the current population configuration. We identify them and create a list to remove them when running tha bayescan analysis
system("awk '$3 == "1" {print $1}' ./BAYESCAN/spis.popgenclust.bayescan > ./BAYESCAN/monomophicsnpsindex2remove.txt") #15 monomorphic markers in all population (genetic clusters)


#### 09a.02. Formating SNP file for BAYPASS (25,318 SNPs, 312 individuals)

# convert the BAYESCAN format to BAYPASS format (remove 15 monomorphic SNPs (90 rows))
system("sed 1,5d ./BAYESCAN/spis.popgenclust.bayescan | awk -F" " '$3 != "1" {print $4, $5}' | sed 1,1d | sed '/^ $/d'> tmp") # 1) remove header information; 2) print only columns with the allele counts per populations excluding monomorphic SNPs; 3) remove empty lines

# create the genotype matrix for BAYPASS (25,318-15 = 25,303 polymorphic markers)
system("awk -v numRows=25303 -f tst.awk tmp  > ./BAYPASS/spis.popgenclust.baypass") 

# create an SNPs id dictionary BAYPASS
system("awk 'NR==FNR{a[$1]++;next} !($4 in a)' monomophicsnpsindex2remove.txt snps.id.dictionary.BAYESCAN.txt | awk 'BEGIN{OFS="\t";}{print $1, $2, $3}' | awk -v OFS='\t' '{print$0, NR}' > snps.id.dictionary.BAYPASS.txt") #CHROM, POS, ID, index

## REMEMBER THAT BAYESCAN AND PAYPASS DICTIONARY ARE DIFFERENT BECAUSE I HAD TO MANUALLY REMOVE MONOMORPHIC SNPS FROM THE BAYPASS GENOTYPE FILE (15 SNPs)

