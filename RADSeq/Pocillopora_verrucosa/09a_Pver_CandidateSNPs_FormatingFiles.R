###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 09a Candidate SNPs for Positive Selection  - Pocillopora verrucosa #####

#### 09a.01. Formating SNP file for BAYESCAN

# libraries
library(hierfstat)
library(vcfR)
library(adegenet)
library(poppr)

# set working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_positive_selection_analyses_20210308")

# reformat the vcffile to contain only the individuals that had over 90% assignment to a single genetic cluster
# ln -s ~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_individuals_split_by_cluster_20210302/K2/pver.genclust.strata.K2.tsv .
# awk '{print $1}' pver.genclust.strata.K2.tsv | sed 1,1d > pver.genclust.ind.noadmix #203 individuals
# vcftools --vcf pver.298.LE.filtered.recode.indnames.vcf --keep pver.genclust.ind.noadmix --recode --recode-INFO-all --out pver.noadmix.clust


# generate a strata file
pver.strata <- read.delim("pver.genclust.strata.K2.tsv", h=T, sep = "")
pver.strata$STRATA <- factor(pver.strata$STRATA, levels = c("PCL1", "PCL2"), ordered = T)

# read vcf file and generate a genind object
pver.vcfR <- read.vcfR(file = "pver.noadmix.clust.recode.vcf")
pver.genind <- vcfR2genind(pver.vcfR,  ind.names= pver.strata$INDIVIDUALS, pop = pver.strata$STRATA)

# convert genind to hierfstat and write a BAYESCAN format
pver.hierfstat <- genind2hierfstat(pver.genind)
write.bayescan(pver.hierfstat, fn = "./BAYESCAN/pver.popgenclust.bayescan") #write a bayescan file

# create an SNPs id dictionary BAYESCAN
system("grep -v "#" pver.noadmix.clust.recode.vcf | awk 'BEGIN{OFS="\t";}{sub(/\_size.*$/,"",$1); print $1,$2,$3}' | sed 's/Pver_Sc//' | awk -v OFS='\t' '{print$0, NR}'  > snps.id.dictionary.BAYESCAN.txt") #CHROM, POS, ID, index
       
# There are some SNPs that turn to be monomorphic in the current population configuration. We identify them and create a list to remove them when running tha bayescan analysis
system("awk '$3 == "1" {print $1}' ./BAYESCAN/pver.popgenclust.bayescan > ./BAYESCAN/monomophicsnpsindex2remove.txt") # 8 monomorphic markers in all population (genetic clusters)
       
###
# convert the BAYESCAN format to BAYPASS format (remove 15 monomorphic SNPs (90 rows))
system("sed 1,5d ./BAYESCAN/pver.popgenclust.bayescan | awk -F" " '$3 != "1" {print $4, $5}' | sed 1,1d | sed '/^ $/d'> tmp") # 1) remove header information; 2) print only columns with the allele counts per populations excluding monomorphic SNPs; 3) remove empty lines
       
# create the genotype matrix for BAYPASS (35,207-4 = 35203 polymorphic markers)
system("awk -v numRows=35203 -f tst.awk tmp  > ./BAYPASS/pver.popgenclust.baypass") # split the temporal file into 25318 that corresponds to the number of snps in the data
       
# create an SNPs id dictionary BAYPASS
system("awk 'NR==FNR{a[$1]++;next} !($4 in a)' monomophicsnpsindex2remove.txt snps.id.dictionary.BAYESCAN.txt | awk 'BEGIN{OFS="\t";}{print $1, $2, $3}' | awk -v OFS='\t' '{print$0, NR}' > snps.id.dictionary.BAYPASS.txt") #CHROM, POS, ID, index
       
## REMEMBER THAT BAYESCAN AND PAYPASS DICTIONARY ARE DIFFERENT BECAUSE I HAD TO MANUALLY REMOVE MONOMORPHIC SNPS FROM THE BAYPASS GENOTYPE FILE (15 SNPs)
       
