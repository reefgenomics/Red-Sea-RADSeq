#!/bin/bash

###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 01 RADseq filtering - Stylophora pistillata #####

####  01.01 RADloci identification and samples genotyping using stacks v2.0 (bash environment) ####
cd ~/RADseq/Samples/Alignments/al-conc

# gstacks program will read in aligned reads to reference genome. Gstacks will assemble and merge paired-end contigs, call variant sites in the population and genotypes in each sample)
stacks-2.0/gstacks -I . -M /home/RADseq/Genotyping/spis-popmap-448samples --max-insert-len 700 -O /home/RADseq/Genotyping/Spis/ -t 24

####  01.02 Initial filtering to generate a VCF file ####
cd ~/RADseq/Genotyping/Spis
stacks-2.0/populations -P ./ -M ../spis-popmap-448samples -O ./p1.mac4.r0.8.448samples -p 1 -r 0.8 --min_mac 4 -e pstI --merge_sites --verbose --vcf --vcf_haplotypes -t 40

####  01.03 VCF filtering using the R package "radiator" (R environment)####
# The following lines were executed in R
# if (!require("devtools")) install.packages("devtools", , dependencies = TRUE)
# library(devtools)
# devtools::install_github("thierrygosselin/radiator", dependencies = TRUE)
library(radiator) #version 1.1.9

setwd("~/RADseq/Genotyping/Spis/p1.mac4.r0.8.448samples")

strata.spis <- read.delim("../../spis-reefmap-448samples.radiator.txt", h=T) # read in a file where each samples is associated to the sampling site
radiator::summary_strata(strata = strata.spis)

spis.data.filtered <- radiator::filter_rad(data = "populations.snps.vcf", strata = "../../spis-reefmap-448samples.radiator.txt", filename = "spis.filtered", pop.levels = c("MAQ_R1", "MAQ_R2", "WAJ_R1", "WAJ_R3", "WAJ_R4", "YAN_R1", "YAN_R3", "YAN_R4", "KAU_R1", "KAU_R2", "KAU_R3", "DOG_R1", "DOG_R2", "DOG_R3", "FAR_R1", "FAR_R2", "FAR_R3", "FAR_R4"))

# radiator::filter_rad is an interactive filtering tool were each filter can be chosen based on data visualization.
# The following filter thresholds were applied:
# Filter individuals with a maximum genotype missingness of 20% (0.2) 
# Filter SNPs based on their MAC of 10  (MAC range: [1 - 367]; MAF range: [0.0014 - 0.4973]) 
# Filter SNPs based on their min mean DP 15 and max mean DP 75 (~5X min meanDP) 
# Filtering markers based on maximum missing proportion, maximum missing proportion allowed 0.1 
# Filtering markers based on position in the read OFF
# Filtering markers based on number of SNPs per locus OFF (less than 20 per locus is ok)
# Filter individuals based on heterozygosity/missingness (0.052 / 0.1025). 
# Filter individuals based on potential clones (0.25)
# Filtering markers based on HWD. If a marker if in HWD in 8 reef or more with a mid - pvalue of 0.0001 then the SNP is blacklisted 

# After filtering
# 367 / 18 / 1735 / 58658 / 287194	==> individuals / strata / chrom / locus / markers	
# Data summary: 
# number of samples: 367
# number of markers: 287194
## The filtering data was outputted to the folder spis_filter_rad_20210221@1444 

# Check missing proportion of each individual after filtering was done
spis.filtered.LD <- read_rad(data = "./spis_filter_rad_20210221@1444/13_filtered/radiator_data_20210221@2011.rad")
spis.filtered.LD.gds <- tidy2gds(spis.filtered.LD) #radiator_20210224@1042.gds.rad
spis.ind.missing.prop <- detect_mixed_genomes(spis.filtered.LD.gds)

# Get diversity Stats of heterozygosity per Reef
spis.LD.data.sum.stats.REEF <- summary_rad(spis.data.filtered$output$tidy.data)

# Get diversity Stats of heterozygosity per Region
colnames(spis.data.filtered$output$tidy.data)[1] <- "POP_ID_REEF"
spis.data.filtered$output$tidy.data$POP_ID <- gsub(pattern = "_.*","", spis.data.filtered$output$tidy.data$POP_ID_REEF)
spis.LD.data.sum.stats.REGION <- summary_rad(spis.data.filtered$output$tidy.data)

# Filtering variants in LD
spis.LE.filtered <- filter_ld(spis.data.filtered$gds) #tidy data didn't work
# Filtering marker in LD at short distance using MAC statistics 0 / 0 / 228536 (/ chrom / locus / SNP) blacklisted
# Filtering marker in LD at long distance using MAC statistics / 0 / 33340 / 33340 (/ chrom / locus / SNP) blacklisted

# After filtering LD
# 367 / 18 / 1735 / 25318 / 25318 (individuals / strata / chrom / locus / SNP)
## The filtering data was outputted to the folder spis_filter_ld_20210221@2117

# whitelist of markers (SNPs) in Linkage Equilibrium within each scaffold
whitelist.markers <- read.delim("./spis_filter_ld_20210221@2117/whitelist.long.ld_0.1.tsv", sep = "")
colnames(whitelist.markers)
whitelist.markers <- whitelist.markers[,c(7,2,3,6)] #MARKER, "CHROM"   "LOCUS"   "POS"  
# create tidy dataset of filtered SNPs (Linkage equilibrium)
spis.LE.filtered.tidy <- radiator::filter_whitelist(data = spis.filtered.LD, whitelist.markers = whitelist.markers)

# Check heterozygosity and missing proportion of each individual after LD filtering was done
spis.LE.filtered.gds <- tidy2gds(spis.LE.filtered.tidy) #radiator_20210224@1056.gds.rad
spis.ind.missing.prop.LE <- detect_mixed_genomes(spis.LE.filtered.gds) 

####  01.04 Generate a filtered VCF file (bash environment) ####
# Generate list of individual to keep (367 individuals)
awk '{print $1}' spis_filter_rad_20210221@1444/13_filtered/strata.filtered.tsv | sed 1,1d | sed 's/-filtered/_filtered/' > spis_vcftools_filtered_vcf/spis.ind.to.keep # folder created to store the filter VCF file (in Linkage Equilibrium)

# Generate list of markers to keep (25,318 SNPs in LE)
awk 'BEGIN{OFS="\t";}{print $2, $6}' spis_filter_ld_20210221@2117/whitelist.long.ld_0.1.tsv | sed 1,1d | sed 's/Spis_/Spis./' > spis_vcftools_filtered_vcf/whitelist.snps.LE

cd ~/RADseq/Genotyping/Spis/p1.mac4.r0.8.448samples/spis_vcftools_filtered_vcf/ 
# Using VCFTOOLS filter the vcf file produced by STACKS for the filtered set of SNPs. This in this vcf file only one position in the chromosome is considered, that is if there were two SNPs called at the same position because they were cut by different restriction enzyme, STACKS will only output one SNP per site
# --ordered-export — if data is reference aligned, exports will be ordered; only a single representative of each overlapping site (stacks::populations)
ln -s ~/RADseq/Genotyping/Spis/p1.mac4.r0.8.448samples/populations.snps.vcf .

vcftools --vcf populations.snps.vcf --positions whitelist.snps.LE --keep spis.ind.to.keep --out spis.LE.filtered --recode --recode-INFO-all #After filtering, kept 25318 out of a possible 287192 Sites

# Order of individuals in vcf file (important for ADMIXTURE and PCA to define the order of individuals
grep "#CHROM" spis.LE.filtered.recode.indnames.vcf | sed 's/\t/\n/g' | sed 1,9d > spis.ind.order.in.vcf.txt #sed 1,9d deletes the first 9th lines of the vcf file

