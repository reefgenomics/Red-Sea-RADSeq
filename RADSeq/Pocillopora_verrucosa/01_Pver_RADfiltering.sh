#!/bin/bash

###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 01 RADseq filtering - Pocillopora verrucosa #####

####  01.01 RADloci identification and samples genotyping using stacks v2.0 ####
cd ~/RADseq/Samples/Alignments/al-conc

# gstacks program will read in aligned reads to reference genome. Gstacks will assemble and merge paired-end contigs, call variant sites in the population and genotypes in each sample)
stacks-2.0/gstacks -I . -M /home/RADseq/Genotyping/pver-popmap-316samples --max-insert-len 700 -O /home/RADseq/Genotyping/Pver/ -t 30

####  01.02 Initial filtering to generate a VCF file ####
cd ~/RADseq/Genotyping/Pver

stacks-2.0/populations -P ./ -M ../pver-popmap-316samples_nopop -O ./p1.mac4.r0.8.316samples/ -p 1 -r 0.8 --min_mac 4 -e pstI --merge_sites --lnl_lim -10 --verbose --vcf --vcf_haplotypes -t 30

####  01.03 VCF filtering using the R package "radiator" (R environment)####
# The following lines were executed in R
# if (!require("devtools")) install.packages("devtools", , dependencies = TRUE)
# library(devtools)
# devtools::install_github("thierrygosselin/radiator", dependencies = TRUE)
library(radiator) #version 1.1.9

setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples")

strata.pver <- read_strata("../../pver-reefmap-316samples.radiator.txt")$strata

radiator::summary_strata(strata = strata.pver)

pver.data.filtered <- radiator::filter_rad(data = "populations.snps.vcf", strata = strata.pver, filename = "pver.filtered", pop.levels = c("MAQ_R1", "MAQ_R2", "WAJ_R1", "WAJ_R2", "WAJ_R3", "YAN_R1", "YAN_R3", "KAU_R1", "KAU_R2", "KAU_R3", "DOG_R1", "DOG_R2", "DOG_R3", "FAR_R3"))

# radiator::filter_rad is an interactive filtering tool were each filter can be chosen based on data visualization.
# The following filter thresholds were applied:
# Filter individuals with a maximum genotype missingness of 20% (0.2) ==> 10 individuals blacklisted
# Filter SNPs based on their MAC of 8  (MAC range: [1 - 306]; MAF range: [0.0016 - 0.5]) ==>  / 286 / 5637 / 249275 (/ chrom / locus / SNP) blacklisted
# Filter SNPs based on their min mean DP 15 and max mean DP 75 (~5X min meanDP) ==>  520 / 18251 / 294222 (/ chrom / locus / SNP) blacklisted
# Filtering markers based on maximum missing proportion, maximum missing proportion allowed 0.1 ==> 391 / 5923 / 33832 (/ chrom / locus / SNP) blacklisted
# Filtering markers based on position in the read OFF
# Filtering markers based on number of SNPs per locus OFF (less than 20 per locus is ok)
# Filter individuals based on heterozygosity/missingness. ==> 6 individuals blacklisted due to outlier heterozygosity associated to missing data
# Filter individuals based on potential clones (there was only strong evidence for two pair of clones from same reef (2 clones in FAR and two clones in MAQ-R1)) PFAR-R3-15-filtered-sorted and PMAQ-R1-8-filtered-sorted 
# Filtering markers based on HWD. If a marker if in HWD in 7 reef or more with a mid - pvalue of 0.0001 then the SNP is blacklisted ==>  / 0 / 0 / 1 (/ chrom / locus / SNP) blacklisted

# After filtering
# 298 / 14 / 3131 / 69457 / 366790	==> individuals / strata / chrom / locus / markers	
# Data summary: 
# number of samples: 298
# number of markers: 366790
## The filtering data was outputted to the folder pver_filter_rad_20210222@1117

# Check missing proportion of each individual after filtering was done
pver.filtered.LD <- read_rad(data = "pver_filter_rad_20210222@1117/13_filtered/radiator_data_20210222@1527.rad")
pver.filtered.LD.gds <- tidy2gds(pver.filtered.LD) #radiator_20210224@1303.gds.rad
pver.ind.missing.prop <- detect_mixed_genomes(pver.filtered.LD.gds)
# closefn.gds("radiator_20210224@1303.gds.rad")

# Get diversity Stats of heterozygosity per Reef
pver.LD.data.sum.stats.REEF <- summary_rad(pver.data.filtered$output$tidy.data)

# Get diversity Stats of heterozygosity per Region
colnames(pver.data.filtered$output$tidy.data)[1] <- "POP_ID_REEF"
pver.data.filtered$output$tidy.data$POP_ID <- gsub(pattern = "_.*","", pver.data.filtered$output$tidy.data$POP_ID_REEF)
pver.LD.data.sum.stats.REGION <- summary_rad(pver.data.filtered$output$tidy.data)

# Filtering variants in LD
pver.LE.filtered <- filter_ld(pver.data.filtered$gds) #tidy data didn't work
# Filtering marker in LD at short distance using MAC statistics 0 / 0 / 297333 (/ chrom / locus / SNP) blacklisted
# Filtering marker in LD at long distance using Missing statistics  / 0 / 34249 / 34249 (/ chrom / locus / SNP) blacklisted

# After filtering LD
# 298 / 14 / 3131 / 35208 / 35208 (individuals / strata / chrom / locus / SNP)
## The filtering data was outputted to the folder pver_filter_ld_20210222@1736

# whitelist of markers (SNPs) in Linkage Equilibrium within each scaffold
whitelist.markers <- read.delim("./pver_filter_ld_20210222@1736/whitelist.long.ld_0.1.tsv", sep = "")
colnames(whitelist.markers)
whitelist.markers <- whitelist.markers[,c(7,2,3,6)] #MARKER, "CHROM"   "LOCUS"   "POS"  
# create tidy dataset of filtered SNPs (Linkage equilibrium)
pver.LE.filtered.tidy <- radiator::filter_whitelist(data = pver.filtered.LD, whitelist.markers = whitelist.markers)
write.table(pver.LE.filtered.tidy, file = "02_filter_ld_20210222@1736/pver.LE.filtered.tidy.tsv", sep = "\t", row.names = F, quote = F)

# Check heterozygosity and missing proportion of each individual after LD filtering was done
pver.LE.filtered.gds <- tidy2gds(pver.LE.filtered.tidy) #radiator_20210224@1422.gds.rad
pver.ind.missing.prop.LE <- detect_mixed_genomes(pver.LE.filtered.gds) 

####  01.04 Generate a filtered VCF file (bash environment) ####
# Generate list of individual to keep (298 individuals)
awk '{print $1}' 00_filter_rad_20210222@1117/13_filtered/strata.filtered.tsv | sed 1,1d | sed 's/-filtered/_filtered/' > 03_vcftools_filtered_vcf/pver.298.ind.to.keep

# Generate list of markers to keep (35,208 SNPs in LE)
awk 'BEGIN{OFS="\t";}{print $2, $6}' 02_filter_ld_20210222@1736/whitelist.long.ld_0.1.tsv | sed 1,1d  > 03_vcftools_filtered_vcf/whitelist.snps.LE

cd ~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_vcftools_filtered_vcf/
# Using VCFTOOLS filter the vcf file produced by STACKS for the filtered set of SNPs. This in this vcf file only one position in the chromosome is considered, that is if there were two SNPs called at the same position because they were cut by different restriction enzyme, STACKS will only output one SNP per site
# --ordered-export — if data is reference aligned, exports will be ordered; only a single representative of each overlapping site (stacks::populations)
ln -s ~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_filter_rad_20210222@1117/genetic_diversity_in_STACKS/For_Pi_STACKS/ByREEF/populations.snps.vcf .

vcftools --vcf populations.snps.vcf --positions whitelist.snps.LE --keep pver.298.ind.to.keep --out pver.298.LE.filtered --recode --recode-INFO-all  # After filtering, kept 35207 out of a possible 366786 Sites. One site got discarded by STACKS

# Remove the suffix "_filtered-sorted"
sed 's/_filtered-sorted//g' pver.298.LE.filtered.recode.vcf > pver.298.LE.filtered.recode.indnames.vcf # This vcf file should be use in SNPrelate

### After running SNPrelate I identify the sample that clusters with the type 7 sample and remove it from the analysis.

# Generate list of individual to keep (296 individuals) Two individuals were removed as they were type 7 
awk '{print $1}' ../03_SNPrelate_20210224/pver.new.strata.296ind.NOtype7.tsv | sed 1,1d | sed 's/-filtered-sorted//' > pver.296.ind.to.keep

# re filter the vcfile to remove the two individual samples (potential type 7)
vcftools --vcf pver.298.LE.filtered.recode.indnames.vcf --keep pver.296.ind.to.keep --out pver.LE.filtered --recode --recode-INFO-all  # After filtering, kept 296 out of 298 Individuals

# Order of individuals in vcf file (important for ADMIXTURE and PCA to define the order of individuals
grep "#CHROM" pver.LE.filtered.recode.vcf | sed 's/\t/\n/g' | sed 1,9d > pver.ind.order.in.vcf.txt

























