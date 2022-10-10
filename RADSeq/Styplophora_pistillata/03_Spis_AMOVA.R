###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 03 Analysis of Molecular Variance (AMOVA) - Stylophora pistillata #####

# Using set of unlinked SNPs (25318 SNPs - 367 ind)

# Libraries
library(adegenet)
library(poppr)
library(vcfR) 

# Working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_AMOVA_20210226")

# Create a strata to test with AMOVA
spis.strata <- read.delim("../spis_vcftools_filtered_vcf/spis.ind.order.in.vcf.txt", header = F)
colnames(spis.strata) <- "INDIVIDUALS"
spis.strata$REEF <- gsub("S","", sub('^([^-]+-[^-]+).*', '\\1', spis.strata$INDIVIDUALS))
spis.strata$REEF <- factor(spis.strata$REEF, levels = c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R3", "WAJ-R4", "YAN-R1", "YAN-R3", "YAN-R4", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R1", "FAR-R2", "FAR-R3", "FAR-R4"), ordered = T)
spis.strata$REGION <- gsub("-.*","", spis.strata$REEF)
spis.strata$REGION <- factor(spis.strata$REGION , levels = c("MAQ", "WAJ", "YAN", "KAU", "DOG", "FAR"))
strata <- spis.strata[,-1] 

# Create a genind object from a vcffile
spis.vcfR <- read.vcfR(file = "../spis_vcftools_filtered_vcf/spis.LE.filtered.recode.indnames.vcf")
spis.genind <- vcfR2genind(spis.vcfR,  ind.names= spis.strata$INDIVIDUALS, pop = spis.strata$REEF, strata = strata)


# AMOVA hierarchical
spis.amova <- poppr.amova(spis.genind, hier = ~REGION/REEF, within = F, threads = 30, method = "pegas")
# 6897 loci contained missing values greater than 5%
# Removing 6897 loci

# Percentage of variance explained by each strata
spis.amova$varcomp/sum(spis.amova$varcomp)*100
# REGION      REEF     Error 
# 4.296672  8.109608 87.593720 

# If the population was panmictic, we would expect to see extremely small variance components for Region and Reef compared to Error (i.e. the variation from individuals within populations). 

# To test for significance
spis.amova.sig <- poppr.amova(spis.genind, hier = ~REGION/REEF, within = F, threads = 30, method = "pegas", nperm = 1000)

