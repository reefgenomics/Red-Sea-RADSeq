###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 03 Analysis of Molecular Variance (AMOVA) - Pocillopora verrucosa #####

# Using set of unlinked SNPs (35,208 SNPs - 296)

# Libraries
library(adegenet)
library(poppr)
library(vcfR) 

# Working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_AMOVA_20210226")

# Create a strata to test with AMOVA
pver.strata <- read.delim("../pver_vcftools_filtered_vcf/pver.ind.order.in.vcf.txt", header = F)
colnames(pver.strata) <- "INDIVIDUALS"
pver.strata$REEF <- gsub("P","", sub('^([^-]+-[^-]+).*', '\\1', pver.strata$INDIVIDUALS))
pver.strata$REEF <- factor(pver.strata$REEF, levels = c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R2", "WAJ-R3", "YAN-R1", "YAN-R3", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R3"), ordered = T)
pver.strata$REGION <- gsub("-.*","", pver.strata$REEF)
pver.strata$REGION <- factor(pver.strata$REGION , levels = c("MAQ", "WAJ", "YAN", "KAU", "DOG", "FAR"))
strata <- pver.strata[,-1] 
  
# Create a genind object from a vcffile
pver.vcfR <- read.vcfR(file = "../pver_vcftools_filtered_vcf/pver.LE.filtered.recode.vcf")
pver.genind <- vcfR2genind(pver.vcfR,  ind.names= pver.strata$INDIVIDUALS, pop = pver.strata$REEF, strata = strata)


# AMOVA hierarchical
pver.amova <- poppr.amova(pver.genind, hier = ~REGION/REEF, within = F, threads = 30, method = "pegas")
# 7203 loci contained missing values greater than 5%
# Removing 7203 loci

# Percentage of variance explained by each strata
pver.amova$varcomp/sum(pver.amova$varcomp)*100
# REGION        REEF       Error 
# 0.28483466  0.02829822 99.68686712

# If the population was panmictic, we would expect to see extremely small variance components for Region and Reef compared to Error (i.e. the variation from individuals within populations). 

# To test for significance
pver.amova.sig <- poppr.amova(pver.genind, hier = ~REGION/REEF, within = F, threads = 30, method = "pegas", nperm = 1000)

