###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 09d Visualization of Candidate loci for Positive Selection  - Stylophora pistillata #####

# Summary of positive selection genes identified with BAYESCAN and BAYPASS

# Library
library(scales)
library(dplyr)
library(ggplot2)
source("/home/buitracn/RADseq/Genotyping/tools/manhattan.plot.function.qman.R") # I've manipulate the original qqman manhattan function to change the line types and colors


# Working directory
setwd("~/RADseq/Genotyping/p1_mac4_r80_448samples/spis_positive_selection_analyses_20210308/summary_BAYESCAN_BAYPASS")

# positive selection lists
spis.bayescan.positive <- read.delim("../BAYESCAN/spis.BAYESCAN.pr100.whitelist.markers.positive.selection.tsv", sep = "") # 76 positive candidates

spis.baypass.positive <- read.delim("../BAYPASS/BAYPASS.whitelist.markers.positive.selection.POD.calibrated.tsv", sep = "") # 111 positive candidates

spis.common.positive.candidates <- intersect(spis.baypass.positive$ID,spis.bayescan.positive$ID) # 59 positive candidates in common  

spis.all.positive.candidates.ID <- data.frame(ID=merge(spis.bayescan.positive, spis.baypass.positive , by = "ID", all = T)$ID)
spis.all.positive.candidates.ID$BAYESCAN <- factor(spis.all.positive.candidates.ID$ID %in% spis.bayescan.positive$ID, levels= c("FALSE", "TRUE"), ordered=T)
spis.all.positive.candidates.ID$BAYPASS <- factor(spis.all.positive.candidates.ID$ID %in% spis.baypass.positive$ID, levels= c("FALSE", "TRUE"), ordered=T)

write.table(spis.all.positive.candidates.ID, "spis.positive.candidates.bayescan.baypass.tsv", quote=F, row.names=F)
  
# baypass dictionary
spis.baypass.dictionary <- read.delim("../snps.id.dictionary.BAYPASS.txt", sep = "", h=F)
colnames(spis.baypass.dictionary) <- c("Scaffold", "POS", "ID", "index")

# baypass ID of all positive SNPs. This is to get the index of the markers identified by bayecan as well
all.positive.candidates.baypass.dict <- merge(spis.baypass.dictionary, spis.all.positive.candidates.ID, by = "ID")

# baypass XtX table
spis.baypass.table <- read.delim("../BAYPASS/Spis_3_summary_pi_xtx.out", sep = "")
# add scaffold info to the baypass XtX table
spis.baypass.table.scaff <- merge(spis.baypass.table, spis.baypass.dictionary, by.x = "MRK", by.y = "index")


# XtX table including both candidates identified by baypass and bayescan
spis.xtx.baypass.bayescan <- spis.baypass.table.scaff[spis.baypass.table.scaff$Scaffold %in% all.positive.candidates.baypass.dict$Scaffold, ]
spis.xtx.baypass.bayescan$Scaffold <- as.numeric(spis.xtx.baypass.bayescan$Scaffold)
spis.xtx.baypass.bayescan$POS <- as.numeric(spis.xtx.baypass.bayescan$POS)

snps.bayescan2highlight <- all.positive.candidates.baypass.dict[grep("TRUE", all.positive.candidates.baypass.dict$BAYESCAN), 4]

# Plot SNPs per scaffold (only those scaffolds with candidate SNPs for selection)
# library(qqman)
pdf("./Spis.XtX.onlyCandidateSNPsBAYPASSBAYESCAN_20210313.pdf", width = 10, height = 6)
manhattan(spis.xtx.baypass.bayescan,chr="Scaffold", 
          p="M_XtX",bp="POS",snp="MRK",logp=FALSE, 
          ylab="XtX", xlab="Scaffold ID", ylim = c(0, 30), 
          suggestiveline = 14.51123, genomewideline = F,
          highlight = snps.bayescan2highlight,
          col = alpha(c("grey47", "black"), 0.6))
dev.off()

