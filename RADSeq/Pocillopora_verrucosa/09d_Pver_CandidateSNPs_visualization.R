###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 09d Visualization of Candidate loci for Positive Selection  - Pocillopora verrucosa #####

# Summary of positive selection genes identified with BAYESCAN and BAYPASS

# Library
library(scales)
library(dplyr)
library(ggplot2)
source("/home/buitracn/RADseq/Genotyping/tools/manhattan.plot.function.qman.R") # I've manipulate the original qqman manhattan function to change the line types and colors


# Working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_positive_selection_analyses_20210308/summary_BAYESCAN_BAYPASS")

# positive selection lists
pver.bayescan.positive <- read.delim("../BAYESCAN/pver.BAYESCAN.pr100.whitelist.markers.positive.selection.tsv", sep = "") # 22 positive candidates

pver.baypass.positive <- read.delim("../BAYPASS/BAYPASS.whitelist.markers.positive.selection.POD.calibrated.tsv", sep = "") # 85 positive candidates (in 81scaffolds)

pver.common.positive.candidates <- intersect(pver.baypass.positive$ID,pver.bayescan.positive$ID) # 22 positive candidates in common  

pver.all.positive.candidates.ID <- data.frame(ID=merge(pver.bayescan.positive, pver.baypass.positive , by = "ID", all = T)$ID)
pver.all.positive.candidates.ID$BAYESCAN <- factor(pver.all.positive.candidates.ID$ID %in% pver.bayescan.positive$ID, levels= c("FALSE", "TRUE"), ordered=T)
pver.all.positive.candidates.ID$BAYPASS <- factor(pver.all.positive.candidates.ID$ID %in% pver.baypass.positive$ID, levels= c("FALSE", "TRUE"), ordered=T)

write.table(pver.all.positive.candidates.ID, "pver.positive.candidates.bayescan.baypass.tsv", quote=F, row.names=F)
  
# baypass dictionary
pver.baypass.dictionary <- read.delim("../snps.id.dictionary.BAYPASS.txt", sep = "", h=F)
colnames(pver.baypass.dictionary) <- c("Scaffold", "POS", "ID", "index")

# baypass ID of all positive SNPs. This is to get the index of the markers identified by bayecan as well
all.positive.candidates.baypass.dict <- merge(pver.baypass.dictionary, pver.all.positive.candidates.ID, by = "ID")

# baypass XtX table
pver.baypass.table <- read.delim("../BAYPASS/Pver_3_summary_pi_xtx.out", sep = "")
# add scaffold info to the baypass XtX table
pver.baypass.table.scaff <- merge(pver.baypass.table, pver.baypass.dictionary, by.x = "MRK", by.y = "index")


# XtX table including both candidates identified by baypass and bayescan
pver.xtx.baypass.bayescan <- pver.baypass.table.scaff[pver.baypass.table.scaff$Scaffold %in% all.positive.candidates.baypass.dict$Scaffold, ]
pver.xtx.baypass.bayescan$Scaffold <- as.numeric(gsub("Pver_xfSc", "", gsub("Pver_xpSc", "", pver.xtx.baypass.bayescan$Scaffold))) # define the vector of scaffold IDs as numeric
pver.xtx.baypass.bayescan$POS <- as.numeric(pver.xtx.baypass.bayescan$POS)

snps.bayescan2highlight <- all.positive.candidates.baypass.dict[grep("TRUE", all.positive.candidates.baypass.dict$BAYESCAN), 4]

# Plot SNPs per scaffold (only those scaffolds with candidate SNPs for selection)
# library(qqman)
pdf("./pver.XtX.onlyCandidateSNPsBAYPASSBAYESCAN_20210313.pdf", width = 10, height = 6)
manhattan(pver.xtx.baypass.bayescan,chr="Scaffold", 
          p="M_XtX",bp="POS",snp="MRK",logp=FALSE, 
          ylab="XtX", xlab="Scaffold ID", ylim = c(0, 30), 
          suggestiveline = 5.515877, genomewideline = F,
          highlight = snps.bayescan2highlight,
          col = alpha(c("grey47", "black"), 0.6))
dev.off()

