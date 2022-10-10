###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 04 Pairwise FST between reefs - Pocillopora verrucosa #####
# Using set of unlinked SNPs (35,208 SNPs - 296 ind)

library(dartR)
library(StAMPP)
library(reshape2)
library(dplyr)
library(ggplot2)
library(scales)

# Working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_PairwiseFST_20210228")

# Genind to genlight (The genind object generated for AMOVa was used here)
pver.genlight <- gi2gl(pver.genind)

# Pairwise FSt using StAMPP
pver.pairwise.fst <- stamppFst(pver.genlight, nboots = 1000, percent = 95, nclusters = 40)
class(pver.pairwise.fst$Fsts)

pver.int.conf <- pver.pairwise.fst$Bootstraps[,c(1,2,1003,1004,1005,1006)]
pver.int.conf$CI.95 <- paste(format(round(pver.int.conf$`Lower bound CI limit`,3),3), format(round(pver.int.conf$`Upper bound CI limit`, 3), 3), sep = "-")

pver.int.conf$Fst <- format(round(pver.int.conf$Fst, 3), 3)

pver.int.conf$Population1 <- factor(pver.int.conf$Population1, levels=c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R2", "WAJ-R3", "YAN-R1", "YAN-R3", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R3"), ordered = T)

pver.int.conf$Population2 <- factor(pver.int.conf$Population2, levels=c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R2", "WAJ-R3", "YAN-R1", "YAN-R3", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R3"), ordered = T)

pver.ci.matrix <- acast(pver.int.conf, Population1~Population2, value.var="CI.95")
write.table(pver.ci.matrix, file = "pver.fst.ci95.1000perm.tsv", quote = F, sep = "\t")

pver.fst.matrix <- t(acast(pver.int.conf, Population1~Population2, value.var="Fst"))
write.table(pver.fst.matrix, file = "pver.fst.value.tsv", quote = F, sep = "\t")



# complete the symetric matrix (better visual representation)
pver.pairwise.fst.low.tri <- pver.pairwise.fst$Fsts
pver.pairwise.fst.upper.tri <- t(pver.pairwise.fst.low.tri)

pver.pairwise.fst.complete <- matrix(NA, nrow = 14, ncol = 14)
pver.pairwise.fst.complete[upper.tri(pver.pairwise.fst.complete)] <- pver.pairwise.fst.upper.tri[upper.tri(pver.pairwise.fst.upper.tri)]
pver.pairwise.fst.complete[lower.tri(pver.pairwise.fst.complete)] <- pver.pairwise.fst.low.tri[lower.tri(pver.pairwise.fst.low.tri)]

diag(pver.pairwise.fst.complete) <- 0 # fill the diagonal with 0
rownames(pver.pairwise.fst.complete) <- rownames(pver.pairwise.fst.low.tri)
colnames(pver.pairwise.fst.complete) <- colnames(pver.pairwise.fst.low.tri)

# Format data for heatmap
pver.melted.Allmarkers.fst <- melt(pver.pairwise.fst.complete, na.rm =TRUE)
pver.melted.Allmarkers.fst$value <- as.numeric(levels(pver.melted.Allmarkers.fst$value))[pver.melted.Allmarkers.fst$value]

#reverese the order of the populations to preverse north to south order
pver.melted.Allmarkers.fst.rev <- pver.melted.Allmarkers.fst %>%
  # convert state to factor and reverse order of levels
  mutate(Var1=factor(Var1,levels=rev(sort(unique(Var1)))))


mid <- median(pver.melted.Allmarkers.fst.rev$value)


pdf('Pver_FST_35kSNPs_heatmap_samespisscale.pdf', width = 6.5, height = 6)
ggplot(data = pver.melted.Allmarkers.fst.rev, aes(Var1, Var2, fill = value))+ geom_tile(color = "white")+ 
  scale_fill_gradient(low = "#0EC4F9", high = "#a9cfdc", limits=range(pver.melted.Allmarkers.fst.rev$value), name="FST")  +
  ggtitle(expression(atop("Pocillopora verrucosa - Pairwise FST, WC (1984)", atop(italic("N = 296, L = 35,207"), ""))))+
  theme(plot.background=element_blank(), panel.border=element_blank(), axis.line = element_line())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),axis.text.y = element_text(size = 13), axis.title = element_blank()) + 
  theme(legend.text = element_text(size =11), legend.title = element_text(size =12)) +
  theme(plot.title = element_text(size = 14)) +
  coord_flip() 
dev.off()

