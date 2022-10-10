###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 04 Pairwise FST between reefs - Stylophora pistillata #####
# Using set of unlinked SNPs (25,318 SNPs - 367 ind)

library(dartR)
library(StAMPP)
library(reshape2)
library(dplyr)
library(ggplot2)
library(scales)

# Working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_PairwiseFST_20210228")

# Genind to genlight
# genind object was loaded from the AMOVA folder
spis.genind <- load("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_AMOVA_20210226/spis.genind.RData")
spis.genlight <- gi2gl(spis.genind)

# Pairwise FSt using StAMPP
spis.pairwise.fst <- stamppFst(spis.genlight, nboots = 1000, percent = 95, nclusters = 40)
class(spis.pairwise.fst$Fsts)

spis.int.conf <- spis.pairwise.fst$Bootstraps[,c(1,2,1003,1004,1005,1006)]
spis.int.conf$CI.95 <- paste(format(round(spis.int.conf$`Lower bound CI limit`,3),3), format(round(spis.int.conf$`Upper bound CI limit`, 3), 3), sep = "-")

spis.int.conf$Fst <- format(round(spis.int.conf$Fst, 3), 3)

spis.int.conf$Population1 <- factor(spis.int.conf$Population1, levels=c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R3", "WAJ-R4", "YAN-R1", "YAN-R3", "YAN-R4", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R1", "FAR-R2", "FAR-R3", "FAR-R4"), ordered = T)

spis.int.conf$Population2 <- factor(spis.int.conf$Population2, levels=c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R3", "WAJ-R4", "YAN-R1", "YAN-R3", "YAN-R4", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R1", "FAR-R2", "FAR-R3", "FAR-R4"), ordered = T)

spis.ci.matrix <- acast(spis.int.conf, Population1~Population2, value.var="CI.95")
write.table(spis.ci.matrix, file = "spis.fst.ci95.1000perm.tsv", quote = F, sep = "\t")

spis.fst.matrix <- t(acast(spis.int.conf, Population1~Population2, value.var="Fst"))
write.table(spis.fst.matrix, file = "spis.fst.value.tsv", quote = F, sep = "\t")


# complete the symetric matrix (better visual representation)
spis.pairwise.fst.low.tri <- spis.pairwise.fst$Fsts
spis.pairwise.fst.upper.tri <- t(spis.pairwise.fst.low.tri)

spis.pairwise.fst.complete <- matrix(NA, nrow = 18, ncol = 18)
spis.pairwise.fst.complete[upper.tri(spis.pairwise.fst.complete)] <- spis.pairwise.fst.upper.tri[upper.tri(spis.pairwise.fst.upper.tri)]
spis.pairwise.fst.complete[lower.tri(spis.pairwise.fst.complete)] <- spis.pairwise.fst.low.tri[lower.tri(spis.pairwise.fst.low.tri)]

diag(spis.pairwise.fst.complete) <- 0 # fill the diagonal with 0
rownames(spis.pairwise.fst.complete) <- rownames(spis.pairwise.fst.low.tri)
colnames(spis.pairwise.fst.complete) <- colnames(spis.pairwise.fst.low.tri)

# Format data for heatmap
spis.melted.Allmarkers.fst <- melt(spis.pairwise.fst.complete, na.rm =TRUE)
spis.melted.Allmarkers.fst$value <- as.numeric(levels(spis.melted.Allmarkers.fst$value))[spis.melted.Allmarkers.fst$value]

#reverese the order of the populations to preverse north to south order
spis.melted.Allmarkers.fst.rev <- spis.melted.Allmarkers.fst %>%
  # convert state to factor and reverse order of levels
  mutate(Var1=factor(Var1,levels=rev(sort(unique(Var1)))))


mid <- median(spis.melted.Allmarkers.fst.rev$value)

pdf('Spis_FST_25kSNPs_heatmap_3colors.pdf', width = 6.5, height = 6)
ggplot(data = spis.melted.Allmarkers.fst.rev, aes(Var1, Var2, fill = value))+ geom_tile(color = "white")+ 
  scale_fill_gradient2(midpoint = 0.09, low="#0EC4F9", mid="#FFD700", high="#FF6347", limits=range(spis.melted.Allmarkers.fst.rev$value), name="FST")  +
  ggtitle(expression(atop("Stylophora pistillata - Pairwise FST, WC (1984)", atop(italic("N = 367, L = 25,318"), ""))))+
  theme(plot.background=element_blank(), panel.border=element_blank(), axis.line = element_line())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),axis.text.y = element_text(size = 13), axis.title = element_blank()) + 
  theme(legend.text = element_text(size =11), legend.title = element_text(size =12)) +
  theme(plot.title = element_text(size = 14)) +
  coord_flip() 
dev.off()

