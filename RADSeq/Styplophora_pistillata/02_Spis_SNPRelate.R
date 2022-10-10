###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 02 Clustering analysis of filtered data with SNPrelate - Stylophora pistillata #####

# library(devtools)
# devtools::install_github("thierrygosselin/radiator", dependencies = TRUE)
library(radiator)
library(SNPRelate)
library(SeqArray)
library(gdsfmt)
library(scales)
library(dplyr) 
library(dendextend) #easy manypulation of dendograms

# Set working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_SNPrelate_20210222")

# Check the filtered strata file
filtered.strata.spis <- radiator::read_strata(strata = "../spis_filter_rad_20210221@1444/13_filtered/strata.filtered.tsv")
radiator::summary_strata(strata = filtered.strata.spis$strata)
# Number of populations: 18
# Number of individuals: 367

# The VCF file, using the example in the SNPRelate package
vcf.fn <- "../spis_vcftools_filtered_vcf/spis.LE.filtered.recode.indnames.vcf"

# Reformat
snpgdsVCF2GDS(vcf.fn, "spis.LE.gds", method="biallelic.only")

# Generating SNPRelate object/file...
# SNPRelate GDS: spis.filtered.18reef_snprelate_20200305@1132.gds.rad
# once the gds file has been created one can open it by typing
spis.18reef.snprelate <- SNPRelate::snpgdsOpen("spis.LE.gds", readonly = FALSE)

# Inspect the structure of the GDS file
spis.18reef.snprelate

# Check the order of the samples in the gds file
sample.id.18reef <- read.gdsn(index.gdsn(spis.18reef.snprelate, "sample.id"))
# gds file was created with samples in the order of the vcf file

# Add samples annotation to the newly created GDS file
# create element to add to the gds file. 
samp.annot <- data.frame(reef.group = filtered.strata.spis$strata$STRATA,
                         loc.group = gsub("_.*", "", filtered.strata.spis$strata$STRATA))


#add samp.annot element to the gds object
add.gdsn(spis.18reef.snprelate, "sample.annot", samp.annot)


# Read population information
reef.code.18reef <- read.gdsn(index.gdsn(spis.18reef.snprelate, path="sample.annot/reef.group"))
loc.code.18reef <- read.gdsn(index.gdsn(spis.18reef.snprelate, path="sample.annot/loc.group"))
table(reef.code.18reef)
table(loc.code.18reef)

# cbind.data.frame(sample.id.18reef, reef.code.18reef, loc.code.18reef)

####  02.01 Principal Component Analysis
pca.spis.18reef <- snpgdsPCA(spis.18reef.snprelate, autosome.only = FALSE, num.thread=2)

# variance proportion explained (%)
pc.percent.18reef <- pca.spis.18reef$varprop*100
head(round(pc.percent.18reef, 2))
# [1] 3.79 2.73 2.14 2.01 1.87 1.77
# This is considerably very little variance explained by each of the PCA axis

# Create a data.frame
sampID.annot <-data.frame(sample.id = sample.id.18reef[match(pca.spis.18reef$sample.id, sample.id.18reef)],
                          reef = factor(reef.code.18reef, levels = c("MAQ_R1", "MAQ_R2", "WAJ_R1", "WAJ_R3", "WAJ_R4", "YAN_R1", "YAN_R3", "YAN_R4", "KAU_R1", "KAU_R2", "KAU_R3", "DOG_R1", "DOG_R2", "DOG_R3", "FAR_R1", "FAR_R2", "FAR_R3", "FAR_R4"))[match(pca.spis.18reef$sample.id, sample.id.18reef)],
                          location = factor(loc.code.18reef, levels = c("MAQ", "WAJ", "YAN", "KAU", "DOG", "FAR"))[match(pca.spis.18reef$sample.id, sample.id.18reef)])


# Visualize the first 4 PCs

# Set color pallete for the 18 reefs
cols.18reef <- c("#0066CC","#3399FF","#009999", "#00CCCC","#33FFFF","#009900","#66CC00","#80FF00","#CCCC00","#FFFF00","#FFFF66","#D2691E","#FF8000","#FFB266","#CC0000","#FF3333","#FF6666","#FF9999" )
show_col(cols.18reef)

# Create a color vector corresponding to levels in the reef vector
cols.18reef_t1 <- cols.18reef[sampID.annot$reef]

# R plot by reefs
pdf(file = "PCA_367samplesbyreef.pdf", width = 7, height = 7)
lbls.18reef <- paste("PC", 1:4, "\n", format(pc.percent.18reef[1:4], digits=2), "%", sep="")
pairs(pca.spis.18reef$eigenvect[,1:4], labels=lbls.18reef, col=cols.18reef_t1, oma=c(3,3,3,7))
par(xpd = T, mar = par()$mar + c(0,0,0,7))
legend("topright",legend=levels(sampID.annot$reef), pch=19, pt.cex=1, col=cols.18reef, cex = 0.5, inset=c(-0.02,0), box.lty=0, y.intersp=2, x.intersp =0.5, bty='n', box.col = "white", bg = "white")
par(mar=c(3, 3, 3, 3) + 0.1)
dev.off()

# Set a color palette for locations
cols.6loc <- c("#0066CC", "#00CCCC", "#66CC00", "#FFFF00", "#FF8000", "#FF3333")

# Create a color vector corresponding to levels in the reef vector
cols.6loc_t1 <- cols.6loc[sampID.annot$location]

# R plot by latitudinal location
pdf(file = "PCA_367samplesbylocation.pdf", width = 7, height = 7)
lbls.18reef <- paste("PC", 1:4, "\n", format(pc.percent.18reef[1:4], digits=2), "%", sep="")
pairs(pca.spis.18reef$eigenvect[,1:4], labels=lbls.18reef, col=cols.6loc_t1, oma=c(3,3,3,7))
par(xpd = T, mar = par()$mar + c(0,0,0,7))
legend("topright",legend=levels(sampID.annot$location), pch=19, pt.cex=1, col=cols.6loc, cex = 0.5, inset=c(-0.02,0), box.lty=0, y.intersp=2, x.intersp =0.5, bty='n', box.col = "white", bg = "white")
par(mar=c(3, 3, 3, 3) + 0.1)
dev.off()

# There is not a real pattern observed in the PCA. Diagonals clines can be interpret as a signature of a hierarchical clustering (by pcaadapt authors)

####  02.02 Identity by State Analysis
spis.ibs.18reef <- SNPRelate::snpgdsIBS(
  gdsobj = spis.18reef.snprelate,
  autosome.only = FALSE,
  num.thread = 7,
  remove.monosnp = TRUE,
  verbose = TRUE)

spis.matrix <- spis.ibs.18reef$ibs
spis.ibs.18reef[["sample.id"]]=gsub("-filtered-sorted", "",spis.ibs.18reef[["sample.id"]])
spis.ibs.18reef[["sample.id"]]=gsub("-", "_",spis.ibs.18reef[["sample.id"]])
colnames(spis.matrix )=spis.ibs.18reef[["sample.id"]]
rownames(spis.matrix )=spis.ibs.18reef[["sample.id"]]
write.table(spis.matrix, file = "spis.ind.distance.matrix.IBS.snprelate.tsv", row.names = T, col.names = T, quote = F)
  
# multidimensional scaling analysis on the n×n matrix of genome-wide IBS pairwise distances
loc.18reef <- cmdscale(1 - spis.ibs.18reef$ibs, k = 2)
x <- loc.18reef[, 1]; y <- loc.18reef[, 2]

# R plot by sampling site (reef)
pdf(file = "MDS_367samplesbyreef.pdf", width = 7, height = 7)
par(mar=c(3, 3, 3, 5), xpd=TRUE)
plot(x, y, col=cols.18reef_t1, xlab = "", ylab = "",
     main = "Stylophora pistillata - Multidimensional Scaling Analysis (MDS)",
     cex.main=1)
legend("topright",legend=levels(sampID.annot$reef), pch=19, pt.cex=1, col=cols.18reef, cex = 0.6, inset=c(-0.15,0), box.lty=0, y.intersp=2, x.intersp =1, bty='n', box.col = "white",bg = "white")
dev.off()

# R plot by latitudinal location (region)
pdf(file = "MDS_367samplesbylocation.pdf", width = 7, height = 7)
par(mar=c(3, 3, 3, 5), xpd=TRUE)
plot(x, y, col=cols.6loc_t1, xlab = "", ylab = "",
     main = "Stylophora pistillata - Multidimensional Scaling Analysis (MDS)",
     cex.main=1)
legend("topright",legend=levels(sampID.annot$location), pch=19, pt.cex=1, col=cols.6loc, cex = 0.6, inset=c(-0.15,0), box.lty=0, y.intersp=2, x.intersp =1, bty='n', box.col = "white",bg = "white")
dev.off()

# Clusters are not clearly defined by the reefs nor regions

# To perform cluster analysis on the n×n matrix of genome-wide IBS pairwise distances, and determine the groups by a permutation score
# snpgdsHCluster calls the function hclust to perform hierarchical cluster analysis, using method="average". The distance between two groups is defined as the average distance between each of their members
set.seed(100)
spis.ibs.18reef.hc <- snpgdsHCluster(
  snpgdsIBS(
    spis.18reef.snprelate,
    autosome.only = FALSE,
    num.thread=2))

# Determine groups of individuals by population information
rv2 <- snpgdsCutTree(spis.ibs.18reef.hc, samp.group=factor(reef.code.18reef, levels = c("MAQ_R1", "MAQ_R2", "WAJ_R1", "WAJ_R3", "WAJ_R4", "YAN_R1", "YAN_R3", "YAN_R4", "KAU_R1", "KAU_R2", "KAU_R3", "DOG_R1", "DOG_R2", "DOG_R3", "FAR_R1", "FAR_R2", "FAR_R3", "FAR_R4")))
# Create 18 groups.
names(rv2)
# [1] "sample.id"   "z.threshold" "outlier.n"   "samp.order"  "samp.group" 
# [6] "dmat"        "dendrogram" 

# Define the IBS based dendrogram as an object 
spis.ibs.18reef.dend <- as.dendrogram(rv2$dendrogram)

# Assigning the labels of the dendrogram object with new colors:
labels_colors(spis.ibs.18reef.dend) <- cols.18reef[sampID.annot$reef][order.dendrogram(spis.ibs.18reef.dend)]

# Dendrogram with labels for each individual 
pdf(file = "Identity-By-State-analysis-genotypes-18reefwithIDnames.pdf", width = 20, height = 6)
par(mar=c(6, 3, 3, 3), xpd=TRUE)
spis.ibs.18reef.dend %>%
  set("labels_cex", 0.4)%>%
  set("labels_col", value=alpha(cols.18reef[sampID.annot$reef][order.dendrogram(spis.ibs.18reef.dend)], alpha = 1))%>%
  set("leaves_pch", 19)%>% #symbol of the end of each leaf
  set("leaves_col", value=alpha(cols.18reef[sampID.annot$reef][order.dendrogram(spis.ibs.18reef.dend)], alpha = 1))%>% #assign color to each leaf
  #set("branches_k_color", k=5)%>% #assign colors to n number of branches
  hang.dendrogram(hang = -0.1)%>%
  plot(main="Hierarchical clustering of Stylophora pistillata from the Red Sea based on identity by state (IBS)", cex.main=1, ylab = "IBS based distance", ylim = c(0,0.15))
legend("topright",legend=levels(sampID.annot$reef), pch=19, col=cols.18reef, ncol = 1, cex = 0.7, inset=c(-0,0), box.lty=0, y.intersp=0.8, x.intersp =0.5, bty='n', box.col = "white",bg = "white")
dev.off()

# dendrograms often suggest a correct number of clusters when there is no real evidence to support the conclusion (https://www.displayr.com/what-is-dendrogram/)

# The reefs with fewer individuals are KAU_R1 (8), WAJ_R1 (12), DOG_R1 (14), but they cluster together with nearby reef, therefore it is not worth removing them from the analysis

# To close the GDS file
# SNPRelate::snpgdsClose or
closefn.gds(spis.18reef.snprelate)

