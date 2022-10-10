###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 02 Clustering analysis of filtered data with SNPrelate - Pocillopora verrucosa #####

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
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_SNPrelate_20210224")

# Check the filtered strata file
filtered.strata.pver <- radiator::read_strata(strata = "../00_filter_rad_20210222@1117/13_filtered/strata.filtered.tsv")
radiator::summary_strata(strata = filtered.strata.pver$strata)
# Number of populations: 18
# Number of individuals: 298

# The VCF file, using the example in the SNPRelate package
vcf.fn <- "../pver_vcftools_filtered_vcf/pver.298.LE.filtered.recode.indnames.vcf"

# Reformat
snpgdsVCF2GDS(vcf.fn, "pver.LE.gds", method="biallelic.only")


# Generating SNPRelate object/file...
# SNPRelate GDS: pver.filtered.14reef_snprelate_20200305@1132.gds.rad
# once the gds file has been created one can open it by typing
pver.14reef.snprelate <- SNPRelate::snpgdsOpen("pver.LE.gds", readonly = FALSE)

# Inspect the structure of the GDS file
pver.14reef.snprelate
SNPRelate::snpgdsSummary(pver.14reef.snprelate)
# The file name: /home/buitracn/RADseq-Big-project/pver/p1.mac4.r0.8.316samples.Pver_20200218/03_SNPrelate_20210224/pver.LE.gds 
# The total number of samples: 298 
# The total number of SNPs: 35208 
# SNP genotypes are stored in SNP-major mode (Sample X SNP).

# Check the order of the samples in the gds file
sample.id.14reef <- read.gdsn(index.gdsn(pver.14reef.snprelate, "sample.id"))
# gds file was created with samples in the order of the vcf file

# Add samples annotation to the newly created GDS file
# create element to add to the gds file. 
samp.annot <- data.frame(reef.group = filtered.strata.pver$strata$STRATA,
                         loc.group = factor(gsub("_.*", "", filtered.strata.pver$strata$STRATA), levels = c("MAQ", "WAJ", "YAN", "KAU", "DOG", "FAR"), ordered = T))


#add samp.annot element to the gds object
add.gdsn(pver.14reef.snprelate, "sample.annot", samp.annot)


# Read population information
reef.code.14reef <- read.gdsn(index.gdsn(pver.14reef.snprelate, path="sample.annot/reef.group"))
loc.code.14reef <- read.gdsn(index.gdsn(pver.14reef.snprelate, path="sample.annot/loc.group"))
table(reef.code.14reef)
table(loc.code.14reef)

# cbind.data.frame(sample.id.14reef, reef.code.14reef, loc.code.14reef)


####  02.01 Principal Component Analysis
pca.pver.14reef <- snpgdsPCA(pver.14reef.snprelate, autosome.only = FALSE, num.thread=2)

# variance proportion explained (%)
pc.percent.14reef <- pca.pver.14reef$varprop*100
head(round(pc.percent.14reef, 2))
# [1] 0.62 0.42 0.42 0.42 0.41 0.41
# This is considerably very little variance explained by each of the PCA axis

# Create a data.frame
sampID.annot <-data.frame(sample.id = sample.id.14reef[match(pca.pver.14reef$sample.id, sample.id.14reef)],
                          reef = factor(reef.code.14reef, levels = c("MAQ_R1", "MAQ_R2", "WAJ_R1", "WAJ_R2", "WAJ_R3", "YAN_R1", "YAN_R3", "KAU_R1", "KAU_R2", "KAU_R3", "DOG_R1", "DOG_R2", "DOG_R3", "FAR_R3"))[match(pca.pver.14reef$sample.id, sample.id.14reef)],
                          location = factor(loc.code.14reef, levels = c("MAQ", "WAJ", "YAN", "KAU", "DOG", "FAR"))[match(pca.pver.14reef$sample.id, sample.id.14reef)])

# Visualize the first 4 PCs

# Set up a color pallete for the 14 reefs
cols.14reef <- c("#0066CC","#3399FF","#009999", "#27C3AE", "#00CCCC", "#009900","#66CC00", "#CCCC00","#FFFF00","#FFFF66","#D2691E","#FF8000","#FFB266","#FF6666")
show_col(cols.14reef)
# Create a color vector corresponding to levels in the reef vector
cols.14reef_t1 <- cols.14reef[sampID.annot$reef]

# R plot by reefs
pdf(file = "PCA_298samplesbyreef.pdf", width = 7, height = 7)
lbls.14reef <- paste("PC", 1:4, "\n", format(pc.percent.14reef[1:4], digits=2), "%", sep="")
pairs(pca.pver.14reef$eigenvect[,1:4], labels=lbls.14reef, col=cols.14reef_t1, oma=c(3,3,3,7))
par(xpd = T, mar = par()$mar + c(0,0,0,7))
legend("topright",legend=levels(sampID.annot$reef), pch=19, pt.cex=1, col=cols.14reef, cex = 0.5, inset=c(-0.02,0), box.lty=0, y.intersp=2, x.intersp =0.5, bty='n', box.col = "white", bg = "white")
par(mar=c(3, 3, 3, 3) + 0.1)
dev.off()

# Set a color palette for locations
cols.6loc <- c("#0066CC", "#00CCCC", "#66CC00", "#FFFF00", "#FF8000", "#FF3333")

# Create a color vector corresponding to levels in the reef vector
cols.6loc_t1 <- cols.6loc[sampID.annot$location]

# R plot by latitudinal location
pdf(file = "PCA_298samplesbylocation.pdf", width = 7, height = 7)
lbls.14reef <- paste("PC", 1:4, "\n", format(pc.percent.14reef[1:4], digits=2), "%", sep="")
pairs(pca.pver.14reef$eigenvect[,1:4], labels=lbls.14reef, col=cols.6loc_t1, oma=c(3,3,3,7))
par(xpd = T, mar = par()$mar + c(0,0,0,7))
legend("topright",legend=levels(sampID.annot$location), pch=19, pt.cex=1, col=cols.6loc, cex = 0.5, inset=c(-0.02,0), box.lty=0, y.intersp=2, x.intersp =0.5, bty='n', box.col = "white", bg = "white")
par(mar=c(3, 3, 3, 3) + 0.1)
dev.off()

# Individuals from all populations are indistinguishable except in the southern most reef)

####  02.02 Identity by State Analysis
pver.ibs.14reef <- SNPRelate::snpgdsIBS(
  gdsobj = pver.14reef.snprelate,
  autosome.only = FALSE,
  num.thread = 7,
  remove.monosnp = TRUE,
  verbose = TRUE)

pver.matrix <- pver.ibs.14reef$ibs
pver.ibs.14reef[["sample.id"]]=gsub("-filtered-sorted", "",pver.ibs.14reef[["sample.id"]])
pver.ibs.14reef[["sample.id"]]=gsub("-", "_",pver.ibs.14reef[["sample.id"]])
colnames(pver.matrix )=pver.ibs.14reef[["sample.id"]]
rownames(pver.matrix )=pver.ibs.14reef[["sample.id"]]
write.table(pver.matrix, file = "pver.ind.distance.matrix.IBS.snprelate.tsv", row.names = T, col.names = T, quote = F)
  
# multidimensional scaling analysis on the n×n matrix of genome-wide IBS pairwise distances
loc.14reef <- cmdscale(1 - pver.ibs.14reef$ibs, k = 2)
x <- loc.14reef[, 1]; y <- loc.14reef[, 2]

# R plot by reef
pdf(file = "MDS_298samplesbyreef.pdf", width = 7, height = 7)
par(mar=c(3, 3, 3, 5), xpd=TRUE)
plot(x, y, col=cols.14reef_t1, xlab = "", ylab = "",
     main = "Stylophora pistillata - Multidimensional Scaling Analysis (MDS)",
     cex.main=1)
legend("topright",legend=levels(sampID.annot$reef), pch=19, pt.cex=1, col=cols.14reef, cex = 0.6, inset=c(-0.15,0), box.lty=0, y.intersp=2, x.intersp =1, bty='n', box.col = "white",bg = "white")
dev.off()

# R plot by latitudinal location
pdf(file = "MDS_298samplesbylocation.pdf", width = 7, height = 7)
par(mar=c(3, 3, 3, 5), xpd=TRUE)
plot(x, y, col=cols.6loc_t1, xlab = "", ylab = "",
     main = "Stylophora pistillata - Multidimensional Scaling Analysis (MDS)",
     cex.main=1)
legend("topright",legend=levels(sampID.annot$location), pch=19, pt.cex=1, col=cols.6loc, cex = 0.6, inset=c(-0.15,0), box.lty=0, y.intersp=2, x.intersp =1, bty='n', box.col = "white",bg = "white")
dev.off()


# Only reef FAR-R3 was different from the rest of the reefs

# To perform cluster analysis on the n×n matrix of genome-wide IBS pairwise distances, and determine the groups by a permutation score
# snpgdsHCluster calls the function hclust to perform hierarchical cluster analysis, using method="average". The distance between two groups is defined as the average distance between each of their members
set.seed(100)
pver.ibs.14reef.hc <- snpgdsHCluster(
  snpgdsIBS(
    pver.14reef.snprelate,
    autosome.only = FALSE,
    num.thread=2))

# Determine groups of individuals by population information
rv2 <- snpgdsCutTree(pver.ibs.14reef.hc, samp.group=factor(reef.code.14reef, levels = c("MAQ_R1", "MAQ_R2", "WAJ_R1", "WAJ_R2", "WAJ_R3", "YAN_R1", "YAN_R3",  "KAU_R1", "KAU_R2", "KAU_R3", "DOG_R1", "DOG_R2", "DOG_R3", "FAR_R3")))
# Create 14 groups.
names(rv2)
# [1] "sample.id"   "z.threshold" "outlier.n"   "samp.order"  "samp.group" 
# [6] "dmat"        "dendrogram" 

# Define the IBS based dendrogram as an object 
pver.ibs.14reef.dend <- as.dendrogram(rv2$dendrogram)
 
# Assigning the labels of the dendrogram object with new colors:
labels_colors(pver.ibs.14reef.dend) <- cols.14reef[sampID.annot$reef][order.dendrogram(pver.ibs.14reef.dend)]

# Dendrogram with labels for each individual 
pdf(file = "Identity-By-State-analysis-genotypes-14reefwithIDnames.pdf", width = 20, height = 6)
par(mar=c(6, 3, 3, 3), xpd=TRUE)
pver.ibs.14reef.dend %>%
  set("labels_cex", 0.4)%>%
  set("labels_col", value=alpha(cols.14reef[sampID.annot$reef][order.dendrogram(pver.ibs.14reef.dend)], alpha = 1))%>%
  set("leaves_pch", 19)%>% #symbol of the end of each leaf
  set("leaves_col", value=alpha(cols.14reef[sampID.annot$reef][order.dendrogram(pver.ibs.14reef.dend)], alpha = 1))%>% #assign color to each leaf
  #set("branches_k_color", k=5)%>% #assign colors to n number of branches
  hang.dendrogram(hang = -0.1)%>%
  plot(main="Hierarchical clustering of Pocillopora verrucosa from the Red Sea based on identity by state (IBS)", cex.main=1,  ylim = c(0,0.15), ylab="IBS-based distance")
legend("topright",legend=levels(sampID.annot$reef), pch=19, col=cols.14reef, ncol = 1, cex = 0.7, inset=c(-0,0), box.lty=0, y.intersp=0.8, x.intersp =0.5, bty='n', box.col = "white",bg = "white")
dev.off()

# dendrograms often suggest a correct number of clusters when there is no real evidence to support the conclusion (https://www.displayr.com/what-is-dendrogram/)

### sample PFAR−R3−05 (type 7) and PFAR−R3−25 (that clusters with type 7) should be removed from further analyses
# load the Linkage equilibrium filtered file (14 reefs - 298 inds) and blacklist PFAR−R3−05 (type 7) and PFAR−R3−25 have to be blacklisted
pver.14reef.LE.filtered.tidy <- read.delim("../pver_filter_ld_20210222@1736/pver.LE.filtered.tidy.tsv", sep = "", header = T)

# remove outlier individuals  (14 reefs - 296 inds)
pver.14reef.LE.filtered.noType7 <-  pver.14reef.LE.filtered.tidy %>%
  dplyr::filter(!INDIVIDUALS %in% c("PFAR-R3-05-filtered-sorted", "PFAR-R3-26-filtered-sorted"))
write.table(pver.14reef.LE.filtered.noType7, file = "pver.LE.filtered.tidy.296ind.tsv", sep = "\t", row.names = F, quote = F)

# generate strata file for the filtered data set
pver.new.strata <- radiator::generate_strata(data = pver.14reef.LE.filtered.noType7, pop.id = TRUE)
pver.new.strata <- pver.new.strata[,c(2,1)]
write.table(pver.new.strata, file = "pver.new.strata.296ind.NOtype7.tsv", quote = F, row.names = F)

# To close the GDS file
# SNPRelate::snpgdsClose or
closefn.gds(pver.14reef.snprelate)

