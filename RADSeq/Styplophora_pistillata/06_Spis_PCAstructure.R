###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 06 PCR population structure - Stylophora pistillata #####

#### 06.01 Evaluate potential patterns of population structure using a multivariate free approach as PCA

# Set working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_PCAplink_20210222")

# Run PCA
# variants aren't ordered in ascending order in each chromosome, so bed file has to be done first to overcome the problem and then PCA can be perform (plink1.9)
system("~/RADseq/Genotyping/tools/plink --vcf ../spis_vcftools_filtered_vcf/spis.LE.filtered.recode.indnames.vcf --allow-extra-chr --make-bed --out spis.LE.filtered.4PCA")
system("~/RADseq/Genotyping/tools/plink --bfile spis.LE.filtered.4PCA --pca --allow-extra-chr --out spis.LE.filtered.PCA")


#### 06.02 PCA visualization (LE dataset)

# Libraries
library(tidyverse)
library(scales)
library(PCAviz)
library(ggplot2)
library(magrittr)

# Read PCA output from PLINK 1.9
spis.pca <- read_table2("./spis.LE.filtered.PCA.eigenvec", col_names = FALSE)
spis.eigenval <- scan("./spis.LE.filtered.PCA.eigenval")

# sort out the pca data
# remove nuisance column
spis.pca <- spis.pca[,-2]
# set names
names(spis.pca)[1] <- "ind"
names(spis.pca)[2:ncol(spis.pca)] <- paste0("PC", 1:(ncol(spis.pca)-1))

# sort out the individual species and pops
# reef
spis.pca$reef <- sub('^([^-]+-[^-]+).*', '\\1', spis.pca$ind)
spis.pca$reef <- gsub("S","", spis.pca$reef)
spis.pca$reef <- factor(spis.pca$reef, levels = c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R3", "WAJ-R4", "YAN-R1", "YAN-R3", "YAN-R4", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R1", "FAR-R2", "FAR-R3", "FAR-R4"), ordered = T)
spis.pca$region <- gsub("-.*","", spis.pca$reef)
spis.pca$region <- factor(spis.pca$region, levels = c("MAQ", "WAJ", "YAN", "KAU", "DOG", "FAR"))

# first convert to percentage variance explained
spis.pve <- data.frame(PC = 1:20, pve = spis.eigenval/sum(spis.eigenval)*100)
spis.sdev <- sqrt(spis.eigenval)
  
#create a bar plot showing the percentage of variance each principal component explains.
ggplot(spis.pve, aes(x=PC, y=pve)) + geom_bar(stat = "identity") + ylab("Percentage variance explained") + theme_light()
ggsave("spis.pca.variance.percentage.explained.pdf", width = 7, height = 6)

# calculate the cumulative sum of the percentage variance explained
cumsum(spis.pve$pve) # the first 3 PCA explain 34.36655% of the cumulative variance

# Color for each reef
spis.reef.color <- c("#0066CC","#3399FF","#009999", "#00CCCC","#33FFFF","#009900","#66CC00","#80FF00","#CCCC00","#FFFF00","#FFFF66","#D2691E","#FF8000","#FFB266","#CC0000","#FF3333","#FF6666","#FF9999" )
show_col(spis.reef.color)

# Visualization using PCAviz
# pca matrix
spis.pcavis <- spis.pca[,c(-1,-22,-23)]

# Create an informative matrix with latitudinal data
spis.pop.latitude.coord <- c(rep(28.52616667, 22), rep(28.4268, 22), rep(26.18751389, 12), rep(26.24119444, 23), rep(26.18505556, 25), rep(23.94741667, 19), rep(23.95533333, 25), rep(23.9115, 25), rep(22.31916667, 8), rep(22.06722222, 18), rep(22.51333333, 20), rep(19.63511111, 14), rep(19.61402778, 21), rep(19.66569167, 25), rep(16.57930556, 25), rep(16.57899444, 19), rep(16.52518056, 18), rep(16.52736111, 26) ) 

spis.info.mat <- data.frame(Reef=spis.pca$reef, Individuals=spis.pca$ind, Region=spis.pca$region, Latitude=spis.pop.latitude.coord) 

# Create a PCAviz class object
source("pcaviz.shape.R") # To have 6 differnt shapes for each region
spis.pcaviz <- pcaviz(x = spis.pcavis, dat = spis.info.mat, sdev = spis.sdev, var= 114.8153) 

# Custom-made theme for the PCAviz plot
My_Theme = theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11)) 



# Plot Discriminant function 1 against latitude with group centroids (summary of information)
pdf("Spis_PC1_vs_Latitude.pdf")
plot.pcaviz(spis.pcaviz, coords = c("PC1","Latitude"), group="Reef", 
     include.with.pc.axes= "pve", 
     show.legend = T,color = "Reef", shape="Region", colors= spis.reef.color,
     draw.points =T, group.summary.labels = F, theme = My_Theme,
     geom.point.params= list(size = 5, alpha = 0.7), 
     geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("Spis_PC2_vs_Latitude.pdf")
plot.pcaviz(spis.pcaviz, coords = c("PC2","Latitude"), group="Reef", 
     include.with.pc.axes= "pve", 
     show.legend = T,color = "Reef", shape="Region", colors= spis.reef.color,
     draw.points =T, group.summary.labels = F, theme = My_Theme,
     geom.point.params= list(size = 5, alpha = 0.7), 
     geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("Spis_PC3_vs_Latitude.pdf")
plot.pcaviz(spis.pcaviz, coords = c("PC3","Latitude"), group="Reef", 
     include.with.pc.axes= "pve", 
     show.legend = T,color = "Reef", shape="Region", colors= spis.reef.color,
     draw.points =T, group.summary.labels = F, theme = My_Theme,
     geom.point.params= list(size = 5, alpha = 0.7), 
     geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("Spis_PC1_vs_PC2.pdf")
plot.pcaviz(spis.pcaviz, coords = c("PC1","PC2"), group="Reef", 
     include.with.pc.axes= "pve", 
     show.legend = T,color = "Reef", shape="Region", colors= spis.reef.color,
     draw.points =T, group.summary.labels = F, theme = My_Theme,
     geom.point.params= list(size = 5, alpha = 0.7), 
     geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("Spis_PC1_vs_PC3.pdf")
plot.pcaviz(spis.pcaviz, coords = c("PC1","PC3"), group="Reef", 
            include.with.pc.axes= "pve", 
            show.legend = T,color = "Reef", shape="Region", colors= spis.reef.color,
            draw.points =T, group.summary.labels = F, theme = My_Theme,
            geom.point.params= list(size = 5, alpha = 0.7), 
            geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("Spis_PC2_vs_PC3.pdf")
plot.pcaviz(spis.pcaviz, coords = c("PC2","PC3"), group="Reef", 
            include.with.pc.axes= "pve", 
            show.legend = T,color = "Reef", shape="Region", colors= spis.reef.color,
            draw.points =T, group.summary.labels = F, theme = My_Theme,
            geom.point.params= list(size = 5, alpha = 0.7), 
            geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

