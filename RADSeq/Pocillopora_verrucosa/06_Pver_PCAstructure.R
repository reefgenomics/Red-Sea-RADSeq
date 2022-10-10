###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 06 PCR population structure - Pocillopora verrucosa #####

#### 06.01 Evaluate potential patterns of population structure using a multivariate free approach as PCA


# Set working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_PCAplink_20210225")


# Run PCA
# variants aren't ordered in ascending order in each chromosome, so bed file has to be done first to overcome the problem and then PCA can be perform (plink1.9)
system("~/RADseq/Genotyping/tools/plink --vcf ../pver_vcftools_filtered_vcf/pver.LE.filtered.recode.vcf --allow-extra-chr --make-bed --out pver.LE.filtered.PCA")
system("~/RADseq/Genotyping/tools/plink --bfile pver.LE.filtered.PCA --pca --allow-extra-chr --out pver.LE.filtered.PCA")

#### 06.02 PCA visualization (LE dataset)
# PCA visualization of Stylophora istillata (LE dataset)

# Libraries
library(tidyverse)
library(scales)
library(PCAviz)
library(ggplot2)
library(magrittr)

# Set working directory
setwd("~/RADseq-Big-project/pver/p1.mac4.r0.8.316samples.Pver_20200218/05_PCAplink_20210225")

# Read PCA output from PLINK 1.9
pver.pca <- read_table2("./pver.LE.filtered.PCA.eigenvec", col_names = FALSE)
pver.eigenval <- scan("./pver.LE.filtered.PCA.eigenval")

# sort out the pca data
# remove nuisance column
pver.pca <- pver.pca[,-2]
# set names
names(pver.pca)[1] <- "ind"
names(pver.pca)[2:ncol(pver.pca)] <- paste0("PC", 1:(ncol(pver.pca)-1))

# sort out the individual species and pops
# reef
pver.pca$reef <- sub('^([^-]+-[^-]+).*', '\\1', pver.pca$ind)
pver.pca$reef <- gsub("P","", pver.pca$reef)
pver.pca$reef <- factor(pver.pca$reef, levels = c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R2", "WAJ-R3", "YAN-R1", "YAN-R3", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R3"), ordered = T)
pver.pca$region <- gsub("-.*","", pver.pca$reef)
pver.pca$region <- factor(pver.pca$region, levels = c("MAQ", "WAJ", "YAN", "KAU", "DOG", "FAR"))

# first convert to percentage variance explained
total.var <- sum(pver.eigenval)
pver.pve <- data.frame(PC = 1:20, pve = pver.eigenval/sum(pver.eigenval)*100)
pver.sdev <- sqrt(pver.eigenval)
  
#create a bar plot showing the percentage of variance each principal component explains.
ggplot(pver.pve, aes(x=PC, y=pve)) + geom_bar(stat = "identity") + ylab("Percentage variance explained") + theme_light()
ggsave("Pver.pca.variance.percentage.explained.pdf", width = 7, height = 6)

# calculate the cumulative sum of the percentage variance explained
cumsum(pver.pve$pve) # the first 3 PCA explain 17.395958% of the cumulative variance

# Color for each reef
pver.reef.color <- c("#0066CC","#3399FF","#009999", "#27C3AE", "#00CCCC", "#009900","#66CC00", "#CCCC00","#FFFF00","#FFFF66","#D2691E","#FF8000","#FFB266","#FF6666")
show_col(pver.reef.color)

# Visualization using PCAviz
# pca matrix
pver.pcavis <- pver.pca[,c(-1,-22,-23)]

# Create an informative matrix with latitudinal data
pver.pop.latitude.coord <- c(rep(28.52616667, 12), rep(28.4268, 23), rep(26.18751389, 25), rep(26.16661111, 25), rep(26.24119444, 20), rep(23.94741667, 17), rep(23.95533333, 17), rep(22.31916667, 25), rep(22.06722222, 23), rep(22.51333333, 25), rep(19.63511111, 22), rep(19.61402778, 20), rep(19.66569167, 20), rep(16.52518056, 22))

pver.info.mat <- data.frame(Reef=pver.pca$reef, Individuals=pver.pca$ind, Region=pver.pca$region, Latitude=pver.pop.latitude.coord) 

source("pcaviz.shape.R") # To have 6 differnt shapes for each region

# Create a PCAviz class object
pver.pcaviz <- pcaviz(x = pver.pcavis, dat = pver.info.mat, sdev = pver.sdev, var= total.var) 

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
pdf("pver_PC1_vs_Latitude.pdf")
plot.pcaviz(pver.pcaviz, coords = c("PC1","Latitude"), group="Reef", 
     include.with.pc.axes= "pve", 
     show.legend = T,color = "Reef", shape="Region", colors= pver.reef.color,
     draw.points =T, group.summary.labels = F, theme = My_Theme,
     geom.point.params= list(size = 5, alpha = 0.7), 
     geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("pver_PC2_vs_Latitude.pdf")
plot.pcaviz(pver.pcaviz, coords = c("PC2","Latitude"), group="Reef", 
     include.with.pc.axes= "pve", 
     show.legend = T,color = "Reef", shape="Region", colors= pver.reef.color,
     draw.points =T, group.summary.labels = F, theme = My_Theme,
     geom.point.params= list(size = 5, alpha = 0.7), 
     geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("pver_PC3_vs_Latitude.pdf")
plot.pcaviz(pver.pcaviz, coords = c("PC3","Latitude"), group="Reef", 
     include.with.pc.axes= "pve", 
     show.legend = T,color = "Reef", shape="Region", colors= pver.reef.color,
     draw.points =T, group.summary.labels = F, theme = My_Theme,
     geom.point.params= list(size = 5, alpha = 0.7), 
     geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("pver_PC1_vs_PC2.pdf")
plot.pcaviz(pver.pcaviz, coords = c("PC1","PC2"), group="Reef", 
     include.with.pc.axes= "pve", 
     show.legend = T,color = "Reef", shape="Region", colors= pver.reef.color,
     draw.points =T, group.summary.labels = F, theme = My_Theme,
     geom.point.params= list(size = 5, alpha = 0.7), 
     geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("pver_PC1_vs_PC3.pdf")
plot.pcaviz(pver.pcaviz, coords = c("PC1","PC3"), group="Reef", 
            include.with.pc.axes= "pve", 
            show.legend = T,color = "Reef", shape="Region", colors= pver.reef.color,
            draw.points =T, group.summary.labels = F, theme = My_Theme,
            geom.point.params= list(size = 5, alpha = 0.7), 
            geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

pdf("pver_PC2_vs_PC3.pdf")
plot.pcaviz(pver.pcaviz, coords = c("PC2","PC3"), group="Reef", 
            include.with.pc.axes= "pve", 
            show.legend = T,color = "Reef", shape="Region", colors= pver.reef.color,
            draw.points =T, group.summary.labels = F, theme = My_Theme,
            geom.point.params= list(size = 5, alpha = 0.7), 
            geom.point.summary.params = list(size =0,alpha = 0), draw.linear.fit = F)
dev.off()

