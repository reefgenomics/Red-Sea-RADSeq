setwd("~/Documents/Bioinformatics_scripts/R_scripts/Spis_Pver/")

library(phyloseq)
library(ggplot2)
library(reshape2)
library(microbiome)
library(patchwork)

asv=read.table("Input_files/SpisPver_ASVs_noContanoOut.txt", header = TRUE, row.names = 1)[, 1:660]
temp=read.table("Input_files/temperature_data.txt", header = T,  sep = "\t")
temp$temp_category=as.factor(round(temp$Temperature, 0))
temp$temp_category2=ifelse(temp$Temperature >= 27 & temp$Temperature < 28, "cool",
                           ifelse(temp$Temperature >= 29 & temp$Temperature < 30, "med-cool",
                                  ifelse(temp$Temperature >= 30 & temp$Temperature < 31, "med-warm",
                                         ifelse(temp$Temperature >= 31 & temp$Temperature < 32, "warm", "other"))))

map=read.table("Input_files/SpisPver_metadata.txt", header = T, row.names = 1, sep = "\t")
map$site=factor(map$site, levels=c("MAQ","WAJ", "YAN","KAU","DOG" , "FAR"))
map$site_reef=paste(map$site, map$Reef, sep = "_")
host=read.table("Input_files/new_host_clusters_April2021.txt", header = T, sep = "\t")
host$INDIVIDUALS=gsub("-", "_", host$INDIVIDUALS)
map$PopID=host$STRATA[match(rownames(map),host$INDIVIDUALS)]
map$Temperature=temp$temp_category2[match(map$site_reef, temp$Reef)]
tax=read.table("Input_files/SpisPver_ASVs_noContanoOut.txt", header = TRUE, row.names = 1)[, 662:667]
asv_filt=asv[,which(colnames(asv) %in% host$INDIVIDUALS)]

otu.t= otu_table(asv_filt, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

##transform data and subset phyloseq objects
phy.t=microbiome::transform(phy.all, transform = "compositional", target = "OTU", shift = 0, scale = 1)
poc=subset_samples(phy.t, Species=="Pocillopora")
sty=subset_samples(phy.t, Species=="Stylophora")

Spis_Pal=c("#FE92CD", "#1E90FF","#FF4500","#C71585",  "#FFD700", "#32CD32")
Pver_Pal=c("#9370DB",  "#00CED1")
P4=c("#2E33D1", "#FFEE32","#D37D47", "#F43535")

#Temperature ordination
poc_ord = ordinate(poc, method = "PCoA", distance = "bray")
poc_plot_reef=plot_ordination(poc,poc_ord, color = "Temperature")  + geom_point(size = 2, alpha = 1) + theme_bw()  + 
  ggtitle("Pocillopora verrucosa") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4) + theme_bw()
poc_plot_clust=plot_ordination(poc,poc_ord, color = "PopID")  + geom_point(size = 2, alpha = 1) + theme_bw()  + 
  ggtitle("Pocillopora verrucosa") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=Pver_Pal) + theme_bw()

sty_ord = ordinate(sty, method = "PCoA", distance = "bray")
sty_plot_reef=plot_ordination(sty,sty_ord, color = "Temperature")  + geom_point(size = 2, alpha = 1) + theme_bw()  + 
  ggtitle("Stylophora pistillata") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4) + theme_bw()
sty_plot_clust=plot_ordination(sty,sty_ord, color = "PopID")  + geom_point(size = 2, alpha = 1) + theme_bw()  + 
  ggtitle("Stylophora pistillata") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=Spis_Pal) + theme_bw()

pdf("outputs/SpisPver_ordination_popGenClusters_temperature.pdf", width=8,height=7, pointsize = 10)
(poc_plot_reef+sty_plot_reef)/(poc_plot_clust+sty_plot_clust)
dev.off() 
