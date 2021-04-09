setwd("~/Documents/Bioinformatics_scripts/R_scripts/Spis_Pver/")

library(phyloseq)
library(ggplot2)
library(reshape2)
library(microbiome)
library(patchwork)

asv=read.table("Input_files/SpisPver_ASVs_noContanoOut.txt", header = TRUE, row.names = 1)[, 1:660]
map=read.table("Input_files/SpisPver_metadata.txt", header = T, row.names = 1, sep = "\t")
map$site=factor(map$site, levels=c("MAQ","WAJ", "YAN","KAU","DOG" , "FAR"))
map$site_reef=paste(map$site, map$Reef, sep = "_")
# host=read.table("Input_files/Host_clusters.txt", header = T, sep = "\t")
host=read.table("Input_files/new_host_clusters_April2021.txt", header = T, sep = "\t")
host$INDIVIDUALS=gsub("-", "_", host$INDIVIDUALS)
map$PopID=host$STRATA[match(rownames(map),host$INDIVIDUALS)]
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

P6=c("#222E50", "#007991", "#BCD8C1", "#E9D985", "#F29469", "#BE3F23") #27, 29, 32, 34
P10=c("#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#ffcc00", "#ADADAD")
P10=c("#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#ffcc00",  "#2077b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df89","#17becf","#9edae5","#e377c2","#f7b6d2","#ADADAD")
P20=c("#fad390", "#f6b93b", "#fa983a", "#e58e26", "#f8c291", "#e55039", "#eb2f06", "#b71540", "#6a89cc", "#4a69bd","#1e3799", "#0c2461", "#82ccdd", "#60a3bc", "#3c6382", "#0a3d62", "#b8e994", "#78e08f", "#38ada9", "#079992", "#C0C0C0")

poc_ord = ordinate(poc, method = "PCoA", distance = "bray")
poc_plot_reef=plot_ordination(poc,poc_ord, color = "site_reef")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Pocillopora") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P20) + theme_bw()
poc_plot_clust=plot_ordination(poc,poc_ord, color = "PopID")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Pocillopora") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P6) + theme_bw()

sty_ord = ordinate(sty, method = "PCoA", distance = "bray")
sty_plot_reef=plot_ordination(sty,sty_ord, color = "site_reef")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Stylophora") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P20) + theme_bw()
sty_plot_clust=plot_ordination(sty,sty_ord, color = "PopID")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("Stylophora") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P6) + theme_bw()

poc_plot_reef+sty_plot_reef
poc_plot_clust+sty_plot_clust
  

pdf("outputs/SpisPver_ordination_popGenClusters.pdf", width=10,height=5, pointsize = 10)
gridExtra::grid.arrange( poc_plot, sty_plot, ncol=2)
dev.off()
