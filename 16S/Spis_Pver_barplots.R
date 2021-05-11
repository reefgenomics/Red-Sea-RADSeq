
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(gridExtra)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Spis_Pver/")
asv=read.table("Input_files/SpisPver_ASVs_noContanoOut.txt", header = TRUE, row.names = 1)
map=read.table("Input_files/SpisPver_metadata.txt", header = T, row.names = 1, sep = "\t")
list=read.table("Input_files/new_barplots_set.txt")
list$V1=gsub("-", "_", list$V1)
subset(list, !list$V1 %in%  colnames(asv[, 1:660]))

######################################################################
#### Taxonomic profiles of the 10 most abundant genera per sample ####
######################################################################

names(asv)
asv.tax.ag=aggregate(asv[, 1:660], by = list(asv[, 666]), FUN =  sum) #genus
topFamilies=asv.tax.ag[order(rowSums(asv.tax.ag[,2:ncol(asv.tax.ag) ]),decreasing = TRUE),][1:10,1] 
asv.tax.ag$Group.1=ifelse(asv.tax.ag$Group.1 %in% topFamilies, as.character(asv.tax.ag$Group.1), "zOthers")
asv.tax.ag2=aggregate(asv.tax.ag[, 2:ncol(asv.tax.ag)], by = list(asv.tax.ag$Group.1), FUN =  sum)
asv.tax.ag2$SMAQ_R1_18=0
asv.tax.ag2$SMAQ_R1_17=0
asv.tax.ag2$SMAQ_R1_30=0
all.l=reshape2::melt(asv.tax.ag2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Genus","Sample","Abundance")

## Add sample information
all.l$Species=map$Species[match(all.l$Sample, rownames(map))]
final=all.l  %>% group_by(Species,Sample, Genus) %>% dplyr::summarise(Abundance=sum(Abundance))
final$Family=asv$Family[match(final$Genus, asv$Genus)]
final$Tax=paste(final$Genus," (",final$Family, ") ",sep ="") 
final$Reef=paste(map$site,map$Reef,sep ="_")[match(final$Sample, rownames(map))] 
final$Reef=factor(final$Reef, levels=c("MAQ_R1","MAQ_R2","WAJ_R1","WAJ_R2","WAJ_R3","WAJ_R4", "YAN_R1","YAN_R3","YAN_R4","KAU_R1","KAU_R2","KAU_R3","DOG_R1" ,"DOG_R2","DOG_R3", "FAR_R1","FAR_R2","FAR_R3","FAR_R4"))

poci=subset(final, Species == "Pocillopora")
poci$Sample=factor(poci$Sample, levels=subset(list, V1 %like% "^P")$V1)
stylo=subset(final, Species == "Stylophora")
stylo$Sample=factor(stylo$Sample, levels=subset(list, V1 %like% "^S")$V1)

P10=c("#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#ffcc00", "#ADADAD")
#pdf("outputs/SpisPver_ASV_barplots_genus.pdf",  width=7,height=7, pointsize = 12)

legend=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Tax), data = poci, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme_minimal()+ theme(axis.text.x=element_text(angle=90,hjust=1), legend.position = "right",  plot.title = element_text(hjust = 0.5),panel.spacing = unit(2, "lines")) + labs( y= "Percentage of 16S rRNA sequences", x="") + scale_fill_manual(values=P10)  + labs(title="") + guides(fill=guide_legend(ncol=1)) #+ scale_x_discrete(limits = rev(levels(final$Species)))
plot_legend=get_legend(legend, position = NULL)

pver_plot=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Tax ), data = poci, stat="identity", position = "fill", width=1) +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1), legend.position = "none",  plot.title = element_text(hjust = 0.5),panel.spacing = unit(0, "lines")) + 
  labs( y= "Percentage of 16S rRNA sequences", x="") + scale_fill_manual(values=P10)  + 
  labs(title="") + guides(fill=guide_legend(ncol=1)) + 
  facet_grid(.~Reef, scales="free", space="free_x") +scale_x_discrete()#+ scale_x_discrete(limits = rev(levels(final$Sample)))

spis_plot=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Tax), data = stylo, stat="identity", position = "fill", width=1) +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1), legend.position = "none",  plot.title = element_text(hjust = 0.5),panel.spacing = unit(0, "lines")) + 
  labs( y= "Percentage of 16S rRNA sequences", x="") + scale_fill_manual(values=P10)  + 
  labs(title="") + guides(fill=guide_legend(ncol=1)) + 
  facet_grid(.~Reef, scales="free", space="free_x") +scale_x_discrete()#+ scale_x_discrete(limits = rev(levels(final$Sample)))

pdf("outputs/SpisPver_ASV_barplots_genus.pdf",  width=20,height=7, pointsize = 12)
pver_plot/spis_plot | plot_legend
dev.off() 



######################################################################
#####Taxonomic profiles of the 10 most abundant ASV per sample#####
######################################################################
names(asv)
asv$ASV=ifelse(rownames(asv) %in% c("ASV0001", "ASV0002","ASV0003","ASV0004","ASV0005","ASV0006","ASV0007","ASV0008","ASV0009","ASV0010","ASV0011","ASV0012","ASV0013","ASV0014","ASV0015","ASV0016","ASV0017","ASV0018","ASV0019","ASV0020" ), as.character(rownames(asv)), "Others")
all.2=aggregate(asv[, 1:660], by = list(asv$ASV), FUN =  sum)
all.2$SMAQ_R1_18=0
all.2$SMAQ_R1_17=0
all.2$SMAQ_R1_30=0
all.l=reshape2::melt(all.2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("ASV","Sample","Abundance")

## Add sample information
all.l$Species=map$Species[match(all.l$Sample, rownames(map))]
final=all.l  %>% group_by(Species,Sample, ASV) %>% dplyr::summarise(Abundance=sum(Abundance))
#final$Family=asv$Family[match(final$ASV, asv$ASV)]
final$Tax=paste(asv$ASV," (",asv$Genus,"-",asv$Family ,") ",sep ="")[match(final$ASV, asv$ASV)]
final$Reef=paste(map$site,map$Reef,sep ="_")[match(final$Sample, rownames(map))] 
final$Reef=factor(final$Reef, levels=c("MAQ_R1","MAQ_R2","WAJ_R1","WAJ_R2","WAJ_R3","WAJ_R4", "YAN_R1","YAN_R3","YAN_R4","KAU_R1","KAU_R2","KAU_R3","DOG_R1" ,"DOG_R2","DOG_R3", "FAR_R1","FAR_R2","FAR_R3","FAR_R4"))

poci=subset(final, Species == "Pocillopora")
poci$Sample=factor(poci$Sample, levels=subset(list, V1 %like% "^P")$V1)
stylo=subset(final, Species == "Stylophora")
stylo$Sample=factor(stylo$Sample, levels=subset(list, V1 %like% "^S")$V1)

P20=c("#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#ffcc00",  "#2077b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df89","#17becf","#9edae5","#e377c2","#f7b6d2","#ADADAD")
#pdf("outputs/SpisPver_ASV_barplots_genus.pdf",  width=7,height=7, pointsize = 12)

legend=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Tax), data = poci, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme_minimal()+ theme(axis.text.x=element_text(angle=90,hjust=1), legend.position = "right",  plot.title = element_text(hjust = 0.5),panel.spacing = unit(2, "lines")) + labs( y= "Percentage of 16S rRNA sequences", x="") + scale_fill_manual(values=P20)  + labs(title="") + guides(fill=guide_legend(ncol=1)) #+ scale_x_discrete(limits = rev(levels(final$Species)))
plot_legend=get_legend(legend, position = NULL)

pver_plot=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Tax ), data = poci, stat="identity", position = "fill", width=1) +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1), legend.position = "none",  plot.title = element_text(hjust = 0.5),panel.spacing = unit(0, "lines")) + 
  labs( y= "Percentage of 16S rRNA sequences", x="") + scale_fill_manual(values=P20)  + 
  labs(title="") + guides(fill=guide_legend(ncol=1)) + 
  facet_grid(.~Reef, scales="free", space="free_x") +scale_x_discrete()#+ scale_x_discrete(limits = rev(levels(final$Sample)))

spis_plot=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Tax), data = stylo, stat="identity", position = "fill", width=1) +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1), legend.position = "none",  plot.title = element_text(hjust = 0.5),panel.spacing = unit(0, "lines")) + 
  labs( y= "Percentage of 16S rRNA sequences", x="") + scale_fill_manual(values=P20)  + 
  labs(title="") + guides(fill=guide_legend(ncol=1)) + 
  facet_grid(.~Reef, scales="free", space="free_x") +scale_x_discrete()#+ scale_x_discrete(limits = rev(levels(final$Sample)))

pdf("outputs/SpisPver_ASV_barplots_ASV.pdf",  width=20,height=7, pointsize = 12)
pver_plot/spis_plot | plot_legend
dev.off() 


### top genera per species ###
Pver=asv.tax.ag[,grepl("^P|^G", colnames(asv.tax.ag))]
Spis=asv.tax.ag[,grepl("^S|^G", colnames(asv.tax.ag))]

Pver_topGenera=Pver[order(rowSums(Pver[,2:ncol(Pver) ]),decreasing = TRUE),][1:10,1] 
Spis_topGenera=Spis[order(rowSums(Spis[,2:ncol(Spis) ]),decreasing = TRUE),][1:10,1] 
intersect(Pver_topGenera,Spis_topGenera)


