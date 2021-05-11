#library(GUniFrac)
library(vegan)

#data in
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Spis_Pver/")
asv=read.table("Input_files/SpisPver_ASVs_noContanoOut.txt", header = TRUE, row.names = 1)[, 1:660]
map=read.table("Input_files/SpisPver_metadata.txt", header = T, row.names = 1, sep = "\t")
temp=read.table("Input_files/temperature_data.txt", header = T,  sep = "\t")
temp$temp_category=as.factor(round(temp$Temperature, 0))
temp$temp_category2=ifelse(temp$Temperature >= 27 & temp$Temperature < 28, "cool",
                           ifelse(temp$Temperature >= 29 & temp$Temperature < 30, "med-cool",
                                  ifelse(temp$Temperature >= 30 & temp$Temperature < 31, "med-warm",
                                         ifelse(temp$Temperature >= 31 & temp$Temperature < 32, "warm", "other"))))

host=read.table("Input_files/new_host_clusters_April2021.txt", header = T, sep = "\t")
host$INDIVIDUALS=gsub("-", "_", host$INDIVIDUALS)
asv_filt2=asv[,which(colnames(asv) %in% host$INDIVIDUALS)]

###rarefying
# cnts=t(asv[,1:asv_number])
# min(rowSums(cnts)) # determine sample with lowest counts
# asv.rar=Rarefy(cnts, 1018)$otu.tab.rff

###alpha diversity stimates
alpha=as.data.frame(t(estimateR(t(asv_filt2),  smallsample = TRUE)))
alpha$Shannon=vegan::diversity(t(asv_filt2), index = "shannon")#$shannon
alpha$Reef=paste(map$site,map$Reef, sep = "_")[match(rownames(alpha),rownames(map))]
alpha$Reef=factor(alpha$Reef, levels=c("MAQ_R1","MAQ_R2","WAJ_R1","WAJ_R2","WAJ_R3","WAJ_R4", "YAN_R1","YAN_R3","YAN_R4","KAU_R1","KAU_R2","KAU_R3","DOG_R1" ,"DOG_R2","DOG_R3", "FAR_R1","FAR_R2","FAR_R3","FAR_R4"))
alpha$Species=ifelse(rownames(alpha) %like% "^S", "Stylophora", "Pocillopora")
alpha$PopID=host$STRATA[match(rownames(alpha),host$INDIVIDUALS)]
alpha$Temperature=temp$temp_category2[match(alpha$Reef,temp$Reef)]

##text observed
#species_pool=ifelse(rownames(t(asv_filt2)) %like% "^S", "Stylophora", "Pocillopora")
#specpool(t(asv_filt2), species_pool) # Pocillopora 15,577 and Stylophora 23,736

# colors

Spis_Pal=c("#FE92CD", "#1E90FF","#FF4500","#C71585",  "#FFD700", "#32CD32")
Pver_Pal=c("#9370DB",  "#00CED1")
temp_Pal=c("#2E33D1", "#FFEE32","#D37D47", "#F43535")


P6=c("#222E50", "#007991", "#BCD8C1", "#E9D985", "#F29469", "#BE3F23") 
P4=c("#2E33D1", "#FFEE32","#D37D47", "#F43535")
P20=c("#fad390", "#f6b93b", "#fa983a", "#e58e26", "#f8c291", "#e55039", "#eb2f06", "#b71540", "#6a89cc", "#4a69bd","#1e3799", "#0c2461", "#82ccdd", "#60a3bc", "#3c6382", "#0a3d62", "#b8e994", "#78e08f", "#38ada9", "#079992", "#C0C0C0")

##################################################
##################### Stats ######################
##################################################

## 1 between species
shapiro.test(alpha$Shannon) # p-value < 0.05 implying we can't assume the normality.
kruskal.test(alpha$Shannon ~ alpha$Species) # Kruskal-Wallis chi-squared = 21.373, df = 1, p-value = 3.781e-06
pdf("outputs/SpisPver_alphaDiversity_species.pdf", width=4,height=4, pointsize = 10)
ggplot(alpha, aes(x=Species, y=Shannon, fill=Species)) + 
  stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) + 
  scale_fill_manual(values=P6) +  theme_classic() + 
  labs( y= "ASV Shannon diversity", x="", title = "") +
  annotate(geom="text", x=1.5, y=5.5, label= "Kruskal-Wallis \n p-value = 3.8e-06")
dev.off()

# 2. Pver per reef and cluster
pver_alpha=subset(alpha, Species == "Pocillopora")

#cluster
kruskal.test(pver_alpha$Shannon ~ pver_alpha$PopID) # Kruskal-Wallis chi-squared = 9.1621, df = 1, p-value = 0.002471
wilcox_alphaD_pver_popid=pairwise.wilcox.test(pver_alpha$Shannon, pver_alpha$PopID, p.adj = "fdr")
#write.table(wilcox_alphaD_pver_popid$p.value,"outputs/pver_wilcox_alphaD_popID.txt", quote = F, sep = "\t")
plot_pver_alpha_popID=ggplot(pver_alpha, aes(x=PopID, y=Shannon, fill=PopID)) + 
  stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) + 
  scale_fill_manual(values=Pver_Pal) +  theme_classic() + 
  labs( y= "ASV Shannon diversity", x="", title = "") +
  annotate(geom="text", x=1, y=5.5, label= "A") + annotate(geom="text", x=2, y=5.5, label= "B") + 
  guides(fill=guide_legend(title="Host genetic\ncluster"))

#reef
kruskal.test(pver_alpha$Shannon ~ pver_alpha$Reef) # Kruskal-Wallis chi-squared = 9.1621, df = 1, p-value = 0.002471
wilcox_alphaD_pver_reef=pairwise.wilcox.test(pver_alpha$Shannon, pver_alpha$Reef, p.adj = "fdr")
#write.table(wilcox_alphaD_pver_reef$p.value,"outputs/pver_wilcox_alphaD_reef.txt", quote = F, sep = "\t")
plot_pver_alpha_reef=ggplot(pver_alpha, aes(x=Reef, y=Shannon, fill=Reef)) + 
  stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) + 
  scale_fill_manual(values=P20) +  theme_classic() + 
  labs( y= "ASV Shannon diversity", x="", title = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#temperature
kruskal.test(pver_alpha$Shannon ~ pver_alpha$Temperature) # Kruskal-Wallis chi-squared = 9.1621, df = 1, p-value = 0.002471
wilcox_alphaD_pver_reef=pairwise.wilcox.test(pver_alpha$Shannon, pver_alpha$Temperature, p.adj = "fdr")
#write.table(wilcox_alphaD_pver_reef$p.value,"outputs/pver_wilcox_alphaD_temp.txt", quote = F, sep = "\t")
plot_pver_alpha_temp=ggplot(pver_alpha, aes(x=Reef, y=Shannon, fill=Temperature)) + 
  stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) + 
  scale_fill_manual(values=temp_Pal) +  theme_classic() + 
  labs( y= "ASV Shannon diversity", x="", title = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  annotate("segment", x = 1, xend = 2, y = 5.5, yend = 5.5, colour="#2E33D1") + annotate(geom="text", x=1.5, y=5.8, label= "A") +
  annotate("segment", x = 3, xend = 7, y = 5.5, yend = 5.5, colour="#FFEE32") + annotate(geom="text", x=5, y=5.8, label= "B") + 
  annotate("segment", x = 8, xend = 10, y = 5.5, yend = 5.5, colour="#D37D47") + annotate(geom="text", x=9, y=5.8, label= "C") +
  annotate("segment", x = 11, xend = 14, y = 5.5, yend = 5.5, colour="#F43535") + annotate(geom="text", x=12.5, y=5.8, label= "C") 
#pdf("outputs/Pver_alphaDiversity.pdf", width=10,height=4, pointsize = 10)
#plot_pver_alpha_popID+plot_pver_alpha_reef
#dev.off()

# 3. Spis per reef and cluster
spis_alpha=subset(alpha, Species == "Stylophora")

#cluster
kruskal.test(spis_alpha$Shannon ~ spis_alpha$PopID) # Kruskal-Wallis chi-squared = 13.645, df = 5, p-value = 0.01803
wilcox_alphaD_spis_popid=pairwise.wilcox.test(spis_alpha$Shannon, spis_alpha$PopID, p.adj = "fdr")
#write.table(wilcox_alphaD_spis_popid$p.value,"outputs/spis_wilcox_alphaD_popID.txt", quote = F, sep = "\t")
plot_spis_alpha_popID=ggplot(spis_alpha, aes(x=PopID, y=Shannon, fill=PopID)) + 
  stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) + 
  scale_fill_manual(values=Spis_Pal) +  theme_classic() + 
  labs( y= "ASV Shannon diversity", x="", title = "") + 
  annotate(geom="text", x=1, y=5.6, label= "A") +
  annotate(geom="text", x=2, y=5.6, label= "B") +
  annotate(geom="text", x=3, y=5.6, label= "AB") +
  annotate(geom="text", x=4, y=5.6, label= "A") +
  annotate(geom="text", x=5, y=5.6, label= "A") +
  annotate(geom="text", x=6, y=5.6, label= "A")  + 
  guides(fill=guide_legend(title="Host genetic\ncluster"))
  

#reef
kruskal.test(spis_alpha$Shannon ~ spis_alpha$Reef) # Kruskal-Wallis chi-squared = 105.06, df = 17, p-value = 1.017e-14
wilcox_alphaD_spis_reef=pairwise.wilcox.test(spis_alpha$Shannon, spis_alpha$Reef, p.adj = "fdr")
#write.table(wilcox_alphaD_spis_reef$p.value,"outputs/spis_wilcox_alphaD_reef.txt", quote = F, sep = "\t")
plot_spis_alpha_reef=ggplot(spis_alpha, aes(x=Reef, y=Shannon, fill=Reef)) + 
  stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) + 
  scale_fill_manual(values=P20) +  theme_classic() + 
  labs( y= "ASV Shannon diversity", x="", title = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#temperature
kruskal.test(spis_alpha$Shannon ~ spis_alpha$Temperature) # Kruskal-Wallis chi-squared = 105.06, df = 17, p-value = 1.017e-14
wilcox_alphaD_spis_reef=pairwise.wilcox.test(spis_alpha$Shannon, spis_alpha$Temperature, p.adj = "fdr")
#write.table(wilcox_alphaD_spis_reef$p.value,"outputs/spis_wilcox_alphaD_temp.txt", quote = F, sep = "\t")
plot_spis_alpha_temp=ggplot(spis_alpha, aes(x=Reef, y=Shannon, fill=Temperature)) + 
  stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) + 
  scale_fill_manual(values=temp_Pal) +  theme_classic() + 
  labs( y= "ASV Shannon diversity", x="", title = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  annotate("segment", x = 1, xend = 2, y = 5.5, yend = 5.5, colour="#2E33D1") + annotate(geom="text", x=1.5, y=5.8, label= "A") +
  annotate("segment", x = 3, xend = 8, y = 5.5, yend = 5.5, colour="#FFEE32") + annotate(geom="text", x=5.5, y=5.8, label= "B") + 
  annotate("segment", x = 9, xend = 11, y = 5.5, yend = 5.5, colour="#D37D47") + annotate(geom="text", x=10, y=5.8, label= "B") +
  annotate("segment", x = 12, xend = 18, y = 5.5, yend = 5.5, colour="#F43535") + annotate(geom="text", x=15, y=5.8, label= "B")


#pdf("outputs/spis_alphaDiversity.pdf", width=10,height=4, pointsize = 10)
#plot_spis_alpha_popID+plot_spis_alpha_reef
#dev.off()

pdf("outputs/alphaDiversity.pdf", width=10,height=8, pointsize = 10)
(plot_pver_alpha_popID+plot_pver_alpha_temp)/(plot_spis_alpha_popID+plot_spis_alpha_temp) + plot_annotation(tag_levels = 'A') 
dev.off()

