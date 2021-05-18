setwd("~/Documents/Bioinformatics_scripts/R_scripts/Spis_Pver/")

library(indicspecies)
library(UpSetR)
library(ComplexHeatmap)
library(patchwork)

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
asv_filt=asv[,which(colnames(asv) %in% host$INDIVIDUALS)]
asv_filt2=asv_filt[rowSums( asv_filt > 0 ) >= 5,]  # at least 5 samples -> 3455

asv.n=as.data.frame(t(apply(asv_filt2, 2,function(x) {x/sum(x)})))
sample_number=ncol(asv.n)

asv.n$PopID=paste(host$STRATA)[match(rownames(asv.n), host$INDIVIDUALS)]
asv.n$Species=ifelse(rownames(asv.n) %like% "^S", "stylphora", "Pocillopora")
asv.n$Reef=paste(map$site, map$Reef, sep = "_")[match(rownames(asv.n), rownames(map))] 
asv.n$Temperature=temp$temp_category2[match(asv.n$Reef, temp$Reef)]

poci=subset(asv.n, Species == "Pocillopora")
styl=subset(asv.n, Species == "stylphora")

######################################################
### Pocillopora : indicator ASVs per temperatures ###
######################################################

poci_indi_temp=multipatt(poci[,1:sample_number],poci$Temperature)
summary(poci_indi_temp)
poci_indi_temp_df=as.data.frame(poci_indi_temp[["sign"]])
poci_indi_temp_df_sig=subset(poci_indi_temp_df, p.value < 0.05 | p.value == 0.05 )
colnames(poci_indi_temp_df_sig)[1:4]=c("cool", "med_cool", "med_warm", "warm")

# make combination matrix
# poci_indi_temp_lt=list(cool=rownames(subset(poci_indi_temp_df_sig, cool == 1)), 
#      med_cool=rownames(subset(poci_indi_temp_df_sig, med_cool == 1)),
#      med_warm=rownames(subset(poci_indi_temp_df_sig,med_warm == 1)),
#      warm=rownames(subset(poci_indi_temp_df_sig,warm == 1)))
# m1=make_comb_mat(poci_indi_temp_lt, mode = "intersect")
# upset1=UpSet(m1, pt_size = unit(4, "mm"), lwd = 2)

######################################################
### stylphora: indicator ASVs per temperatures ###
######################################################

styl_indi_temp=multipatt(styl[,1:sample_number],styl$Temperature)
styl_indi_temp_df=as.data.frame(styl_indi_temp[["sign"]])
styl_indi_temp_df_sig=subset(styl_indi_temp_df, p.value < 0.05 | p.value == 0.05 )
colnames(styl_indi_temp_df_sig)[1:4]=c("cool", "med_cool", "med_warm", "warm")

# # make combination matrix
# styl_indi_temp_lt=list(cool=rownames(subset(styl_indi_temp_df_sig, cool == 1)), 
#                        med_cool=rownames(subset(styl_indi_temp_df_sig, med_cool == 1)),
#                        med_warm=rownames(subset(styl_indi_temp_df_sig,med_warm == 1)),
#                        warm=rownames(subset(styl_indi_temp_df_sig,warm == 1)))
# m2=make_comb_mat(styl_indi_temp_lt, mode = "intersect")
# upset2=UpSet(m2,  pt_size = unit(4, "mm"), lwd = 2)

##################################################
### Pocillopora : indicator ASVs per clusters ###
#################################################

poci_indi_pop=multipatt(poci[,1:sample_number],poci$PopID)
summary(poci_indi_pop)
poci_indi_pop_df=as.data.frame(poci_indi_pop[["sign"]])
poci_indi_pop_df_sig=subset(poci_indi_pop_df, p.value < 0.05 | p.value == 0.05 )
colnames(poci_indi_pop_df_sig)[1:2]=c("CL1", "CL2")

# # make combination matrix
# poci_indi_pop_lt=list(CL1=rownames(subset(poci_indi_pop_df_sig, CL1 == 1)), 
#                        CL2=rownames(subset(poci_indi_pop_df_sig, CL2 == 1)))
# m3=make_comb_mat(poci_indi_pop_lt, mode = "intersect")
# upset3=UpSet(m3,  pt_size = unit(4, "mm"), lwd = 2)

######################################################
### stylphora: indicator ASVs per temperatures ###
######################################################

styl_indi_pop=multipatt(styl[,1:sample_number],styl$PopID)
summary(styl_indi_pop)
styl_indi_pop_df=as.data.frame(styl_indi_pop[["sign"]])
styl_indi_pop_df_sig=subset(styl_indi_pop_df, p.value < 0.05 | p.value == 0.05 )
colnames(styl_indi_pop_df_sig)[1:6]=c("CL1", "CL2","CL3", "CL4","CL5", "CL6")

# # make combination matrix
# styl_indi_pop_lt=list(CL1=rownames(subset(styl_indi_pop_df_sig, CL1 == 1)), 
#                       CL2=rownames(subset(styl_indi_pop_df_sig, CL2 == 1)),
#                       CL3=rownames(subset(styl_indi_pop_df_sig, CL3 == 1)),
#                       CL4=rownames(subset(styl_indi_pop_df_sig, CL4 == 1)),
#                       CL5=rownames(subset(styl_indi_pop_df_sig, CL5 == 1)),
#                       CL6=rownames(subset(styl_indi_pop_df_sig, CL6 == 1)))
# m4=make_comb_mat(styl_indi_pop_lt, mode = "intersect")
# upset4=UpSet(m4,  pt_size = unit(4, "mm"), lwd = 2)


### Plots
Spis_Pal=c("#FE92CD", "#1E90FF","#FF4500","#C71585",  "#FFD700", "#32CD32")
Pver_Pal=c("#9370DB",  "#00CED1")
temp_Pal=c("#2E33D1", "#FFEE32","#D37D47", "#F43535")

upset1=upset(poci_indi_temp_df_sig, 
             sets = c("cool", "med_cool", "med_warm", "warm"), 
             order.by="freq",  point.size=4,
             sets.bar.color=temp_Pal)

upset2=upset(styl_indi_temp_df_sig, 
             sets = c("cool", "med_cool", "med_warm", "warm"), 
             order.by="freq",  point.size=4,
             sets.bar.color=temp_Pal)

upset3=upset(poci_indi_pop_df_sig, 
             sets = c("CL1", "CL2"), 
             order.by="freq",  point.size=3,
             sets.bar.color=Pver_Pal)

upset4=upset(styl_indi_pop_df_sig, 
             sets = c("CL1", "CL2","CL3", "CL4","CL5", "CL6"), 
             order.by="freq",  point.size=3,
             sets.bar.color=Spis_Pal)

pdf("outputs/upset_indic_temp_pver.pdf", height = 4, width = 4)
upset1 
dev.off()

pdf("outputs/upset_indic_temp_spis.pdf",  height = 4, width = 4)
 upset2 
dev.off()


pdf("outputs/upset_indic_pver_clust.pdf",  height = 4, width = 2)
upset3
dev.off()

pdf("outputs/upset_indic_spis_clust.pdf",  height = 4, width = 6)
upset4
dev.off()

# info text
nrow(subset(styl_indi_temp_df_sig,warm == 1 & cool == 0 & med_warm == 0 & med_cool == 0))
nrow(subset(poci_indi_temp_df_sig,warm == 1 & cool == 0 & med_warm == 0 & med_cool == 0))

## donuts

tax=read.table("Input_files/SpisPver_ASVs_noContanoOut.txt", header = TRUE, row.names = 1)[, 661:667]
pastel_pal=c("#81b4e1","#7bd0e2","#a8e6cf","#bedc92","#fff8aa","#fdd389","#fabc97","#f5a3b3","#e194c2","#a486bc")
palA=c( "#969A97","#81b4e1","#7bd0e2","#a8e6cf","#bedc92","#fff8aa","#D1CBC7")

#pver warm
pver_warm=subset(poci_indi_temp_df_sig,index == 4)
pver_warm$Taxa=tax$Family[match(rownames(pver_warm), rownames(tax))]
pver_warm[is.na(pver_warm)] <- "Unclassified"
pver_warm_s = pver_warm %>% group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                            levels = c(Taxa[1:6], "Other"))) %>% tally(n) 

pver_warm_s$fraction = pver_warm_s$n / sum(pver_warm_s$n) # Compute percentages
pver_warm_s$ymax = cumsum(pver_warm_s$fraction) # Compute the cumulative percentages (top of each rectangle)
pver_warm_s$ymin = c(0, head(pver_warm_s$ymax, n=-1)) # Compute the bottom of each rectangle
pver_warm_s$labelPosition <- (pver_warm_s$ymax + pver_warm_s$ymin) / 2
pver_warm_s$label = pver_warm_s$n
pver_warm_don1=ggplot(pver_warm_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

#pver med_warm
pver_med_warm=subset(poci_indi_temp_df_sig,index == 3)
pver_med_warm$Taxa=tax$Family[match(rownames(pver_med_warm), rownames(tax))]
pver_med_warm[is.na(pver_med_warm)] <- "Unclassified"
pver_med_warm_s = pver_med_warm %>%  group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                                      levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
pver_med_warm_s$fraction = pver_med_warm_s$n / sum(pver_med_warm_s$n) # Compute percentages
pver_med_warm_s$ymax = cumsum(pver_med_warm_s$fraction) # Compute the cumulative percentages (top of each rectangle)
pver_med_warm_s$ymin = c(0, head(pver_med_warm_s$ymax, n=-1)) # Compute the bottom of each rectangle
pver_med_warm_s$labelPosition <- (pver_med_warm_s$ymax + pver_med_warm_s$ymin) / 2
pver_med_warm_s$label = pver_med_warm_s$n
pver_med_warm_don2=ggplot(pver_med_warm_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

#pver med_cool
pver_med_cool=subset(poci_indi_temp_df_sig,index == 2)
pver_med_cool$Taxa=tax$Family[match(rownames(pver_med_cool), rownames(tax))]
pver_med_cool[is.na(pver_med_cool)] <- "Unclassified"
pver_med_cool_s = pver_med_cool %>%  group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                                      levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
pver_med_cool_s$fraction = pver_med_cool_s$n / sum(pver_med_cool_s$n) # Compute percentages
pver_med_cool_s$ymax = cumsum(pver_med_cool_s$fraction) # Compute the cumulative percentages (top of each rectangle)
pver_med_cool_s$ymin = c(0, head(pver_med_cool_s$ymax, n=-1)) # Compute the bottom of each rectangle
pver_med_cool_s$labelPosition <- (pver_med_cool_s$ymax + pver_med_cool_s$ymin) / 2
pver_med_cool_s$label = pver_med_cool_s$n
pver_med_cool_don3=ggplot(pver_med_cool_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

#pver cool
pver_cool=subset(poci_indi_temp_df_sig,index == 1)
pver_cool$Taxa=tax$Family[match(rownames(pver_cool), rownames(tax))]
pver_cool[is.na(pver_cool)] <- "Unclassified"
pver_cool_s = pver_cool %>%  group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                              levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
pver_cool_s$fraction = pver_cool_s$n / sum(pver_cool_s$n) # Compute percentages
pver_cool_s$ymax = cumsum(pver_cool_s$fraction) # Compute the cumulative percentages (top of each rectangle)
pver_cool_s$ymin = c(0, head(pver_cool_s$ymax, n=-1)) # Compute the bottom of each rectangle
pver_cool_s$labelPosition <- (pver_cool_s$ymax + pver_cool_s$ymin) / 2
pver_cool_s$label = pver_cool_s$n
pver_cool_don4=ggplot(pver_cool_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

pdf("outputs/pver_donuts_temp.pdf",  height = 1, width = 1, pointsize = 0.1)
pver_warm_don1+pver_med_warm_don2+pver_med_cool_don3+pver_cool_don4
dev.off()


#spis warm
spis_warm=subset(styl_indi_temp_df_sig,index == 4)
spis_warm$Taxa=tax$Family[match(rownames(spis_warm), rownames(tax))]
spis_warm[is.na(spis_warm)] <- "Unclassified"
spis_big_s = spis_warm %>% group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                             levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
spis_big_s$fraction = spis_big_s$n / sum(spis_big_s$n) # Compute percentages
spis_big_s$ymax = cumsum(spis_big_s$fraction) # Compute the cumulative percentages (top of each rectangle)
spis_big_s$ymin = c(0, head(spis_big_s$ymax, n=-1)) # Compute the bottom of each rectangle
spis_big_s$labelPosition <- (spis_big_s$ymax + spis_big_s$ymin) / 2
spis_big_s$label = spis_big_s$n
spis_warm_don1=ggplot(spis_big_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

#spis med_warm
spis_med_warm=subset(styl_indi_temp_df_sig,index == 3)
spis_med_warm$Taxa=tax$Family[match(rownames(spis_med_warm), rownames(tax))]
spis_med_warm[is.na(spis_med_warm)] <- "Unclassified"
spis_med_warm_s = spis_med_warm %>%  group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                                      levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
spis_med_warm_s$fraction = spis_med_warm_s$n / sum(spis_med_warm_s$n) # Compute percentages
spis_med_warm_s$ymax = cumsum(spis_med_warm_s$fraction) # Compute the cumulative percentages (top of each rectangle)
spis_med_warm_s$ymin = c(0, head(spis_med_warm_s$ymax, n=-1)) # Compute the bottom of each rectangle
spis_med_warm_s$labelPosition <- (spis_med_warm_s$ymax + spis_med_warm_s$ymin) / 2
spis_med_warm_s$label = spis_med_warm_s$n
spis_med_warm_don2=ggplot(spis_med_warm_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

#spis med_cool
spis_med_cool=subset(styl_indi_temp_df_sig,index == 2)
spis_med_cool$Taxa=tax$Family[match(rownames(spis_med_cool), rownames(tax))]
spis_med_cool[is.na(spis_med_cool)] <- "Unclassified"
spis_med_cool_s = spis_med_cool %>%  group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                                      levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
spis_med_cool_s$fraction = spis_med_cool_s$n / sum(spis_med_cool_s$n) # Compute percentages
spis_med_cool_s$ymax = cumsum(spis_med_cool_s$fraction) # Compute the cumulative percentages (top of each rectangle)
spis_med_cool_s$ymin = c(0, head(spis_med_cool_s$ymax, n=-1)) # Compute the bottom of each rectangle
spis_med_cool_s$labelPosition <- (spis_med_cool_s$ymax + spis_med_cool_s$ymin) / 2
spis_med_cool_s$label = spis_med_cool_s$n
spis_med_cool_don3=ggplot(spis_med_cool_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

#spis cool
spis_cool=subset(styl_indi_temp_df_sig,index == 1)
spis_cool$Taxa=tax$Family[match(rownames(spis_cool), rownames(tax))]
spis_cool[is.na(spis_cool)] <- "Unclassified"
spis_cool_s = spis_cool %>%  group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                              levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
spis_cool_s$fraction = spis_cool_s$n / sum(spis_cool_s$n) # Compute percentages
spis_cool_s$ymax = cumsum(spis_cool_s$fraction) # Compute the cumulative percentages (top of each rectangle)
spis_cool_s$ymin = c(0, head(spis_cool_s$ymax, n=-1)) # Compute the bottom of each rectangle
spis_cool_s$labelPosition <- (spis_cool_s$ymax + spis_cool_s$ymin) / 2
spis_cool_s$label = spis_cool_s$n
spis_cool_don4=ggplot(spis_cool_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

pdf("outputs/spis_donuts_temp.pdf",  height = 1, width = 1, pointsize = 0.1)
spis_warm_don1+spis_med_warm_don2+spis_med_cool_don3+spis_cool_don4
dev.off()



# white board
lst=list(pver_warm_s$Taxa,pver_med_warm_s$Taxa,pver_med_cool_s$Taxa,pver_cool_s$Taxa,spis_big_s$Taxa,spis_med_warm_s$Taxa,spis_med_cool_s$Taxa,spis_cool_s$Taxa)
unlist(lst)
unique(unlist(lst))


lst=unique(unlist(list(pver_warm_s$Taxa,pver_cool_s$Taxa,spis_big_s$Taxa,spis_cool_s$Taxa)))
sort(lst)




## donuts per cluster

#pver c2
pver_c2=subset(poci_indi_pop_df_sig,index == 1) # c2 because the clusters names are swapped! 
pver_c2$Taxa=tax$Family[match(rownames(pver_c2), rownames(tax))]
pver_c2[is.na(pver_c2)] <- "Unclassified"
pver_c2_s = pver_c2 %>% group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                             levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
pver_c2_s$fraction = pver_c2_s$n / sum(pver_c2_s$n) # Compute percentages
pver_c2_s$ymax = cumsum(pver_c2_s$fraction) # Compute the cumulative percentages (top of each rectangle)
pver_c2_s$ymin = c(0, head(pver_c2_s$ymax, n=-1)) # Compute the bottom of each rectangle
pver_c2_s$labelPosition <- (pver_c2_s$ymax + pver_c2_s$ymin) / 2
pver_c2_s$label = pver_c2_s$n
pver_c2_don=ggplot(pver_c2_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")


#spis big cluster
spis_big=subset(styl_indi_pop_df_sig,index == 5 | index == 48) # c2 because the clusters names are swapped! 
spis_big$Taxa=tax$Family[match(rownames(spis_big), rownames(tax))]
spis_big[is.na(spis_big)] <- "Unclassified"
spis_big_s = spis_big %>% group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                         levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
spis_big_s$fraction = spis_big_s$n / sum(spis_big_s$n) # Compute percentages
spis_big_s$ymax = cumsum(spis_big_s$fraction) # Compute the cumulative percentages (top of each rectangle)
spis_big_s$ymin = c(0, head(spis_big_s$ymax, n=-1)) # Compute the bottom of each rectangle
spis_big_s$labelPosition <- (spis_big_s$ymax + spis_big_s$ymin) / 2
spis_big_s$label = spis_big_s$n
spis_big_don=ggplot(spis_big_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")


#spis c2
spis_c2=subset(styl_indi_pop_df_sig,index == 2 ) # c2 because the clusters names are swapped! 
spis_c2$Taxa=tax$Family[match(rownames(spis_c2), rownames(tax))]
spis_c2[is.na(spis_c2)] <- "Unclassified"
spis_c2_s = spis_c2 %>% group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                                 levels = c(Taxa[1:6], "Other"))) %>% tally(n)
spis_c2_s$fraction = spis_c2_s$n / sum(spis_c2_s$n) # Compute percentages
spis_c2_s$ymax = cumsum(spis_c2_s$fraction) # Compute the cumulative percentages (top of each rectangle)
spis_c2_s$ymin = c(0, head(spis_c2_s$ymax, n=-1)) # Compute the bottom of each rectangle
spis_c2_s$labelPosition <- (spis_c2_s$ymax + spis_c2_s$ymin) / 2
spis_c2_s$label = spis_c2_s$n
spis_c2_don=ggplot(spis_c2_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

#spis c6
spis_c6=subset(styl_indi_pop_df_sig,index == 6 ) # c6 because the clusters names are swapped! 
spis_c6$Taxa=tax$Family[match(rownames(spis_c6), rownames(tax))]
spis_c6[is.na(spis_c6)] <- "Unclassified"
spis_c6_s = spis_c6 %>% group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                         levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
spis_c6_s$fraction = spis_c6_s$n / sum(spis_c6_s$n) # Compute percentages
spis_c6_s$ymax = cumsum(spis_c6_s$fraction) # Compute the cumulative percentages (top of each rectangle)
spis_c6_s$ymin = c(0, head(spis_c6_s$ymax, n=-1)) # Compute the bottom of each rectangle
spis_c6_s$labelPosition <- (spis_c6_s$ymax + spis_c6_s$ymin) / 2
spis_c6_s$label = spis_c6_s$n
spis_c6_don=ggplot(spis_c6_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")


pdf("outputs/spis_donuts_clusters.pdf",  height = 1, width = 1, pointsize = 0.1)
pver_c2_don+spis_big_don+spis_c2_don+spis_c6_don
dev.off()


lst=unique(unlist(list(pver_c2_s$Taxa,spis_big_s$Taxa,spis_c2_s$Taxa,spis_c6_s$Taxa)))
sort(lst)


#spis big cluster
spis_big=subset(styl_indi_pop_df_sig,index == 5 ) # c2 because the clusters names are swapped! 
spis_big$Taxa=tax$Family[match(rownames(spis_big), rownames(tax))]
spis_big[is.na(spis_big)] <- "Unclassified"
spis_big_s = spis_big %>% group_by(Taxa) %>% tally(sort = TRUE) %>% group_by(Taxa = factor(c(Taxa[1:6], rep("Other", n() - 6)),
                                                                                           levels = c(Taxa[1:6], "Other"))) %>% tally(n) 
spis_big_s$fraction = spis_big_s$n / sum(spis_big_s$n) # Compute percentages
spis_big_s$ymax = cumsum(spis_big_s$fraction) # Compute the cumulative percentages (top of each rectangle)
spis_big_s$ymin = c(0, head(spis_big_s$ymax, n=-1)) # Compute the bottom of each rectangle
spis_big_s$labelPosition <- (spis_big_s$ymax + spis_big_s$ymin) / 2
spis_big_s$label = spis_big_s$n
spis_big_don=ggplot(spis_big_s, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Taxa)) +  geom_rect() + 
  geom_text( x=2, aes(y=labelPosition, label=label), size=6) + scale_fill_manual(values = palA)  + 
  coord_polar(theta="y") +  xlim(c(2, 4)) + theme_void() + theme(legend.position = "none")

### text

length(intersect(rownames(pver_c2), rownames(pver_warm)))
