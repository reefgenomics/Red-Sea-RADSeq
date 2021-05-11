setwd("~/Documents/Bioinformatics_scripts/R_scripts/Spis_Pver/")

library(data.table)
library(tidyr)
library(vegan)
library(pairwiseAdonis)
source("~/Documents/Bioinformatics_scripts/R_scripts/orphan/remove_rare.R")

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
#asv_filt2=remove_rare(asv_filt, 0.025) #30 samples (0.05) -> 679 ASVs, # 15 samples (0.025) -> 1317
#asv_filt2=asv_filt
asv_number=nrow(asv_filt)

asv.n=as.data.frame(t(sweep(asv_filt2,2,colSums(asv_filt2),`/`)))
#asv.n=as.data.frame(t(apply(asv_filt2,2, clr)))

asv.n$PopID=paste(host$STRATA)[match(rownames(asv.n), host$INDIVIDUALS)]
asv.n$Species=ifelse(rownames(asv.n) %like% "^S", "Stylophora", "Pocillopora")
asv.n$Reef=paste(map$site, map$Reef, sep = "_")[match(rownames(asv.n), rownames(map))] 
asv.n$Temperature=temp$temp_category2[match(asv.n$Reef, temp$Reef)]

poci=subset(asv.n, Species == "Pocillopora")
styl=subset(asv.n, Species == "Stylophora")

## permanovas overall 
#method_adonis="euclidean"
method_adonis="bray"

adonis(poci[,1:asv_number] ~ poci$Temperature * poci$PopID, method = method_adonis)
adonis(poci[,1:asv_number] ~ poci$PopID * poci$Temperature, method = method_adonis )
adonis(styl[,1:asv_number] ~ styl$Temperature * styl$PopID, method = method_adonis)
adonis(styl[,1:asv_number] ~ styl$PopID * styl$Temperature, method = method_adonis)


## permanovas overall + strata
adonis(poci[,1:asv_number] ~ poci$PopID, strata = poci$Temperature ,method = method_adonis)
adonis(styl[,1:asv_number] ~ styl$PopID, strata = styl$Temperature ,method = method_adonis)


