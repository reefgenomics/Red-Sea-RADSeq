###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 07b ADMIXTURE population structure - Pocillopora verrucosa #####

####  07b.01  Visualization of the cross validation error beteween ADMIXTURE run 
# libraries
library(tidyr)
library(ggplot2)

# set working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_ADMIXTURE_20210223")

# visualization CVE (10 fold) - 10 runs with random seed per each K tested
# for i in 1 2 3 4 5 6 7 8 9 10; do cd run_cve_10fold_"$i" && grep -h CV log*.out && cd .. ; done| awk '{print $4}' | paste -d, - - - - - - - - - - - - -  > pver.crossvalidationerror.runs.txt
k.clust <- c(10, 11, 12, 13, 1, 2, 3, 4, 5, 6, 7, 8, 9)

cv.10.run1 <- c(0.37212,0.38093,0.39032,0.40074,0.29852,0.30518,0.31257,0.31959,0.32685,0.33493,0.34350,0.35256,0.36166)
cv.10.run2 <- c(0.37168,0.38144,0.39009,0.40071,0.29856,0.30530,0.31283,0.31978,0.32724,0.33517,0.34315,0.35190,0.36126)
cv.10.run3 <- c(0.37134,0.38120,0.39063,0.40079,0.29856,0.30532,0.31256,0.31919,0.32639,0.33508,0.34383,0.35238,0.36153)
cv.10.run4 <- c(0.37077,0.38109,0.39039,0.40022,0.29854,0.30536,0.31262,0.31970,0.32713,0.33479,0.34314,0.35190,0.36149)
cv.10.run5 <- c(0.37105,0.38107,0.39000,0.40182,0.29857,0.30534,0.31268,0.31971,0.32679,0.33467,0.34326,0.35192,0.36158)
cv.10.run6 <- c(0.37114,0.38116,0.39156,0.40128,0.29853,0.30539,0.31276,0.31947,0.32697,0.33493,0.34370,0.35280,0.36182)
cv.10.run7 <- c(0.37128,0.38102,0.39177,0.39982,0.29853,0.30535,0.31262,0.31958,0.32728,0.33518,0.34298,0.35266,0.36175)
cv.10.run8 <- c(0.37139,0.38096,0.39072,0.39980,0.29854,0.30525,0.31281,0.31990,0.32664,0.33485,0.34307,0.35299,0.36200)
cv.10.run9 <- c(0.37141,0.38116,0.39124,0.40072,0.29852,0.30524,0.31270,0.31991,0.32634,0.33477,0.34340,0.35243,0.36176)
cv.10.run10 <- c(0.37025,0.38088,0.39123,0.40057,0.29851,0.30524,0.31284,0.31948,0.32695,0.33456,0.34309,0.35257,0.36138)

cv10.df <- data.frame(k.clust, cv.10.run1, cv.10.run2, cv.10.run3, cv.10.run4, cv.10.run5, cv.10.run6, cv.10.run7, cv.10.run8, cv.10.run9,cv.10.run10)
cv10.df$k.clust <- factor(cv10.df$k.clust)
cv10.df.long <- gather(cv10.df, cve.run, value, cv.10.run1:cv.10.run10, factor_key=TRUE)


pdf('pver_ADMIXTURE_bestK_CVerror-based_ALLmarkers_10runs_10fold.pdf', width = 6, height = 4)
ggplot(cv10.df.long, aes(x=k.clust, y=value)) + geom_point() +
  stat_summary(aes(y = value,group=1), fun.y=median, colour="blue", geom="line",group=1) +
  labs(y= "Cross-Validation Error (10 fold)", x="Number of clusters", title="P. verrucosa - ADMIXTURE All Markers \nCVE (10 fold) - Median CV") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 6))
dev.off()

####  07b.02  Visualization of the ADMIXTURE barplots

# After assessment of the CV distribution and value ranges we determine that --cv=10 lead to smaller cross validation values than --cv=5 (default)
# Nevertheless, the ADMIXTURE plots were inspected using PONG

# libraries
library(scales)

# set working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_ADMIXTURE_20210223/Pong_visualization")

# load the filtered strata
pver.14reef.strata <- read.delim("../../pver_vcftools_filtered_vcf/pver.ind.order.in.vcf.txt", header = F)
colnames(pver.14reef.strata) <- "INDIVIDUALS"
pver.14reef.strata$INDIVIDUALS <- gsub("_filtered-sorted", "", pver.14reef.strata$INDIVIDUALS)
pver.14reef.strata$STRATA <- gsub("P", "", gsub("(-[^-]+)-.*","\\1", pver.14reef.strata$INDIVIDUALS))
write.table(pver.14reef.strata$STRATA, "pver.ind2pop", col.names = F, quote = F, row.names = F)

# Get the individuals order based on the major cluster in each reef was attained based on run1 k2
pver.k2.majormode.run1.pong <- read.table("../run_cve_10fold_1/pver.296ind.LE.plink.2.Q")
colnames(pver.k2.majormode.run1.pong) <- c("PCL1", "PCL2")
#add columns of individuals and population strata (plink file and strata file have the individuals in the same order)
pver.k2.majormode.run1.pong$ind <- pver.14reef.strata$INDIVIDUALS
pver.k2.majormode.run1.pong$reef <- pver.14reef.strata$STRATA
pver.k2.majormode.run1.pong$reef <- factor(pver.k2.majormode.run1.pong$reef, levels = c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R2", "WAJ-R3", "YAN-R1", "YAN-R3", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R3"), ordered = T)

reefs <- split(pver.k2.majormode.run1.pong, pver.k2.majormode.run1.pong$reef) # split Q matrix based on reef levels
str(reefs)

# Generate the vector of individuals in the most appropriate order
pver.reefs.bin.list<-lapply(reefs,function(reef){
  ind <- reef[,3] # individuals ids to a vector
  reef.Qmatrix <- reef[,-c(3,4)] # remove ind and reef columns from the Qmatrix
  reef.majormode <- reef.Qmatrix[,names(sort(colSums(reef.Qmatrix), decreasing = TRUE))] # sort the columns in decreasing order of membership in the entire reef
  reef.inds <- cbind(ind, reef.majormode) # add inviduals back to the columns sorted Q matrix
  reef.ind.ord <- reef.inds %>% arrange(across(starts_with("CL"), desc)) # ordered in descending order the individuals based on their membership in each of the columns (clusters)
  return(as.character(reef.ind.ord$ind)) # returns the list of characters of the ordered individuals
})
pver.id.ordered <- Reduce(c,pver.reefs.bin.list) # create a vector by joining the lists of characters produced above

write.table(pver.id.ordered , "pver.ind.ordered.byclusters.txt", quote = F, row.names = F, col.names = F)


# Get the ancestry proportion matrices in the right order (runs CV10)
runs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
k.adm <- c(2,3,4,5,6)
for(run in runs){
  for(k in k.adm){
    s.k.cv10.r<- read.table(paste("../run_cve_10fold_",run,"/pver.296ind.LE.plink.",k,".Q", sep = ""), h=F)
    s.k.cv10.r$v <- factor(pver.14reef.strata$INDIVIDUALS, levels=pver.id.ordered, ordered=T)
    ordered.s.k.cv10.r <- s.k.cv10.r[order(s.k.cv10.r$v),]
    ordered.s.k.cv10.r <- ordered.s.k.cv10.r[,-(k+1)]
    write.table(format(ordered.s.k.cv10.r, scientific=F), file = paste("pver.allmarkers.cv10.r",run,".indordered.",k,".Q", sep = ""), quote = F, col.names = F, row.names = F, sep = " ")
  }
}

# run the pong command top compare plots
# to activate the environment were pong dependencies were installed use the following command
# source /home/buitracn/pong/bin/activate
# to close that environment use `deactivate`

system("python run_pong.py -m pver.filemap.cv10 -i pver.ind2pop -n pver.poporder -l test.col --dist_metric jaccard -v -s 0.99")
system("python run_pong.py -m pver.filemap.run10.cv10 -i pver.ind2pop -n pver.poporder -l test.col --dist_metric jaccard -v -s 0.99") # plotting only run10


######
# Generate the plot outside of PONG
colnames(pver.k2.majormode.run1.pong) <- c("PCL1", "PCL2", "ind", "reef")

# I will use the same color pallete used in pong
pver.clust.colors <- c("#9370DB", "#00CED1")  #"#E1C340"
show_col(pver.clust.colors)

#Exchange the table format from wide to long
library(reshape2)
pver.k2 <- melt(pver.k2.majormode.run1.pong, id.vars=c("ind", "reef"))
pver.k2$reef <- factor(pver.k2$reef, levels=c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R2", "WAJ-R3", "YAN-R1", "YAN-R3", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R3"), ordered = T)
pver.k2$ind <- factor(pver.k2$ind, levels = pver.id.ordered, ordered= T)

k2 <- ggplot(pver.k2, aes(fill=variable, y=value,x=ind))+
  #geom_bar(stat="identity", position="fill", colour="grey33", size=0.0001, width = 1) + # size controls the thicknes of the outline of each bar #geom_col can also be used as a command to plot a stacked bar chart
  geom_bar(stat="identity", position="fill", size=0, width = 1) + # size controls the thicknes of the outline of each bar #geom_col can also be used as a command to plot a stacked bar chart
  facet_grid (~reef, scales = "free", space = "free_x") + #codigo para introducir los grupos. dependiendo de si .~ esta antes o despues del factor la figura es vertical u horizontal
  scale_y_continuous(expand = c(0,0))+ #this command helps to remove gray space below and above plot
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        axis.title.x= element_text(size=8),
        axis.text.x= element_text(size=6, angle = 90),
        axis.ticks.x=element_blank(),
        axis.ticks.x.top=element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size= 8),
        legend.text = element_text(size = 7),
        #legend.position = "none",
        panel.spacing.x= unit(0.3, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6))+
  theme(panel.spacing.x= unit(0.2, "lines"))+
  labs( y= "Membership probability", x="Individuals")+
  scale_fill_manual(values=pver.clust.colors,
                    name= "Clusters\nADMIXTURE")+ #codigo para utilizar una escala de colores personalizada
  guides(fill=guide_legend(ncol=1, keywidth = 0.5, keyheight = 0.5)) # to reduce the size of the legen key
k2
ggsave(filename = "TEST_k2.majormoderun1pong.pdf", plot = k2, width = 70, height = 7,  units = "cm")


####  07b.03  Visualization of the genetic divergence between ancestral populations 
# Libraries
library(ape)

# Working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_ADMIXTURE_20210223/ADMIXTURE_Clusters_divergence_K2_run1")

# read ancestral population divergences (FST) calculated in ADMIXTURE
pver.ancestry.k2 <- read.delim("ADMIXTURE_geneticclusterFST_K2_r1.txt", sep = "\t", header = F)
pver.ancestry.k2 <- as.matrix(pver.ancestry.k2)
rownames(pver.ancestry.k2) <- pver.ancestry.k2[,1]
pver.ancestry.k2 <- pver.ancestry.k2[,-1]
colnames(pver.ancestry.k2) <- rownames(pver.ancestry.k2) 

# Phylogenetic tree
tree.k2 <- nj(pver.ancestry.k2)
class(tree.k2)
## [1] "phylo"
tree.k2 <- ladderize(tree.k2)
tree.k2

write.tree(tree.k2, file = "tree.k2.newick", append = FALSE,
           digits = 10, tree.names = FALSE)

myPal.k2 <- c("#00CED1", "#9370DB")
show_col(myPal.k2)

pdf("pver_ancestral_populations_unrootedtree_k2_FSTdivergence_based.pdf", width = 5, height = 4)
plot(tree.k2, type="unrooted", tip.color = myPal.k2, rotate.tree = 90)
title("Pver Ancestral populations - Unrooted NJ tree (k2)")
dev.off()


plot(tree.k2, type="cladogram", tip.color = myPal.k2)

pdf("Pver_ancestral_populations_unrootedtree_k2_FSTdivergence_based_phylogram.pdf", width = 5, height = 4)
plot(tree.k2, type="phylogram", tip.color = myPal.k2)
title("Pver Ancestral populations - Unrooted NJ tree (k2)")
dev.off()

####  07b.04  Identify individuals assigned to a unique genetic cluster (≥90 asignment probability to a unique genetic cluster) 

# Libraries
library(dplyr)
library(scales)

# Working directory K2
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_ADMIXTURE_20210223/pver_Individuals_split_by_cluster_20210302/K2")

 K=2 => pong representative run run7_allmarkers_cv10_K2 (this result is supported by 10 runs that converge to the same answer with an average similarity of 97.9%)
# The disrupt perm file shows the following order for the .Q column, I will asign the cluster ID based on the order of the columns in the file
# 2 #909595
# 1 #4FDEDD

pver.clust.colors <- c("#00CED1","#9370DB")
show_col(pver.clust.colors)

# Read the ancestry memebership matrix
#ln -s /home/buitracn/RADseq-Big-project/pver/p1.mac4.r0.8.316samples.Pver_20200218/04_ADMIXTURE_20210223/Pong_visualization/pver.allmarkers.cv10.r10.indordered.2.Q .
pver.admix.k2.r10.rep <- read.delim("pver.allmarkers.cv10.r10.indordered.2.Q", sep = "", header = F)

colnames(pver.admix.k2.r10.rep) <- paste("PCL", seq(1:2), sep = "") # bear in mind the match between color and column number

# ln -s /home/buitracn/RADseq-Big-project/pver/p1.mac4.r0.8.316samples.Pver_20200218/04_ADMIXTURE_20210223/Pong_visualization/pver.ind.ordered.byclusters.txt .
pver.id.ordered <- read.delim("pver.ind.ordered.byclusters.txt", header = F)
rownames(pver.admix.k2.r10.rep) <- pver.id.ordered$V1

pver.genclust <- mutate(pver.admix.k2.r10.rep, major.clust = case_when(PCL1 >= 0.9 ~ "PCL1",
                                                                   PCL2 >= 0.9 ~ "PCL2",
                                                                   TRUE ~ "Admix"))

pver.clustersID <- paste("PCL", seq(1:2), sep = "" )
for (i in pver.clustersID){
  cluster <- rownames(subset(pver.genclust, major.clust == i))
  write.table(cluster, paste(i,".txt", sep = ""), quote = F, col.names = F, row.names = F)
}

table(pver.genclust$major.clust)
# Admix   PCL1   PCL2 
# 93    24   179 

# Individuals classified as admixed
# pver.genclust[grep("Admix", pver.genclust$major.clust), ]

# create strata with the new clusters
pver.genclust.strata <- data.frame(INDIVIDUALS = rownames(pver.genclust),
                                   STRATA = pver.genclust$major.clust) 
pver.genclust.strata <- pver.genclust.strata[ grep("Admix", pver.genclust.strata$STRATA, invert = TRUE) , ] #203 individuals

write.table(pver.genclust.strata, file = "pver.genclust.strata.K2.tsv", sep = "\t", quote = F, row.names = F)

