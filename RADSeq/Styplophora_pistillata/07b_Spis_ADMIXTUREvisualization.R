###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 07b ADMIXTURE population structure - Stylophora pistillata #####

####  07b.01  Visualization of the cross validation error beteween ADMIXTURE run 

# libraries
library(tidyr)
library(ggplot2)

# set working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_ADMIXTURE_20210222")

# visualization CVE (10 fold) - 10 runs with random seed per each K tested
# for i in 1 2 3 4 5 6 7 8 9 10; do cd run_cve_10fold_"$i" && grep -h CV log*.out && cd .. ; done| awk '{print $4}' | paste -d, - - - - - - - - - - - - - - - - - > spis.crossvalidationerror.runs.txt
k.clust <- c(10, 11, 12, 13, 14, 15, 16, 17, 1, 2, 3, 4, 5, 6, 7, 8, 9)
cv.10.run1 <- c(0.21361,0.22120,0.21713,0.21818,0.21943,0.21970,0.22234,0.22928,0.26978,0.24828,0.23601,0.23000,0.22629,0.22346,0.21751,0.21745,0.21692)
cv.10.run2 <- c(0.21699,0.21511,0.21643,0.21864,0.21951,0.22335,0.22416,0.22628,0.26980,0.24818,0.23613,0.23156,0.22628,0.22091,0.21936,0.21725,0.21341)
cv.10.run3 <- c(0.21496,0.21497,0.21725,0.21923,0.21933,0.22157,0.22577,0.22769,0.26981,0.24791,0.23607,0.23060,0.22891,0.22143,0.22036,0.21550,0.21355)
cv.10.run4 <- c(0.21493,0.21489,0.21655,0.21800,0.21892,0.22015,0.22580,0.22583,0.26980,0.24791,0.24073,0.23071,0.22569,0.22120,0.21718,0.21541,0.21365)
cv.10.run5 <- c(0.21490,0.21691,0.21729,0.21724,0.21948,0.22433,0.22395,0.22586,0.26981,0.24795,0.23601,0.23069,0.22694,0.22153,0.21934,0.21714,0.21358)
cv.10.run6 <- c(0.21705,0.21736,0.21739,0.21784,0.21886,0.22290,0.22555,0.22874,0.26982,0.24793,0.23609,0.23067,0.22669,0.22158,0.22253,0.21568,0.21701)
cv.10.run7 <- c(0.21624,0.21490,0.21835,0.21904,0.21888,0.22049,0.22512,0.22705,0.26977,0.24788,0.24001,0.23001,0.22781,0.22313,0.21941,0.21545,0.21344)
cv.10.run8 <- c(0.21365,0.21542,0.21744,0.21769,0.22102,0.22101,0.22724,0.22761,0.26983,0.24799,0.24036,0.23152,0.22537,0.22199,0.21714,0.21537,0.21557)
cv.10.run9 <- c(0.21375,0.21554,0.21690,0.21810,0.21972,0.22084,0.22406,0.22678,0.26981,0.24789,0.23602,0.23005,0.22655,0.22421,0.21994,0.21821,0.21689)
cv.10.run10 <- c(0.21470,0.21718,0.21554,0.21923,0.21887,0.22329,0.22467,0.22665,0.26978,0.24831,0.23610,0.23005,0.22919,0.22390,0.21965,0.21561,0.21700)

cv10.df <- data.frame(k.clust, cv.10.run1, cv.10.run2, cv.10.run3, cv.10.run4, cv.10.run5, cv.10.run6, cv.10.run7, cv.10.run8, cv.10.run9,cv.10.run10)
cv10.df$k.clust <- factor(cv10.df$k.clust)
cv10.df.long <- gather(cv10.df, cve.run, value, cv.10.run1:cv.10.run10, factor_key=TRUE)


pdf('spis_ADMIXTURE_bestK_CVerror-based_ALLmarkers_10runs_10fold.pdf', width = 6, height = 4)
ggplot(cv10.df.long, aes(x=k.clust, y=value)) + geom_point() +
  stat_summary(aes(y = value,group=1), fun.y=median, colour="blue", geom="line",group=1) +
  labs(y= "Cross-Validation Error (10 fold)", x="Number of clusters", title="S. pistillata - ADMIXTURE All Markers \nCVE (10 fold) - Median CV") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = 6))
dev.off()

####  07b.02  Visualization of the ADMIXTURE barplots

# After assement of the CV distribution and value ranges we determine that --cv=10 lead to smaller cross validation values than --cv=5 (default)
# Nevetheless the ADMIXTURE plots were inspected using PONG

# libraries
library(scales)

# set working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_ADMIXTURE_20210222/Pong_visualization")

# load the filtered strata
spis.strata <- read.delim("../../spis_vcftools_filtered_vcf/spis.ind.order.in.vcf.txt", header = F)
colnames(spis.strata) <- "INDIVIDUALS"
spis.strata$REEF <- gsub("S","", sub('^([^-]+-[^-]+).*', '\\1', spis.strata$INDIVIDUALS))
spis.strata$REEF <- factor(spis.strata$REEF, levels = c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R3", "WAJ-R4", "YAN-R1", "YAN-R3", "YAN-R4", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R1", "FAR-R2", "FAR-R3", "FAR-R4"), ordered = T)
write.table(spis.strata$REEF, "spis.ind2pop", col.names = F, quote = F, row.names = F)

# id order genetic assignment (sample SKAU-R2-13 was discarded in comparison with previous run)
# Get the individuals order based on the major cluster in each reef was attained based on run5 K7
spis.k7.majormode.run5.pong <- read.table("../run_cve_10fold_5/spis.367ind.LE.plink.7.Q")
colnames(spis.k7.majormode.run5.pong) <- c("SCL1", "SCL2", "SCL3", "SCL4", "SCL5", "SCL6", "SCL7")
spis.k7.majormode.run5.pong$ind <- spis.strata$INDIVIDUALS
spis.k7.majormode.run5.pong$reef <- spis.strata$REEF
spis.k7.majormode.run5.pong$reef <- factor(spis.k7.majormode.run5.pong$reef, levels = c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R3", "WAJ-R4", "YAN-R1", "YAN-R3", "YAN-R4", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R1", "FAR-R2", "FAR-R3", "FAR-R4"), ordered = T)

reefs <- split(spis.k7.majormode.run5.pong, spis.k7.majormode.run5.pong$reef) # split Q matrix based on reef levels
str(reefs)

# Generate the vector of individuals in the most appropriate order
spis.reefs.bin.list<-lapply(reefs,function(reef){
  ind <- reef[,8] # individuals ids to a vector
  reef.Qmatrix <- reef[,-c(8,9)] # remove ind and reef columns from the Qmatrix
  reef.majormode <- reef.Qmatrix[,names(sort(colSums(reef.Qmatrix), decreasing = TRUE))] # sort the columns in decreasing order of membership in the entire reef
  reef.inds <- cbind(ind, reef.majormode) # add individuals back to the columns sorted Q matrix
  reef.ind.ord <- reef.inds %>% arrange(across(starts_with("CL"), desc)) # ordered in descending order the individuals based on their membership in each of the columns (clusters)
  return(as.character(reef.ind.ord$ind)) # returns the list of characters of the ordered individuals
})
spis.id.ordered <- Reduce(c,spis.reefs.bin.list) # create a vector by joining the lists of characters produced above

# some individuals are still in and undesired order. Use the following command to print the vector itself ==> dput(spis.id.ordered)
# I will manually re arrange the 1 individuals that seem odd (SYAN-R1-16)
spis.id.ordered <- c("SMAQ-R1-30","SMAQ-R1-1", "SMAQ-R1-2", "SMAQ-R1-3", "SMAQ-R1-5", "SMAQ-R1-7", "SMAQ-R1-8", "SMAQ-R1-10","SMAQ-R1-12", "SMAQ-R1-14","SMAQ-R1-15","SMAQ-R1-18","SMAQ-R1-23","SMAQ-R1-29","SMAQ-R1-11","SMAQ-R1-19","SMAQ-R1-25","SMAQ-R1-6","SMAQ-R1-17","SMAQ-R1-4", "SMAQ-R1-24","SMAQ-R1-28",
                     "SMAQ-R2-33","SMAQ-R2-7", "SMAQ-R2-9", "SMAQ-R2-13","SMAQ-R2-14", "SMAQ-R2-17","SMAQ-R2-18","SMAQ-R2-37","SMAQ-R2-28","SMAQ-R2-22","SMAQ-R2-16","SMAQ-R2-6", "SMAQ-R2-3", "SMAQ-R2-40", "SMAQ-R2-8", "SMAQ-R2-20","SMAQ-R2-29","SMAQ-R2-30","SMAQ-R2-31","SMAQ-R2-10","SMAQ-R2-11","SMAQ-R2-27",
                     "SWAJ-R1-23", "SWAJ-R1-25","SWAJ-R1-26","SWAJ-R1-27","SWAJ-R1-30","SWAJ-R1-33","SWAJ-R1-34","SWAJ-R1-38","SWAJ-R1-42","SWAJ-R1-43", "SWAJ-R1-45","SWAJ-R1-29",
                     "SWAJ-R3-1", "SWAJ-R3-2", "SWAJ-R3-3", "SWAJ-R3-4", "SWAJ-R3-5", "SWAJ-R3-6", "SWAJ-R3-7","SWAJ-R3-8", "SWAJ-R3-9", "SWAJ-R3-10","SWAJ-R3-11","SWAJ-R3-12","SWAJ-R3-13","SWAJ-R3-14","SWAJ-R3-15","SWAJ-R3-16", "SWAJ-R3-17","SWAJ-R3-21","SWAJ-R3-23","SWAJ-R3-24","SWAJ-R3-18","SWAJ-R3-19","SWAJ-R3-22",
                     "SWAJ-R4-4", "SWAJ-R4-6", "SWAJ-R4-7", "SWAJ-R4-15","SWAJ-R4-20","SWAJ-R4-21","SWAJ-R4-26","SWAJ-R4-19","SWAJ-R4-24","SWAJ-R4-23","SWAJ-R4-18", "SWAJ-R4-25","SWAJ-R4-14","SWAJ-R4-17","SWAJ-R4-1", "SWAJ-R4-3", "SWAJ-R4-5", "SWAJ-R4-8", "SWAJ-R4-9", "SWAJ-R4-10", "SWAJ-R4-11","SWAJ-R4-13","SWAJ-R4-27","SWAJ-R4-30","SWAJ-R4-2", 
                     "SYAN-R1-1", "SYAN-R1-9", "SYAN-R1-10","SYAN-R1-11", "SYAN-R1-12","SYAN-R1-14","SYAN-R1-19","SYAN-R1-20","SYAN-R1-13","SYAN-R1-16","SYAN-R1-15","SYAN-R1-6", "SYAN-R1-30","SYAN-R1-2", "SYAN-R1-7", "SYAN-R1-5", "SYAN-R1-3", "SYAN-R1-29","SYAN-R1-18",
                     "SYAN-R3-2", "SYAN-R3-3", "SYAN-R3-4", "SYAN-R3-7", "SYAN-R3-9", "SYAN-R3-11","SYAN-R3-13","SYAN-R3-16","SYAN-R3-18","SYAN-R3-24","SYAN-R3-6", "SYAN-R3-25", "SYAN-R3-1", "SYAN-R3-12","SYAN-R3-26","SYAN-R3-20","SYAN-R3-23","SYAN-R3-19","SYAN-R3-15","SYAN-R3-21","SYAN-R3-10", "SYAN-R3-14","SYAN-R3-8", "SYAN-R3-17","SYAN-R3-22",
                     "SYAN-R4-2", "SYAN-R4-3", "SYAN-R4-4", "SYAN-R4-5", "SYAN-R4-6", "SYAN-R4-10","SYAN-R4-12","SYAN-R4-18","SYAN-R4-21","SYAN-R4-23","SYAN-R4-25","SYAN-R4-13","SYAN-R4-19","SYAN-R4-15", "SYAN-R4-8", "SYAN-R4-17","SYAN-R4-24","SYAN-R4-20","SYAN-R4-14","SYAN-R4-9", "SYAN-R4-16","SYAN-R4-22","SYAN-R4-1", "SYAN-R4-7", "SYAN-R4-11",
                     "SKAU-R1-4", "SKAU-R1-6", "SKAU-R1-8", "SKAU-R1-21","SKAU-R1-27","SKAU-R1-28","SKAU-R1-23", "SKAU-R1-16",
                     "SKAU-R2-3", "SKAU-R2-4", "SKAU-R2-7", "SKAU-R2-33","SKAU-R2-28","SKAU-R2-32","SKAU-R2-26","SKAU-R2-23", "SKAU-R2-8", "SKAU-R2-10","SKAU-R2-25","SKAU-R2-2", "SKAU-R2-1", "SKAU-R2-5", "SKAU-R2-22","SKAU-R2-31","SKAU-R2-16", "SKAU-R2-19",
                     "SKAU-R3-33","SKAU-R3-37","SKAU-R3-39","SKAU-R3-40","SKAU-R3-41","SKAU-R3-43","SKAU-R3-44","SKAU-R3-47", "SKAU-R3-48","SKAU-R3-50","SKAU-R3-51","SKAU-R3-56","SKAU-R3-58","SKAU-R3-32","SKAU-R3-59","SKAU-R3-36","SKAU-R3-35", "SKAU-R3-52","SKAU-R3-54","SKAU-R3-55",
                     "SDOG-R1-2", "SDOG-R1-26","SDOG-R1-20","SDOG-R1-18","SDOG-R1-22","SDOG-R1-16", "SDOG-R1-21","SDOG-R1-7", "SDOG-R1-24","SDOG-R1-28","SDOG-R1-30","SDOG-R1-29","SDOG-R1-25","SDOG-R1-19",
                     "SDOG-R2-1", "SDOG-R2-6", "SDOG-R2-14","SDOG-R2-16","SDOG-R2-17","SDOG-R2-22","SDOG-R2-24","SDOG-R2-28","SDOG-R2-3", "SDOG-R2-5", "SDOG-R2-8", "SDOG-R2-18","SDOG-R2-26","SDOG-R2-11","SDOG-R2-13","SDOG-R2-4", "SDOG-R2-21","SDOG-R2-29","SDOG-R2-19", "SDOG-R2-30","SDOG-R2-15",
                     "SDOG-R3-9", "SDOG-R3-13","SDOG-R3-15","SDOG-R3-25","SDOG-R3-28","SDOG-R3-23","SDOG-R3-19", "SDOG-R3-21","SDOG-R3-2", "SDOG-R3-5", "SDOG-R3-11","SDOG-R3-12","SDOG-R3-4", "SDOG-R3-7", "SDOG-R3-26","SDOG-R3-17", "SDOG-R3-24","SDOG-R3-6", "SDOG-R3-10","SDOG-R3-18","SDOG-R3-3", "SDOG-R3-8", "SDOG-R3-14","SDOG-R3-16","SDOG-R3-20", 
                     "SFAR-R1-2", "SFAR-R1-4", "SFAR-R1-5", "SFAR-R1-9", "SFAR-R1-10","SFAR-R1-15","SFAR-R1-16","SFAR-R1-18","SFAR-R1-19", "SFAR-R1-20","SFAR-R1-24","SFAR-R1-25","SFAR-R1-17","SFAR-R1-3", "SFAR-R1-7", "SFAR-R1-14","SFAR-R1-29","SFAR-R1-22", "SFAR-R1-27","SFAR-R1-8", "SFAR-R1-28","SFAR-R1-23","SFAR-R1-1", "SFAR-R1-6", "SFAR-R1-11",
                     "SFAR-R2-1", "SFAR-R2-5", "SFAR-R2-10","SFAR-R2-12","SFAR-R2-16","SFAR-R2-23","SFAR-R2-26","SFAR-R2-2", "SFAR-R2-27","SFAR-R2-7", "SFAR-R2-14", "SFAR-R2-15","SFAR-R2-20","SFAR-R2-25","SFAR-R2-13","SFAR-R2-17","SFAR-R2-18","SFAR-R2-21","SFAR-R2-3", 
                     "SFAR-R3-1", "SFAR-R3-5", "SFAR-R3-10","SFAR-R3-15","SFAR-R3-17","SFAR-R3-2", "SFAR-R3-8", "SFAR-R3-3", "SFAR-R3-6", "SFAR-R3-11", "SFAR-R3-13","SFAR-R3-14","SFAR-R3-18","SFAR-R3-9", "SFAR-R3-16","SFAR-R3-20","SFAR-R3-21","SFAR-R3-12",
                     "SFAR-R4-1", "SFAR-R4-5", "SFAR-R4-10","SFAR-R4-11","SFAR-R4-28","SFAR-R4-26","SFAR-R4-14","SFAR-R4-15","SFAR-R4-16","SFAR-R4-21", "SFAR-R4-23","SFAR-R4-2", "SFAR-R4-19","SFAR-R4-3", "SFAR-R4-9", "SFAR-R4-20","SFAR-R4-22","SFAR-R4-24","SFAR-R4-30", "SFAR-R4-4", "SFAR-R4-6", "SFAR-R4-17","SFAR-R4-18","SFAR-R4-7", "SFAR-R4-8", "SFAR-R4-25")

write.table(spis.id.ordered , "spis.ind.ordered.byclusters.txt", quote = F, row.names = F, col.names = F)

# Get the ancestry proportion matrices in the right order (runs CV10)
runs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
k.adm <- c(2,3,4,5,6,7,8,9,10, 11, 12)
for(run in runs){
  for(k in k.adm){
    s.k.cv10.r<- read.table(paste("../run_cve_10fold_",run,"/spis.367ind.LE.plink.",k,".Q", sep = ""), h=F)
    s.k.cv10.r$v <- factor(spis.strata$INDIVIDUALS, levels=spis.id.ordered, ordered=T)
    ordered.s.k.cv10.r <- s.k.cv10.r[order(s.k.cv10.r$v),]
    ordered.s.k.cv10.r <- ordered.s.k.cv10.r[,-(k+1)]
    write.table(format(ordered.s.k.cv10.r, scientific=F), file = paste("spis.allmarkers.cv10.r",run,".indordered.",k,".Q", sep = ""), quote = F, col.names = F, row.names = F, sep = " ")
  }
}

# run the pong command top compare plots
# to activate the environment were pong dependencies were installed use the following command
# source /home/buitracn/pong/bin/activate
# to close that environment use `deactivate`

# I replaced grey for cyan
system("python run_pong.py -m spis.filemap.cv10 -i spis.ind2pop -n spis.poporder -l color.list.k2.k11.cv10 --dist_metric jaccard -v -s 0.95")
system("python run_pong.py -m spis.filemap.run9.cv10 -i spis.ind2pop -n spis.poporder -l color.list.k2.k6.cv10.run9 --dist_metric jaccard -v -s 0.95") # plotting only run1 (only K2 to K6)

######
# Figuring out the best order of the samples checking at run5 K7
spis.k7.majormode.run5.ordered.pong <- read.table("spis.allmarkers.cv10.r5.indordered.7.Q")
spis.k7.majormode.run5.ordered.pong$ind <- spis.id.ordered
spis.k7.majormode.run5.ordered.pong$reef <- gsub("S","", sub('^([^-]+-[^-]+).*', '\\1', spis.k7.majormode.run5.ordered.pong$ind))
spis.k7.majormode.run5.ordered.pong$reef <- factor(spis.k7.majormode.run5.ordered.pong$reef, levels = c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R3", "WAJ-R4", "YAN-R1", "YAN-R3", "YAN-R4", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R1", "FAR-R2", "FAR-R3", "FAR-R4"), ordered = T)

colnames(spis.k7.majormode.run5.ordered.pong) <- c("CL1", "CL2", "CL3", "CL4", "CL5", "CL6", "CL7", "ind", "reef")

# I will use the same color pallete used in pong
spis.clust.colors <- c("#FE92CD", "#C71585", "#32CD32", "#FFD700", "#1E90FF" , "#FF4500","#FD9202")
#check.colors.pong <- c("#FE92CD", "#C71585", "#32CD32", "#FFD700", "#1E90FF" , "#FF4500","#FD9202", "#11E8F2","#FFFAF0","#7B68EE","#000000", "#67BDF9")
show_col(spis.clust.colors)

# Change the table format from wide to long
library(reshape2)
spis.k7 <- melt(spis.k7.majormode.run5.ordered.pong, id.vars=c("ind", "reef"))
spis.k7$reef <- factor(spis.k7$reef, levels=c("MAQ-R1", "MAQ-R2", "WAJ-R1", "WAJ-R3", "WAJ-R4", "YAN-R1", "YAN-R3", "YAN-R4", "KAU-R1", "KAU-R2", "KAU-R3", "DOG-R1", "DOG-R2", "DOG-R3", "FAR-R1", "FAR-R2", "FAR-R3", "FAR-R4"), ordered = T)
spis.k7$ind <- factor(spis.k7$ind, levels = spis.id.ordered, ordered= T)

k7 <- ggplot(spis.k7, aes(fill=variable, y=value,x=ind))+
  #geom_bar(stat="identity", position="fill", colour="grey33", size=0.0001, width = 1) + # size controls the thickness of the outline of each bar #geom_col can also be used as a command to plot a stacked bar chart
  geom_bar(stat="identity", position="fill", size=0, width = 1) + # size controls the thickness of the outline of each bar #geom_col can also be used as a command to plot a stacked bar chart
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
  scale_fill_manual(values=spis.clust.colors,
                    name= "Clusters\nADMIXTURE")+ #codigo para utilizar una escala de colores personalizada
  guides(fill=guide_legend(ncol=1, keywidth = 0.5, keyheight = 0.5)) # to reduce the size of the legen key
k7
ggsave(filename = "TEST_k7.majormoderun5pong.pdf", plot = k7, width = 70, height = 7,  units = "cm")


####  07b.03  Visualization of the genetic divergence between ancestral populations 

# Libraries
library(ape)

# Working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_ADMIXTURE_20210222/ADMIXTURE_Clusters_divergence_K6_run9")

# read ancestral population divergences (FST) calculated in ADMIXTURE
spis.ancestry.k6 <- read.delim("ADMIXTURE_geneticclusterFST_K6_r9.txt", sep = "\t", header = F)
spis.ancestry.k6 <- as.matrix(spis.ancestry.k6)
rownames(spis.ancestry.k6) <- spis.ancestry.k6[,1]
spis.ancestry.k6 <- spis.ancestry.k6[,-1]
colnames(spis.ancestry.k6) <- rownames(spis.ancestry.k6) 

# Phylogenetic tree
tree.k6 <- nj(spis.ancestry.k6)
class(tree.k6)
## [1] "phylo"
tree.k6 <- ladderize(tree.k6)
tree.k6

write.tree(tree.k6, file = "tree.k6.newick", append = FALSE,
           digits = 10, tree.names = FALSE)

myPal.k6 <- c("#FE92CD", "#1E90FF","#FF4500", "#C71585", "#FFD700", "#32CD32")
show_col(myPal.k6)

pdf("Spis_ancestral_populations_unrootedtree_k6_FSTdivergence_based.pdf", width = 5, height = 4)
plot(tree.k6, type="unrooted", tip.color = myPal.k6, rotate.tree = 90)
title("Spis Ancestral populations - Unrooted NJ tree (k6)")
dev.off()


plot(tree.k6, type="cladogram", tip.color = myPal.k6)

pdf("Spis_ancestral_populations_unrootedtree_k6_FSTdivergence_based_phylogram.pdf", width = 5, height = 4)
plot(tree.k6, type="phylogram", tip.color = myPal.k6)
title("Spis Ancestral populations - Unrooted NJ tree (k6)")
dev.off()


####  07b.04  Identify individuals assigned to a unique genetic cluster (≥90 asignment probability to a unique genetic cluster) 

# Libraries
library(dplyr)
library(scales)

# Working directory K6
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_Individuals_split_by_cluster_20210302/K6")

# K=6 => pong representative run run9_allmarkers_cv10_K6 (this result is supported by 2 runs that converge to the same answer with an average similarity of 99.9%)
# The disrupt perm file shows the following order for the .Q column, I will asign the cluster ID based on the order of the columns in the file
#2 #1E90FF blue
#5 #FFD700 yellow
#4 #C71585 fuchsia
#1 #FE92CD pink
#3 #FF4500 red
#6 #32CD32 green

spis.clust.colors <- c("#FE92CD", "#1E90FF","#FF4500", "#C71585", "#FFD700", "#32CD32")
show_col(spis.clust.colors)

# color used in pong 
pong.col <- c("#1E90FF", "#32CD32", "#FFD700", "#FE92CD", "#C71585", "#FF4500", "#A5AEB4", "#FD9202", "#FFFAF0", "#7B68EE", "#000000") 


# Read the ancestry memebership matrix
system("ln -s /home/buitracn/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_ADMIXTURE_20210222/Pong_visualization/spis.allmarkers.cv10.r9.indordered.7.Q .")
spis.admix.K6.r9.rep <- read.delim("spis.allmarkers.cv10.r9.indordered.6.Q", sep = "", header = F)

colnames(spis.admix.K6.r9.rep) <- paste("SCL", seq(1:6), sep = "") # bear in mind the match between color and column number

# ln -s /home/buitracn/RADseq-Big-project/spis/p1_mac4_r80_448samples/04_ADMIXTURE_20210222/Pong_visualization/spis.ind.ordered.byclusters.txt .
spis.id.ordered <- read.delim("spis.ind.ordered.byclusters.txt", header = F)
rownames(spis.admix.K6.r9.rep) <- spis.id.ordered$V1

spis.genclust <- mutate(spis.admix.K6.r9.rep, major.clust = case_when(SCL1 >= 0.9 ~ "CL1",
                                                                      SCL2 >= 0.9 ~ "CL2",
                                                                      SCL3 >= 0.9 ~ "CL3",
                                                                      SCL4 >= 0.9 ~ "CL4",
                                                                      SCL5 >= 0.9 ~ "CL5",
                                                                      SCL6 >= 0.9 ~ "CL6",
                                                                      TRUE ~ "Admix"))

spis.clustersID <- paste("SCL", seq(1:6), sep = "" )
for (i in spis.clustersID){
  cluster <- rownames(subset(spis.genclust, major.clust == i))
  write.table(cluster, paste(i,".txt", sep = ""), quote = F, col.names = F, row.names = F)
}

table(spis.genclust$major.clust)
# Admix   SCL1   SCL2   SCL3   SCL4   SCL5   SCL6 
#  55     47      50     51     71     34      59 

# Within the group classified as admixed individuals (more than 10% of a secundary cluster) there are SKAU-R2=4, SDOG-R1=7, SFAR-R1=3, SFAR-R2=7, SFAR-R3=5, SFAR-R4=4
# spis.genclust[grep("Admix", spis.genclust$major.clust), ]

# create strata with the new clusters
spis.genclust.strata <- data.frame(INDIVIDUALS = rownames(spis.genclust),
                                   STRATA = spis.genclust$major.clust) 
spis.genclust.strata <- spis.genclust.strata[ grep("Admix", spis.genclust.strata$STRATA, invert = TRUE) , ] #312 individuals

write.table(spis.genclust.strata, file = "spis.genclust.strata.K6.tsv", sep = "\t", quote = F, row.names = F)

