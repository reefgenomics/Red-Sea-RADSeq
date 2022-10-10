###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 09b Candidate SNPs for Positive Selection with BAYPASS  - Stylophora pistillata #####

#### 09b.01 BAYPASS run

# BAYPASS - covariace dependent model to identify loci under selection
# a scaled covariance matrix is estimated to account for the historical relationship of the populations based on their allele frequencies.
# The variance–covariance matrix is used to control for evolutionary history in the calculation of XTX for each SNP

# Export baypass to the PATH
# export PATH="/home/buitracn/RADseq/Genotyping/tools/baypass_2.2/sources:$PATH"


# Better suited than covariance-free models for population with complex hierarchical structure

# working directory
setwd("/home/buitracn/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_positive_selection_analyses_20210308/BAYPASS")

# Run the 3 different MCMC chains (using the core model) and compare the covariance matrix "omega" between runs to assure convergence
system("g_baypass -gfile pver.popgenclust.baypass -outprefix Pver_1 -nthreads 10 -seed 658913 2>&1 | tee Pver_1_baypass_coremodel_defaultparm_log.txt")

system("g_baypass -gfile pver.popgenclust.baypass -outprefix Pver_2 -nthreads 10 -seed 895137 2>&1 | tee Pver_2_baypass_coremodel_defaultparm_log.txt")

system("g_baypass -gfile pver.popgenclust.baypass -outprefix Pver_3 -nthreads 10 -seed 59426 2>&1 | tee Pver_3_baypass_coremodel_defaultparm_log.txt")


#### 09b.02 Checking MCMC convergence on three different runs of BAYPASS under the core model
source("/home/buitracn/RADseq/Genotyping/tools/baypass_2.2/utils/baypass_utils.R")
library(geigen)

# Compare Omega covariance matrix from three differnt runs (using 3 different random seed)
pver_omega1 <- as.matrix(read.table("Pver_1_mat_omega.out"))
pver_omega2 <- as.matrix(read.table("Pver_2_mat_omega.out"))
pver_omega3 <- as.matrix(read.table("Pver_3_mat_omega.out"))

fmd.dist(pver_omega1,pver_omega2) #0.008613247
fmd.dist(pver_omega1,pver_omega3) #0.002696452
fmd.dist(pver_omega2,pver_omega3) #0.006676777

# Check differences in candidate loci for selection (done in terminal) - $7 = log10(1/pvalue) p.value(0.05)
# awk '{if ($7 >= "1.32") print $1}' Pver_1_summary_pi_xtx.out | sort  > candidate_loci_pver1
# awk '{if ($7 >= "1.32") print $1}' Pver_2_summary_pi_xtx.out | sort  > candidate_loci_pver2
# awk '{if ($7 >= "1.32") print $1}' Pver_3_summary_pi_xtx.out | sort  > candidate_loci_pver3

# wc -l candidate_loci_pver*
# 2234 candidate_loci_pver1
# 2214 candidate_loci_pver2
# 2234 candidate_loci_pver3

# How many loci where detected in common among those selected as candidates under selection after each MCMC run
# comm candidate_loci_pver1 candidate_loci_pver2 | awk -F"\t" '{print $3}' | sed '/^$/d' | wc -l #2110
# comm candidate_loci_pver1 candidate_loci_pver3 | awk -F"\t" '{print $3}' | sed '/^$/d' | wc -l #2126
# comm candidate_loci_pver2 candidate_loci_pver3 | awk -F"\t" '{print $3}' | sed '/^$/d' | wc -l #2112


#### 09b.03 POD calibration
library(geigen)
library(corrplot)
library(wesanderson)
library(ggplot2)
library(dplyr)
source("/home/buitracn/RADseq/Genotyping/tools/baypass_2.2/utils/baypass_utils.R")
source("/home/buitracn/RADseq/Genotyping/tools/manhattan.plot.function.qman.R") # I've manipulate the original qqman manhattan function to change the line types and colors

# set the working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_positive_selection_analyses_20210308/BAYPASS")


### Omega matrix  Pver3 ##
pver_omega3 <- as.matrix(read.table("Pver_3_mat_omega.out"))

pop.names <- c("PCL1", "PCL2")
dimnames(pver_omega3) <- list(pop.names,pop.names)

### Omega matrix Visualization 
# using SVD decomposition
pdf("./pver_SVD_Omega_matrix_defaultBaypass.pdf", width = 6, height = 6, pointsize = 20)
plot.omega(omega=pver_omega3,pop.names=pop.names, pos = 1, col = c("#9370DB", "#00CED1")) # this corresponds to the colors assigned by pong to the major mode selected
dev.off()


# as a correlation plot
cor.mat <- cov2cor(pver_omega3)
pal <- wes_palette("Zissou1", 10, type = "continuous")
pdf("./pver_Correlation_Omega_matrix_defaultBaypass.pdf", width = 6, height = 6)
corrplot(cor.mat,method="color", col=pal, mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)), tl.cex = 0.7, tl.col = "black")
dev.off()

# as a heatmap and hierarchical clustering tree (using the average agglomeration method)
pdf("./pver_Heatmap_ordered_Omega_matrix_defaultBaypass.pdf", width = 6, height = 6)
hclust.ave <- function(x) hclust(x, method="average")
heatmap(1-cor.mat,hclustfun = hclust.ave,
        main = expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"),
        col=pal)
dev.off()


### XtX statistics (Fst analog)
# Gunther and Coop calculate an FST analog called XtX, which has been standardized by the covariance among populations.
# Calibrating XtX statistics (FST-like) with the simulation and analysis of PODs
# pver_3 run

# XtX stats 
pver.snp.xtx <- read.table("Pver_3_summary_pi_xtx.out",h=T)

#get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef <- read.table("Pver_3_summary_beta_params.out",h=T)$Mean
#upload the original data to obtain total allele count
pver.data <- geno2YN("pver.popgenclust.baypass")


### Simulate POD data to calibrate pver3 run
# Create the POD #better to simulate a higher nuber of SNPs e.g. 100,000
simu.pver <- simulate.baypass(omega.mat=pver_omega3,nsnp=100000,sample.size=pver.data$NN,
                              beta.pi=pi.beta.coef,pi.maf=0,suffix="pverpods") #pi.maf= 0 inactivates MAF filtering.

# Analyze the newly created POD (data file named“G.pverpods” in the example) giving another prefix for the output files:
system("/home/buitracn/RADseq/Genotyping/tools/baypass_2.2/sources/g_baypass -gfile G.pverpods -outprefix POD_pver3 -nthreads 10 -seed 31689 2>&1 | tee POD_pver_3_baypass_coremodel_defaultparm_log.txt")

### Sanity Check: Compare POD and original data estimates
# get estimate of omega from the POD analysis
pod.omega <- as.matrix(read.table("POD_pver3_mat_omega.out"))
plot(pod.omega,pver_omega3) ; abline(a=0,b=1)
fmd.dist(pod.omega,pver_omega3) # 0.9820601
#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution from the POD analysis
pod.pi.beta.coef <- read.table("POD_pver3_summary_beta_params.out",h=T)$Mean
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)
#plot the omega matrix based on POD (to be compared with the real data)
pdf("./PODpver_SVD_Omega_matrix_defaultBaypass.pdf", width = 6, height = 6)
plot.omega(omega=pod.omega,pop.names=pop.names, pos = 1, col = c("#9370DB", "#00CED1"))
dev.off()


## XtX calibration ##
# The XtX is indeed defined as the variance of the standardized allele frequencies of the SNP across the populations
# and is thus analogous to a SNP-specific FST that would account for the overall covariance structure of the population allele frequencies (Gunther and Coop 2013)

# get the pod XtX
pod.xtx <- read.table("POD_pver3_summary_pi_xtx.out",h=T)$M_XtX
# compute the 1% threshold
pod.thresh.sel <- quantile(pod.xtx,probs=0.999) 
# 99.9% 
# 5.515877
pod.thresh.bal <- quantile(pod.xtx,probs=0.001)
# 0.1% 
# 1.395626

# add the thresh to the actual XtX plot
plot(pver.snp.xtx$M_XtX)
abline(h=pod.thresh.sel,lty=2, col="red")

# plot the markers, then select the scaffolds that have candidate outliers for positive and balancing selection 
pver.snp.scaff <- read.table("../snps.id.dictionary.BAYPASS.txt", h=F)
colnames(pver.snp.scaff) <- c("Scaffold", "POS", "ID", "index")
pver.snp.xtx.scaff <- merge(pver.snp.xtx, pver.snp.scaff, by.x = "MRK", by.y = "index")

# Plot all the SNPs regardless on their position in the genome and add the thresholds for balancing and positive selection
ggplot(pver.snp.xtx.scaff, aes(x=MRK, y = M_XtX))+
  geom_point(alpha=0.5, fill="black", color="black", shape=21, size=3, stroke=0.3)+ 
  labs(y = "XtX") +
  labs(x = "SNPs") +
  geom_hline(yintercept = pod.thresh.sel, color = "dodgerblue2", linetype="dashed", size=0.8) + #threshold for positive selection
  geom_hline(yintercept = pod.thresh.bal, color = "darkorange1", linetype="dashed", size=0.8) + #threshold for balancing selection
  theme(
    axis.line = element_line(),
    #axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = 12, family = "Helvetica",face = "bold"),
    axis.text.y=element_text(size=11),
    legend.position = "none")

ggsave("pver_Baypass_XtX_POD_calibrated.pdf", width = 10, height = 6)

# How many SNPs are candidates for selection
nrow(pver.snp.xtx.scaff[pver.snp.xtx.scaff$M_XtX > pod.thresh.sel,]) # 85 candidates for positive selection
nrow(pver.snp.xtx.scaff[pver.snp.xtx.scaff$M_XtX < pod.thresh.bal,]) # 12 candidates for balancing selection

# Extract the scaffold IDs containing candidate SNPs for either positive or balancing selection
pver.scaff <- subset(pver.snp.xtx.scaff, M_XtX > pod.thresh.sel, select=Scaffold) # Only candidates for positive selection

# Extract information only for the scaffold IDs previously identified
pver.snp.scaff.sel <- merge(pver.snp.xtx.scaff, pver.scaff, by="Scaffold") # Only 85 Scaffolds are selected

pver.snp.scaff.sel$Scaffold <- as.numeric(gsub("Pver_xfSc", "", gsub("Pver_xpSc", "", pver.snp.scaff.sel$Scaffold))) # define the vector of scaffold IDs as numeric
pver.snp.scaff.sel$POS <- as.numeric(pver.snp.scaff.sel$POS)

write.table(pver.snp.scaff.sel, file = "BAYPASS.scaffol.with.pos.sel.SNPs.POD.calibrated.tsv", col.names = T, row.names = F, quote = F) 
# pver.snp.scaff.sel.sort <- pver.snp.scaff.sel[order(pver.snp.scaff.sel$Scaffold),]  
# scaff.sel.order <- sort(unique(pver.snp.scaff.sel.sort$Scaffold))
# pver.snp.scaff.sel.sort$Scaffold <- factor(pver.snp.scaff.sel$Scaffold, levels = scaff.sel.order, ordered = T)

# Plot SNPs per scaffold (only those scaffolds with candidate SNPs for selection)
# library(qqman)
pdf("./pver_Baypass_XtX_POD_calibrated_onlyCandidateSNPs.pdf", width = 10, height = 6)
manhattan(pver.snp.scaff.sel,chr="Scaffold", 
          p="M_XtX",bp="POS",snp="MRK",logp=FALSE, 
          ylab="XtX", xlab="Scaffold ID", ylim = c(0, 30), suggestiveline = pod.thresh.sel, genomewideline = F,
          col = alpha(c("grey47", "black"), 0.6))
dev.off()


#### 09b.04 POD Generate the list of markers under positive and balacing selection

# SNPs BAYPASS number identified as candidates for selection
pver.snps.positive.sel.numb <- pver.snp.xtx[pver.snp.xtx$M_XtX > pod.thresh.sel, 1]
pver.snps.balancing.sel.numb <- pver.snp.xtx[pver.snp.xtx$M_XtX < pod.thresh.bal, 1]

# SNPs ID under positive selection  
pver.snps.positive.sel.ID <- data.frame(pver.snp.scaff[pver.snps.positive.sel.numb, ])
write.table(pver.snps.positive.sel.ID, file = "BAYPASS.whitelist.markers.positive.selection.POD.calibrated.tsv", col.names = T, row.names = F, quote = F) 

# SNPs ID under balancing selection  
pver.snps.balancing.sel.ID <- data.frame(pver.snp.scaff[pver.snps.balancing.sel.numb, ])
write.table(pver.snps.balancing.sel.ID, file = "BAYPASS.whitelist.markers.balancing.selection.POD.calibrated.tsv", col.names = T, row.names = F, quote = F) 

