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
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_positive_selection_analyses_20210308/BAYPASS")

# Run the 3 different MCMC chains (using the core model) and compare the covariance matrix "omega" between runs to assure convergence
system("g_baypass -gfile spis.popgenclust.baypass -outprefix Spis_1 -nthreads 10 -seed 742128 2>&1 | tee Spis_1_baypass_coremodel_defaultparm_log.txt")

system("g_baypass -gfile spis.popgenclust.baypass -outprefix Spis_2 -nthreads 10 -seed 820586 2>&1 | tee Spis_2_baypass_coremodel_defaultparm_log.txt")

system("g_baypass -gfile spis.popgenclust.baypass -outprefix Spis_3 -nthreads 10 -seed 31689 2>&1 | tee Spis_3_baypass_coremodel_defaultparm_log.txt")


#### 09b.02 Checking MCMC convergence on three different runs of BAYPASS under the core model

# libraries
source("/home/buitracn/RADseq/Genotyping/tools/baypass_2.2/utils/baypass_utils.R")
library(geigen)

# Compare Omega covariance matrix from three differnt runs (using 3 different random seed)
spis_omega1 <- as.matrix(read.table("Spis_1_mat_omega.out"))
spis_omega2 <- as.matrix(read.table("Spis_2_mat_omega.out"))
spis_omega3 <- as.matrix(read.table("Spis_3_mat_omega.out"))

fmd.dist(spis_omega1,spis_omega2) #0.02999398
fmd.dist(spis_omega1,spis_omega3) #0.02917885
fmd.dist(spis_omega2,spis_omega3) #0.01638651

# Check differences in candidate loci for selection (done in terminal) - $7 = log10(1/pvalue) p.value(0.05)
# awk '{if ($7 >= "1.32") print $1}' Spis_1_summary_pi_xtx.out | sort  > candidate_loci_Spis1
# awk '{if ($7 >= "1.32") print $1}' Spis_2_summary_pi_xtx.out | sort  > candidate_loci_Spis2
# awk '{if ($7 >= "1.32") print $1}' Spis_3_summary_pi_xtx.out | sort  > candidate_loci_Spis3

# wc -l candidate_loci_Spis*
# 1357 candidate_loci_Spis1
# 1357 candidate_loci_Spis2
# 1366 candidate_loci_Spis3

# How many loci where detected in common among those selected as candidates under selection after each MCMC run
# comm candidate_loci_Spis1 candidate_loci_Spis2 | awk -F"\t" '{print $3}' | sed '/^$/d' | wc -l #1305
# comm candidate_loci_Spis1 candidate_loci_Spis3 | awk -F"\t" '{print $3}' | sed '/^$/d' | wc -l #1315
# comm candidate_loci_Spis2 candidate_loci_Spis3 | awk -F"\t" '{print $3}' | sed '/^$/d' | wc -l #1307


#### 09b.03 POD calibration

library(corrplot)
library(wesanderson)
library(ggplot2)
library(dplyr)
source("/home/buitracn/RADseq/Genotyping/tools/manhattan.plot.function.qman.R") # I've manipulate the original qqman manhattan function to change the line types and colors


### Omega matrix Spis3
Spis_omega3 <- as.matrix(read.table("Spis_3_mat_omega.out"))

pop.names <- c("SCL1", "SCL2", "SCL3", "SCL4", "SCL5", "SCL6")
dimnames(Spis_omega3) <- list(pop.names,pop.names)

### Omega matrix Visualization 
# using SVD decomposition
pdf("./Spis_SVD_Omega_matrix_defaultBaypass.pdf", width = 6, height = 6, pointsize = 20)
plot.omega(omega=Spis_omega3,pop.names=pop.names, pos = 1, col = c("#FE92CD", "#1E90FF", "#FF4500", "#C71585", "#FFD700", "#32CD32")) # this corresponds to the colors assigned by pong to the major mode selected
dev.off()

### XtX statistics (Fst analog)
# Gunther and Coop calculate an FST analog called XtX, which has been standardized by the covariance among populations.
# Calibrating XtX statistics (FST-like) with the simulation and analysis of PODs
# Spis_3 run

# XtX stats 
Spis.snp.xtx <- read.table("Spis_3_summary_pi_xtx.out",h=T)

# get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef <- read.table("Spis_3_summary_beta_params.out",h=T)$Mean
# upload the original data to obtain total allele count
Spis.data <- geno2YN("spis.popgenclust.baypass")


### Simulate POD data to calibrate Spis3 run
# Create the POD #better to simulate a higher nuber of SNPs e.g. 100,000
simu.Spis <- simulate.baypass(omega.mat=Spis_omega3,nsnp=100000,sample.size=Spis.data$NN,
                              beta.pi=pi.beta.coef,pi.maf=0,suffix="Spispods") #pi.maf= 0 inactivates MAF filtering.

# Analyze the newly created POD (data file named“G.Spispods” in the example) giving another prefix for the output files:
system("/home/buitracn/RADseq/Genotyping/tools/baypass_2.2/sources/g_baypass -gfile G.Spispods -outprefix POD_Spis3 -nthreads 10 -seed 31689 2>&1 | tee POD_Spis_3_baypass_coremodel_defaultparm_log.txt")


### Sanity Check: Compare POD and original data estimates
# get estimate of omega from the POD analysis
pod.omega <- as.matrix(read.table("POD_Spis3_mat_omega.out"))
plot(pod.omega,Spis_omega3) ; abline(a=0,b=1)
fmd.dist(pod.omega,Spis_omega3) # 0.5220297
# get estimates (post. mean) of both the a_pi and b_pi parameters of
# the Pi Beta distribution from the POD analysis
pod.pi.beta.coef <- read.table("POD_Spis3_summary_beta_params.out",h=T)$Mean
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)
# plot the omega matrix based on POD (to be compared with the real data)
pdf("./PODSpis_SVD_Omega_matrix_defaultBaypass.pdf", width = 6, height = 6)
plot.omega(omega=pod.omega,pop.names=pop.names, pos = 1, col = c("#FE92CD", "#1E90FF", "#FF4500", "#C71585", "#FFD700", "#32CD32"))
dev.off()


### XtX calibration 
# The XtX is indeed defined as the variance of the standardized allele frequencies of the SNP across the populations
# and is thus analogous to a SNP-specific FST that would account for the overall covariance structure of the population allele frequencies (Gunther and Coop 2013)

# get the pod XtX
pod.xtx <- read.table("POD_Spis3_summary_pi_xtx.out",h=T)$M_XtX
# compute the 1% threshold
pod.thresh.sel <- quantile(pod.xtx,probs=0.999) 
# 99.9% 
# 14.51123 
pod.thresh.bal <- quantile(pod.xtx,probs=0.001)
# 0.1% 
# 3.11101

# add the thresh to the actual XtX plot
plot(Spis.snp.xtx$M_XtX)
abline(h=pod.thresh.sel,lty=2, col="red")

# plot the markers, then select the scaffolds that have candidate outliers for positive and balancing selection 
Spis.snp.scaff <- read.table("../snps.id.dictionary.BAYPASS.txt", h=F)
colnames(Spis.snp.scaff) <- c("Scaffold", "POS", "ID", "index")
Spis.snp.xtx.scaff <- merge(Spis.snp.xtx, Spis.snp.scaff, by.x = "MRK", by.y = "index")

# Plot all the SNPs regardless on their position in the genome and add the thresholds for balancing and positive selection
ggplot(Spis.snp.xtx.scaff, aes(x=MRK, y = M_XtX))+
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

ggsave("Spis_Baypass_XtX_POD_calibrated.pdf", width = 10, height = 6)

# How many SNPs are candidates for selection
nrow(Spis.snp.xtx.scaff[Spis.snp.xtx.scaff$M_XtX > pod.thresh.sel,]) # 111 candidates for positive selection
nrow(Spis.snp.xtx.scaff[Spis.snp.xtx.scaff$M_XtX < pod.thresh.bal,]) # 159 candidates for balancing selection

# Extract the scaffold IDs containing candidate SNPs for either positive or balancing selection
# Spis.scaff <- subset(Spis.snp.xtx.scaff, M_XtX > 14.51123 | M_XtX < 3.11101, select=Scaffold) 
Spis.scaff <- subset(Spis.snp.xtx.scaff, M_XtX > pod.thresh.sel, select=Scaffold) # Only candidates for positive selection

# Extract information only for the scaffold IDs previously identified
Spis.snp.scaff.sel <- merge(Spis.snp.xtx.scaff, Spis.scaff, by="Scaffold") # Only 104 Scaffolds are selected

Spis.snp.scaff.sel$Scaffold <- as.numeric(Spis.snp.scaff.sel$Scaffold) # define the vector of scaffold IDs as numeric

Spis.snp.scaff.sel$Scaffold <- as.numeric(Spis.snp.scaff.sel$Scaffold)

write.table(Spis.snp.scaff.sel, file = "BAYPASS.scaffol.with.pos.sel.SNPs.POD.calibrated.tsv", col.names = T, row.names = F, quote = F) 
# Spis.snp.scaff.sel.sort <- Spis.snp.scaff.sel[order(Spis.snp.scaff.sel$Scaffold),]  
# scaff.sel.order <- sort(unique(Spis.snp.scaff.sel.sort$Scaffold))
# Spis.snp.scaff.sel.sort$Scaffold <- factor(Spis.snp.scaff.sel$Scaffold, levels = scaff.sel.order, ordered = T)

# Plot SNPs per scaffold (only those scaffolds with candidate SNPs for selection)
pdf("./Spis_Baypass_XtX_POD_calibrated_onlyCandidateSNPs.pdf", width = 10, height = 6)
manhattan(Spis.snp.scaff.sel,chr="Scaffold", 
          p="M_XtX",bp="POS",snp="MRK",logp=FALSE, 
          ylab="XtX", xlab="Scaffold ID", ylim = c(0, 30), suggestiveline = pod.thresh.sel, genomewideline = F,
          col = alpha(c("grey47", "black"), 0.6))
dev.off()


#### 09b.04 POD Generate the list of markers under positive and balacing selection

# SNPs BAYPASS number identified as candidates for selection
Spis.snps.positive.sel.numb <- Spis.snp.xtx[Spis.snp.xtx$M_XtX > pod.thresh.sel, 1]
Spis.snps.balancing.sel.numb <- Spis.snp.xtx[Spis.snp.xtx$M_XtX < pod.thresh.bal, 1]

# SNPs ID under positive selection  
Spis.snps.positive.sel.ID <- data.frame(Spis.snp.scaff[Spis.snps.positive.sel.numb, ])
write.table(Spis.snps.positive.sel.ID, file = "BAYPASS.whitelist.markers.positive.selection.POD.calibrated.tsv", col.names = T, row.names = F, quote = F) 

# SNPs ID under balancing selection  
Spis.snps.balancing.sel.ID <- data.frame(Spis.snp.scaff[Spis.snps.balancing.sel.numb, ])
write.table(Spis.snps.balancing.sel.ID, file = "BAYPASS.whitelist.markers.balancing.selection.POD.calibrated.tsv", col.names = T, row.names = F, quote = F) 

