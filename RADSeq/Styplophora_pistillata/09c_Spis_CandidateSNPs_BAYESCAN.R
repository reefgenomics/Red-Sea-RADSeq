###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 09c Candidate SNPs for Positive Selection with BAYESCAN  - Stylophora pistillata #####

#### 09c.01 BAYESCAN run

# BAYESCAN - aims at identifying candidate loci under natural selection from genetic data, using differences in allele frequencies between populations.
# Using a prior odds of 1 would mean that for every locus, one assume a priori that the model including selection is as likely as the neutral model. 
# This can lead to false positives when testing a large number of markers. For example with uninformative data, one would conclude that half of the markers show a signal in favour of selection. 
# Higher prior odds will tend to eliminate false positives, but at the cost of reducing the power to detect any marker under selection. 
# A value of 10 seems reasonable for the identification of candidate loci within a few hundreds of markers, whereas values up to 10 000 are generally used in the context of genome wide association studies with millions of SNPs when people want to identify only the top candidates.
# The default value for prior odds in the program is 10 (for every 10 neutral loci in the data set, odds are that 1 locus is under selection).

# To be run in the bash environment
# export PATH="/home/buitracn/RADseq/Genotyping/tools/BayeScan:$PATH"

# working directory
cd /home/buitracn/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_positive_selection_analyses_20210308/BAYESCAN

# prior odds 100
bayescan_2.1  spis.popgenclust.bayescan -d monomophicsnpsindex2remove.txt -n 5000 -thin 20 -nbp 50 -pilot 5000 -burn 50000 -pr_odds 100 -threads 40 -od pr_odds100

# prior odds 500
bayescan_2.1  spis.popgenclust.bayescan -d monomophicsnpsindex2remove.txt -n 5000 -thin 20 -nbp 50 -pilot 5000 -burn 50000 -pr_odds 500 -threads 40 -od pr_odds500

# prior odds 1000       
bayescan_2.1  spis.popgenclust.bayescan -d monomophicsnpsindex2remove.txt -n 5000 -thin 20 -nbp 50 -pilot 5000 -burn 50000 -pr_odds 1000 -threads 40 -od pr_odds1000             


#### 09c.02 Identify the candidates SNPs for positive selection at each run (R environment)

# Libraries
library(dplyr)
library(ggplot2)

# Working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_positive_selection_analyses_20210308/BAYESCAN")

# Bayescan SNPs dictionary
spis.bayescan.dictionary <- read.delim("../snps.id.dictionary.BAYESCAN.txt", header = F)
colnames(spis.bayescan.dictionary) <- c("Scaffold", "POS", "ID", "index")


#### pr_100 ####
spis.bayescan.pr100 <- suppressWarnings(readr::read_table2(
  file = "./pr_odds100/spis.popgenclust.baye_fst.txt",
  skip = 1,
  col_names = c("BAYESCAN_MARKERS", "POST_PROB", "LOG10_PO", "Q_VALUE", "ALPHA", "FST"),
  col_types = c("iddddd"))) %>%
  dplyr::mutate(
    Q_VALUE = dplyr::if_else(Q_VALUE <= 0.0001, 0.0001, Q_VALUE),
    Q_VALUE = round(Q_VALUE, 4),
    POST_PROB = round(POST_PROB, 4),
    LOG10_PO = round(LOG10_PO, 4),
    ALPHA = round(ALPHA, 4),
    FST = round(FST, 6),
    SELECTION = factor(
      dplyr::if_else(ALPHA >= 0 & Q_VALUE <= 0.05, "diversifying",
                     dplyr::if_else(ALPHA >= 0 & Q_VALUE > 0.05, "neutral", "balancing"))),
    LOG10_Q = log10(Q_VALUE)
  )

# get a list of the snp index with evidence of positive, neutral and balancing selection
pr100.positive <- spis.bayescan.pr100[spis.bayescan.pr100$SELECTION=="diversifying", 1] # 76
pr100.neutral <- spis.bayescan.pr100[spis.bayescan.pr100$SELECTION=="neutral", 1] 
pr100.balancing <- spis.bayescan.pr100[spis.bayescan.pr100$SELECTION=="balancing", 1]

# get the ID of the SNPs identified as positive selection
Spis.snps.positive.sel.ID.pr100 <- data.frame(spis.bayescan.dictionary[pr100.positive$BAYESCAN_MARKERS, ])
write.table(Spis.snps.positive.sel.ID.pr100, file = "spis.BAYESCAN.pr100.whitelist.markers.positive.selection.tsv", col.names = T, row.names = F, quote = F) 


ggplot(data = spis.bayescan.pr100, aes(x = LOG10_Q, y = FST)) +
  geom_point(aes(fill=SELECTION),  pch=21, alpha=0.25, size=3) +
  scale_fill_manual(name="Selection", values=c("white","red","orange")) +
  labs(x = "Log10(Q_VALUE)") +
  labs(y = "Fst") +
  geom_vline(xintercept = c(log10(0.05)),color = "red") +
  theme(
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.title = element_text(size = 12, family = "Helvetica",face = "bold"),
    axis.text=element_text(size=11),
  )
ggsave("Spis_Bayescan_prodds_100.pdf", width = 8, height = 5)


#### pr_500 ####
spis.bayescan.pr500 <- suppressWarnings(readr::read_table2(
  file = "./pr_odds500/spis.popgenclust.baye_fst.txt",
  skip = 1,
  col_names = c("BAYESCAN_MARKERS", "POST_PROB", "LOG10_PO", "Q_VALUE", "ALPHA", "FST"),
  col_types = c("iddddd"))) %>%
  dplyr::mutate(
    Q_VALUE = dplyr::if_else(Q_VALUE <= 0.0001, 0.0001, Q_VALUE),
    Q_VALUE = round(Q_VALUE, 4),
    POST_PROB = round(POST_PROB, 4),
    LOG10_PO = round(LOG10_PO, 4),
    ALPHA = round(ALPHA, 4),
    FST = round(FST, 6),
    SELECTION = factor(
      dplyr::if_else(ALPHA >= 0 & Q_VALUE <= 0.05, "diversifying",
                     dplyr::if_else(ALPHA >= 0 & Q_VALUE > 0.05, "neutral", "balancing"))),
    LOG10_Q = log10(Q_VALUE)
  )

# get a list of the snp index with evidence of positive, neutral and balancing selection
pr500.positive <- spis.bayescan.pr500[spis.bayescan.pr500$SELECTION=="diversifying", 1] # 40
pr500.neutral <- spis.bayescan.pr500[spis.bayescan.pr500$SELECTION=="neutral", 1] 
pr500.balancing <- spis.bayescan.pr500[spis.bayescan.pr500$SELECTION=="balancing", 1]

# get the ID of the SNPs identified as positive selection
Spis.snps.positive.sel.ID.pr500 <- data.frame(spis.bayescan.dictionary[pr500.positive$BAYESCAN_MARKERS, ])
write.table(Spis.snps.positive.sel.ID.pr500 , file = "spis.BAYESCAN.pr500.whitelist.markers.positive.selection.tsv", col.names = T, row.names = F, quote = F) 


ggplot(data = spis.bayescan.pr500, aes(x = LOG10_Q, y = FST)) +
  geom_point(aes(fill=SELECTION),  pch=21, alpha=0.25, size=3) +
  scale_fill_manual(name="Selection", values=c("white","red","orange")) +
  labs(x = "Log10(Q_VALUE)") +
  labs(y = "Fst") +
  geom_vline(xintercept = c(log10(0.05)),color = "red") +
  theme(
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.title = element_text(size = 12, family = "Helvetica",face = "bold"),
    axis.text=element_text(size=11),
  )
ggsave("Spis_Bayescan_prodds_500.pdf", width = 8, height = 5)


#### pr_1000 ####
spis.bayescan.pr1000 <- suppressWarnings(readr::read_table2(
  file = "./pr_odds1000/spis.popgenclust.baye_fst.txt",
  skip = 1,
  col_names = c("BAYESCAN_MARKERS", "POST_PROB", "LOG10_PO", "Q_VALUE", "ALPHA", "FST"),
  col_types = c("iddddd"))) %>%
  dplyr::mutate(
    Q_VALUE = dplyr::if_else(Q_VALUE <= 0.0001, 0.0001, Q_VALUE),
    Q_VALUE = round(Q_VALUE, 4),
    POST_PROB = round(POST_PROB, 4),
    LOG10_PO = round(LOG10_PO, 4),
    ALPHA = round(ALPHA, 4),
    FST = round(FST, 6),
    SELECTION = factor(
      dplyr::if_else(ALPHA >= 0 & Q_VALUE <= 0.05, "diversifying",
                     dplyr::if_else(ALPHA >= 0 & Q_VALUE > 0.05, "neutral", "balancing"))),
    LOG10_Q = log10(Q_VALUE)
  )

# get a list of the snp index with evidence of positive, neutral and balancing selection
pr1000.positive <- spis.bayescan.pr1000[spis.bayescan.pr1000$SELECTION=="diversifying", 1] # 30
pr1000.neutral <- spis.bayescan.pr1000[spis.bayescan.pr1000$SELECTION=="neutral", 1] 
pr1000.balancing <- spis.bayescan.pr1000[spis.bayescan.pr1000$SELECTION=="balancing", 1]

# get the ID of the SNPs identified as positive selection
Spis.snps.positive.sel.ID.pr1000 <- data.frame(spis.bayescan.dictionary[pr1000.positive$BAYESCAN_MARKERS, ])
write.table(Spis.snps.positive.sel.ID.pr1000 , file = "spis.BAYESCAN.pr1000.whitelist.markers.positive.selection.tsv", col.names = T, row.names = F, quote = F) 


ggplot(data = spis.bayescan.pr1000, aes(x = LOG10_Q, y = FST)) +
  geom_point(aes(fill=SELECTION),  pch=21, alpha=0.25, size=3) +
  scale_fill_manual(name="Selection", values=c("white","red","orange")) +
  labs(x = "Log10(Q_VALUE)") +
  labs(y = "Fst") +
  geom_vline(xintercept = c(log10(0.05)),color = "red") +
  theme(
    axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.title = element_text(size = 12, family = "Helvetica",face = "bold"),
    axis.text=element_text(size=11),
  )
ggsave("Spis_Bayescan_prodds_1000.pdf", width = 8, height = 5)


### check common SNPs between the different runs ###
length(intersect(intersect(pr100.positive$BAYESCAN_MARKERS,pr500.positive$BAYESCAN_MARKERS), pr1000.positive$BAYESCAN_MARKERS)) #30
length(intersect(pr100.positive$BAYESCAN_MARKERS,pr500.positive$BAYESCAN_MARKERS)) #40
length(intersect(pr100.positive$BAYESCAN_MARKERS,pr1000.positive$BAYESCAN_MARKERS)) #30
length(intersect(pr500.positive$BAYESCAN_MARKERS,pr1000.positive$BAYESCAN_MARKERS)) #30

