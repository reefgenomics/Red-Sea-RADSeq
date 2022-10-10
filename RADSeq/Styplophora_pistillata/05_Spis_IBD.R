###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 05 Isolation by distance - Stylophora pistillata #####
# Using set of unlinked SNPs (25,318 SNPs - 367 ind)

# Libraries
library(dartR)

# Working directory
setwd("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_IBD_20210301")

####### Isolation by distance (IBS) using DartR
# genind object was loaded from the AMOVA folder
spis.genind <- load("~/RADseq/Genotyping/Spis/p1_mac4_r80_448samples/spis_AMOVA_20210226/spis.genind.RData")

spis.pop.coord <- data.frame( lat=c(rep(28.52616667, 22), rep(28.4268, 22), rep(26.18751389, 12), rep(26.24119444, 23), rep(26.18505556, 25), rep(23.94741667, 19), rep(23.95533333, 25), rep(23.9115, 25), rep(22.31916667, 8), rep(22.06722222, 18), rep(22.51333333, 20), rep(19.63511111, 14), rep(19.61402778, 21), rep(19.66569167, 25), rep(16.57930556, 25), rep(16.57899444, 19), rep(16.52518056, 18), rep(16.52736111, 26) ),
                              long= c(rep(34.80397222, 22), rep(34.75171389, 22), rep(36.3492, 12), rep(36.44036111, 23), rep(36.38302778, 25), rep(38.1755, 19), rep(38.20444444, 25), rep(38.15233333, 25), rep(38.85444444, 8), rep(38.76916667, 18), rep(38.92138889, 20), rep(40.57536111, 14), rep(40.63819444, 21), rep(40.62266389, 25), rep(42.14930556, 25), rep(42.23651944, 19), rep(42.03253056, 18), rep(42.03191667, 26) ))


spis.genind@other$latlong<-spis.pop.coord # add coordinates to the genind file

spis.gl <- gi2gl(spis.genind) # convert from genind to genlight

# Mantel test between the two distance matrixes
spis.ibd <- gl.ibd(gl = spis.gl,
                   projected = FALSE,
                   permutations = 999,
                   plot = TRUE
)
write.table(as.matrix(spis.ibd$Dgen), file = "spis.FSTtransfrormed.matrix.dartR.IBD.txt", sep = "\t", quote = F)

write.table(as.matrix(spis.ibd$Dgeo), file = "spis.natural.log.geodist.matrix.dartR.IBD.txt", sep = "\t", quote = F)

# dimensions of figure (5.20 x 4.19)


###########
## Test IBD by removing WAJ-R1
removeInd <- c("SWAJ-R1-23", "SWAJ-R1-25", "SWAJ-R1-26", "SWAJ-R1-27", "SWAJ-R1-29", "SWAJ-R1-30", "SWAJ-R1-33", "SWAJ-R1-34", "SWAJ-R1-38", "SWAJ-R1-42","SWAJ-R1-43", "SWAJ-R1-45")

spis.gi.no.waj.r1 <- spis.genind[!row.names(spis.genind@tab) %in% removeInd]
spis.gl.no.waj.r1 <- gi2gl(spis.gi.no.waj.r1)

# Mantel
spis.ibd.no.waj.r1 <- gl.ibd(gl = spis.gl.no.waj.r1,
                   projected = FALSE,
                   permutations = 999,
                   plot = TRUE
)

