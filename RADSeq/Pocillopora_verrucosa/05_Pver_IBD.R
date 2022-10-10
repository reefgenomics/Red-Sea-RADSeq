###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 05 Isolation by distance - Pocillopora verrucosa #####
# Using set of unlinked SNPs (35,208 SNPs - 296 ind)

# Libraries
library(dartR)

# Working directory
setwd("~/RADseq/Genotyping/Pver/p1.mac4.r0.8.316samples/pver_IBD_20210228")

####### Isolation by distance (IBS) using DartR
# genind object was loaded from the AMOVA folder
pver.pop.coord <- data.frame( lat=c(rep(28.52616667, 12), rep(28.4268, 23), rep(26.18751389, 25), rep(26.16661111, 25), rep(26.24119444, 20), rep(23.94741667, 17), rep(23.95533333, 17), rep(22.31916667, 25), rep(22.06722222, 23), rep(22.51333333, 25), rep(19.63511111, 22), rep(19.61402778, 20), rep(19.66569167, 20), rep(16.52518056, 22)),
                              long= c(rep(34.80397, 12), rep(34.75171, 23), rep(36.34920, 25), rep(36.39228, 25), rep(36.44036, 20), rep(38.17550, 17), rep(38.20444, 17), rep(38.85444, 25), rep(38.76917, 23), rep(38.92139, 25), rep(40.57536, 22), rep(40.63819, 20), rep(40.62266, 20), rep(42.03253, 22)))

pver.genind@other$latlong<-pver.pop.coord # add coordinates to the genind file


pver.gl <- gi2gl(pver.genind) # convert from genind to genlight

# Mantel test between the two distance matrixes
pver.ibd <- gl.ibd(gl = pver.gl,
  #Dgeo = pver.reef.geodist, Distances were calculated by dartR itself (Coordinates transformed to Mercator (google) projection to calculate distances in meters)
  projected = FALSE,
  permutations = 999,
  plot = TRUE
)
write.table(as.matrix(pver.ibd$Dgen), file = "pver.FSTtransfrormed.matrix.dartR.IBD.txt", sep = "\t", quote = F)

write.table(as.matrix(pver.ibd$Dgeo), file = "pver.natural.log.geodist.matrix.dartR.IBD.txt", sep = "\t", quote = F)

# dimensions of figure (5.20 x 4.19)

