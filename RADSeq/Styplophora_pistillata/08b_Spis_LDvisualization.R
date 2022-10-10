###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 08a Linkage Disequilibrium decay visualization - Stylophora pistillata #####

# Libraries
library(binr)
library(ggplot2)
library(scales)

spis.SCL1.10Scaff <- read.delim("spis.SCL1.LD.r2.plink.10largestScaff.ld", sep="", h=T)
spis.SCL2.10Scaff <- read.delim("spis.SCL2.LD.r2.plink.10largestScaff.ld", sep="", h=T)
spis.SCL3.10Scaff <- read.delim("spis.SCL3.LD.r2.plink.10largestScaff.ld", sep="", h=T)
spis.SCL4.10Scaff <- read.delim("spis.SCL4.LD.r2.plink.10largestScaff.ld", sep="", h=T)
spis.SCL5.10Scaff <- read.delim("spis.SCL5.LD.r2.plink.10largestScaff.ld", sep="", h=T)
spis.SCL6.10Scaff <- read.delim("spis.SCL6.LD.r2.plink.10largestScaff.ld", sep="", h=T)
spis.mean.10Scaff <- rbind(spis.SCL1.10Scaff, spis.SCL2.10Scaff, spis.SCL3.10Scaff, spis.SCL4.10Scaff, spis.SCL5.10Scaff, spis.SCL6.10Scaff)


spis.mean.10Scaff <- rbind(spis.SCL1.10Scaff, spis.SCL2.10Scaff, spis.SCL3.10Scaff, spis.SCL4.10Scaff, spis.SCL5.10Scaff, spis.SCL6.10Scaff)

spis.clust.colors <- c("#FE92CD", "#1E90FF", "#FF4500", "#C71585", "#FFD700", "#32CD32" ,"black") # this corresponds to the colors assigned by pong to the major mode selected
show_col(spis.clust.colors)

objects.10Scaff <- list(spis.SCL1.10Scaff, spis.SCL2.10Scaff, spis.SCL3.10Scaff, spis.SCL4.10Scaff, spis.SCL5.10Scaff, spis.SCL6.10Scaff, spis.mean.10Scaff)


spis.cluster <- c("spis.CL1", "spis.CL2", "spis.CL3", "spis.CL4", "spis.CL5", "spis.CL6", "spis.mean")

spis.10Scaff.objects.binr.list<-lapply(objects.10Scaff ,function(object){
  object$dist <- object$BP_B-object$BP_A
  object.summ <- data.frame(object$dist, object$R2) 
  colnames(object.summ) <- c("dist", "rsq")
  object.100Kb <- subset(object.summ, dist <= 100000)
  #object.100Kb <- object.100Kb[!(object.100Kb$dist == 0),]
  bins <- bins(object.100Kb$dist, target.bins = 100, minpts = 100)
  object.100Kb$distc <- cut(object.100Kb$dist, bins.getvals(bins), labels = names(bins$binct))
  object.100Kb.bin <- object.100Kb %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
  object.100Kb.bin <- object.100Kb.bin %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                                end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                                mid=start+((end-start)/2))
  object.100Kb.bin$start <- as.numeric(object.100Kb.bin$start)
  return(object.100Kb.bin)
  
})
names(spis.10Scaff.objects.binr.list)<-spis.cluster
reg<-unlist(lapply(spis.10Scaff.objects.binr.list,nrow))
spis.10Scaff.objects.binr<-do.call(rbind,spis.10Scaff.objects.binr.list)
spis.10Scaff.objects.binr$gen.clust<-factor(sub('spis.','',rep(names(reg),each=reg[1])),levels=c("CL1", "CL2", "CL3", "CL4", "CL5", "CL6", "mean"))

ggplot(spis.10Scaff.objects.binr,aes(start,mean,color=gen.clust,group=gen.clust))+
  geom_point(alpha=0.5)+
  geom_line()+
  scale_color_manual(values = spis.clust.colors)+
  labs(x="Distance (kilobases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0, 2.5*10^4, 5*10^4, 7.5*10^4, 1*10^5),labels=c("0","25","50","75","100"))+
  ylim(0, 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("LD decay (mean r2) - Stylophora pistillata ")

ggsave(filename = "spis.LDdecay.binr100.K6.10largestScaffolds.pdf", width = 7, height = 4)


##### using evenly spaced bins of 100 bps #####

spis.10Scaff.objects.list<-lapply(objects.10Scaff ,function(object){
  object$dist <- object$BP_B-object$BP_A
  object.summ <- data.frame(object$dist, object$R2) 
  colnames(object.summ) <- c("dist", "rsq")
  object.100Kb <- subset(object.summ, dist <= 100000)
  #object.100Kb <- object.100Kb[!(object.100Kb$dist == 0),]
  object.100Kb$distc <- cut(object.100Kb$dist,breaks=c(seq(from=min(0),to=max(100000),by=100)),include.lowest = T, right = T, ordered_result = T)
  object.100Kb.bin <- object.100Kb %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
  object.100Kb.bin <- object.100Kb.bin %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                                end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                                mid=start+((end-start)/2))
  object.100Kb.bin$start <- as.numeric(object.100Kb.bin$start)
  return(object.100Kb.bin)
  
})
names(spis.10Scaff.objects.list)<-spis.cluster
reg<-unlist(lapply(spis.10Scaff.objects.list,nrow))
spis.10Scaff.objects.binr<-do.call(rbind,spis.10Scaff.objects.list)
spis.10Scaff.objects.binr$gen.clust<-factor(sub('spis.','',rep(names(reg),each=reg[1])),levels=c("CL1", "CL2", "CL3", "CL4", "CL5", "CL6", "mean"))

ggplot(spis.10Scaff.objects.binr,aes(start,mean,color=gen.clust,group=gen.clust))+
  geom_point(alpha=0.3)+
  geom_line(size=0.1)+
  scale_color_manual(values = spis.clust.colors)+
  labs(x="Distance (kilobases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0, 2.5*10^4, 5*10^4, 7.5*10^4, 1*10^5),labels=c("0","25","50","75","100"))+
  scale_y_continuous(breaks = round(seq(min(0), max(0.5), by = 0.05),2))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("LD decay (mean r2) - Stylophora pistillata - Ten largest Scaffolds ")

ggsave(filename = "spis.LDdecay.regularbins100.K6.10largestScaffolds.pdf", width = 5, height = 5)

