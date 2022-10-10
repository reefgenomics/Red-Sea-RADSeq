###################################################
##### RAD-seq - Red Sea - Pocilloporid corals #####
###################################################

##### 08a Linkage Disequilibrium decay visualization - Pocillopora verrucosa #####

library(binr)
library(ggplot2)
library(scales)

pver.PCL1.10Scaff <- read.delim("../pver.PCL1.LD.r2.plink.10largestScaff.ld", sep="", h=T)
pver.PCL2.10Scaff <- read.delim("../pver.PCL2.LD.r2.plink.10largestScaff.ld", sep="", h=T)
pver.mean.10Scaff <- rbind(pver.CL1.10Scaff, pver.CL2.10Scaff)

pver.clust.colors <- c("#00CED1", "#9370DB", "black") 
show_col(pver.clust.colors)

objects.10Scaff <- list(pver.PCL1.10Scaff, pver.PCL2.10Scaff, pver.mean.10Scaff)

pver.cluster <- c("pver.PCL1", "pver.PCL2", "pver.mean")

pver.10Scaff.objects.binr.list<-lapply(objects.10Scaff,function(object){
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
names(pver.10Scaff.objects.binr.list)<-pver.cluster
reg<-unlist(lapply(pver.10Scaff.objects.binr.list,nrow))
pver.10Scaff.objects.binr<-do.call(rbind,pver.10Scaff.objects.binr.list)
pver.10Scaff.objects.binr$gen.clust<-factor(sub('pver.','',rep(names(reg),each=reg[1])),levels=c("CL1", "CL2", "mean"))


ggplot(pver.10Scaff.objects.binr,aes(start,mean,color=gen.clust,group=gen.clust))+
  geom_point(alpha=0.5)+
  geom_line()+
  scale_color_manual(values = pver.clust.colors)+
  labs(x="Distance (kilobases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0, 2.5*10^4, 5*10^4, 7.5*10^4, 1*10^5),labels=c("0","25","50","75","100"))+
  ylim(0, 0.6)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("LD decay (mean r2) - Pocillopora verrucosa - Ten largest Scaffolds")

ggsave(filename = "pver.LDdecay.binr100.10largestScaffolds.pdf", width = 7, height = 4)


##### using evenly spaced bins of 100 bps #####

pver.10Scaff.objects.list<-lapply(objects.10Scaff,function(object){
  object$dist <- object$BP_B-object$BP_A
  object.summ <- data.frame(object$dist, object$R2) 
  colnames(object.summ) <- c("dist", "rsq")
  object.100Kb <- subset(object.summ, dist <= 100000)
  #object.100Kb$distc <- cut(object.100Kb$dist,breaks=c(seq(from=min(0),to=max(1000),by=10), seq(from=min(1001),to=max(100000),by=100)),include.lowest = T, right = T, ordered_result = T)
  object.100Kb$distc <- cut(object.100Kb$dist,breaks=c(seq(from=min(0),to=max(100000),by=100)),include.lowest = T, right = T, ordered_result = T)
  object.100Kb.bin <- object.100Kb %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
  object.100Kb.bin <- object.100Kb.bin %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                                end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                                mid=start+((end-start)/2))
  object.100Kb.bin$start <- as.numeric(object.100Kb.bin$start)
  return(object.100Kb.bin)
})

names(pver.10Scaff.objects.list)<-pver.cluster
reg<-unlist(lapply(pver.10Scaff.objects.list,nrow))
pver.10Scaff.objects.binr<-do.call(rbind,pver.10Scaff.objects.list)
pver.10Scaff.objects.binr$gen.clust<-factor(sub('pver.','',rep(names(reg),each=reg[1])),levels=c("CL1", "CL2", "mean"))


ggplot(pver.10Scaff.objects.binr,aes(start,mean,color=gen.clust,group=gen.clust))+
  geom_point(alpha=0.3)+
  geom_line(size=0.1)+
  scale_color_manual(values = pver.clust.colors)+
  labs(x="Distance (kilobases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0, 2.5*10^4, 5*10^4, 7.5*10^4, 1*10^5),labels=c("0","25","50","75","100"))+
  scale_y_continuous(breaks = round(seq(min(0), max(0.5), by = 0.05),2))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("LD decay (mean r2) - Pocillopora verrucosa - Ten largest Scaffolds")

ggsave(filename = "pver.LDdecay.regularbins100.10largestScaffolds.pdf", width = 5, height = 5)

