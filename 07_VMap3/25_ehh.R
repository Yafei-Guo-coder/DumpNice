################################################ hard versus soft figure 4 #############################################
library(ggplot2)
#library(gridExtra)
library(RColorBrewer)
#library(fitdistrplus)
library(ggrepel)

setwd("/Users/guoyafei/Documents/02_VmapIII/19_ehh")
name <- read.table("/Users/guoyafei/Documents/02_VmapIII/11_piratio/pi-ratio/gene.txt", header=F, stringsAsFactors = F)
rownames(name) <- name$V1

cultivar_hap75_up <- 146.13706477215
cultivar_hap75_down <- 130.43189412738
cultivar_hap100_up <- 139.44228330749
cultivar_hap100_down <- 129.97111189885

landrace_hap75_up <- 115.43296111507
landrace_hap75_down <- 102.02953568688
landrace_hap100_up <- 106.79791980148
landrace_hap100_down <- 103.28076413671

#####
#hap75_landrace
all1 <- read.table("landrace.AB.frq5.ehh.hap75", header=F,stringsAsFactors = F)
all1$V5 <- all1$V2/572
all1$name <- name[all1$V1,2]
all1$Frequency <- all1$V5*100
pdf("landrace_hap75_up.pdf",height=4.8, width=8.6)
ggplot(all1, aes(name, log(V3),size=Frequency)) +
  geom_point(colour = "orange",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_hline(yintercept=log(landrace_hap100_up),color='red',linetype = "dotted")+
  xlab("Genes")+
  ylab("log(iHH)")+
  theme_bw()+
  theme(legend.title=element_text(size = 12),axis.text.x = element_text(size = 10,angle=45,hjust=.5, vjust=.5), axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

all1 <- read.table("landrace.AB.frq5.ehh.hap75", header=F,stringsAsFactors = F)
all1$V5 <- all1$V2/572
all1$name <- name[all1$V1,2]
all1$Frequency <- all1$V5*100
pdf("landrace_hap75_down.pdf",height=4.8, width=8.6)
ggplot(all1, aes(name, log(V4),size=Frequency)) +
  geom_point(colour = "orange",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_hline(yintercept=log(landrace_hap100_down),color='red',linetype = "dotted")+
  xlab("Genes")+
  ylab("log(iHH)")+
  theme_bw()+
  theme(legend.title=element_text(size = 12),axis.text.x = element_text(size = 10,angle=45,hjust=.5, vjust=.5), axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()
hap75_landrace
#####
all1 <- read.table("landrace.AB.frq5.ehh.hap100", header=F,stringsAsFactors = F)
all1$V5 <- all1$V2/572
all1$name <- name[all1$V1,2]
all1$Frequency <- all1$V5*100
pdf("landrace_hap100_up.pdf",height=4.8, width=8.6)
ggplot(all1, aes(name, log(V3),size=Frequency)) +
  geom_point(colour = "orange",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_hline(yintercept=log(landrace_hap100_up),color='red',linetype = "dotted")+
  xlab("Genes")+
  ylab("log(iHH)")+
  theme_bw()+
  theme(legend.title=element_text(size = 12),axis.text.x = element_text(size = 10,angle=45,hjust=.5, vjust=.5), axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

all1 <- read.table("landrace.AB.frq5.ehh.hap100", header=F,stringsAsFactors = F)
all1$V5 <- all1$V2/572
all1$name <- name[all1$V1,2]
all1$Frequency <- all1$V5*100
pdf("landrace_hap100_down.pdf",height=4.8, width=8.6)
ggplot(all1, aes(name, log(V4),size=Frequency)) +
  geom_point(colour = "orange",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_hline(yintercept=log(landrace_hap100_down),color='red',linetype = "dotted")+
  xlab("Genes")+
  ylab("log(iHH)")+
  theme_bw()+
  theme(legend.title=element_text(size = 12),axis.text.x = element_text(size = 10,angle=45,hjust=.5, vjust=.5), axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

#####
all1 <- read.table("cultivar.AB.frq5.ehh.hap75", header=F,stringsAsFactors = F)
all1$V5 <- all1$V2/479
all1$name <- name[all1$V1,2]
all1 <- all1[!is.na(all1$name),]
all1$Frequency <- all1$V5*100
pdf("cultivar_hap75_up.pdf",height=4.8, width=8.6)
ggplot(all1, aes(name, log(V3),size=Frequency)) +
  geom_point(colour = "orange",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_hline(yintercept=log(landrace_hap100_up),color='red',linetype = "dotted")+
  xlab("Genes")+
  ylab("log(iHH)")+
  theme_bw()+
  theme(legend.title=element_text(size = 12),axis.text.x = element_text(size = 10,angle=45,hjust=.5, vjust=.5), axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

all1 <- read.table("cultivar.AB.frq5.ehh.hap75", header=F,stringsAsFactors = F)
all1$V5 <- all1$V2/479
all1$name <- name[all1$V1,2]
all1 <- all1[!is.na(all1$name),]
all1$Frequency <- all1$V5*100
pdf("cultivar_hap75_down.pdf",height=4.8, width=8.6)
ggplot(all1, aes(name, log(V4),size=Frequency)) +
  geom_point(colour = "orange",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_hline(yintercept=log(landrace_hap100_down),color='red',linetype = "dotted")+
  xlab("Genes")+
  ylab("log(iHH)")+
  theme_bw()+
  theme(legend.title=element_text(size = 12),axis.text.x = element_text(size = 10,angle=45,hjust=.5, vjust=.5), axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

#####
all1 <- read.table("cultivar.AB.frq5.ehh.hap100", header=F,stringsAsFactors = F)
all1$V5 <- all1$V2/479
all1$name <- name[all1$V1,2]
all1 <- all1[!is.na(all1$name),]
all1$Frequency <- all1$V5*100
pdf("cultivar_hap100_up.pdf",height=4.8, width=8.6)
ggplot(all1, aes(name, log(V3),size=Frequency)) +
  geom_point(colour = "orange",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_hline(yintercept=log(landrace_hap100_up),color='red',linetype = "dotted")+
  xlab("Genes")+
  ylab("log(iHH)")+
  theme_bw()+
  theme(legend.title=element_text(size = 12),axis.text.x = element_text(size = 10,angle=45,hjust=.5, vjust=.5), axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

all1 <- read.table("cultivar.AB.frq5.ehh.hap100", header=F,stringsAsFactors = F)
all1$V5 <- all1$V2/479
all1$name <- name[all1$V1,2]
all1 <- all1[!is.na(all1$name),]
all1$Frequency <- all1$V5*100
pdf("cultivar_hap100_down.pdf",height=4.8, width=8.6)
ggplot(all1, aes(name, log(V4),size=Frequency)) +
  geom_point(colour = "orange",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_hline(yintercept=log(landrace_hap100_down),color='red',linetype = "dotted")+
  xlab("Genes")+
  ylab("log(iHH)")+
  theme_bw()+
  theme(legend.title=element_text(size = 12),axis.text.x = element_text(size = 10,angle=45,hjust=.5, vjust=.5), axis.title.x = element_text(size = 12),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()
