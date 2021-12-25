#working directory: /Users/guoyafei/Documents/07_105XD
setwd("/Users/guoyafei/Documents/07_105XD")
Asite <- read.table("/Users/guoyafei/Documents/07_105XD/A_siteQCfile_3000.txt",header=T,stringsAsFactors = F)
Bsite <- read.table("/Users/guoyafei/Documents/07_105XD/B_siteQCfile_3000.txt",header=T,stringsAsFactors = F)
Dsite <- read.table("/Users/guoyafei/Documents/07_105XD/D_siteQCfile_3000.txt",header=T,stringsAsFactors = F)
allsite <- rbind(Asite,Bsite,Dsite)
ggplot(allsite, aes(x = HeterozygousProportion)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  theme_classic()
ggplot(allsite, aes(x = MissingRate)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  theme_classic()
ggplot(allsite, aes(x = Maf)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  theme_classic()

Ataxa <- read.table("/Users/guoyafei/Documents/07_105XD/A_taxaQCfile.txt",header=T,stringsAsFactors = F)
Btaxa <- read.table("/Users/guoyafei/Documents/07_105XD/B_taxaQCfile.txt",header=T,stringsAsFactors = F)
Dtaxa <- read.table("/Users/guoyafei/Documents/07_105XD/D_taxaQCfile.txt",header=T,stringsAsFactors = F)
alltaxa <- cbind(Ataxa,Btaxa,Dtaxa)
colnames(alltaxa) <- c("Taxa1","HeterozygousProportion1","MissRate1","Taxa2","HeterozygousProportion2","MissRate2","Taxa3","HeterozygousProportion3","MissRate3")
alltaxa$Heterozygous_Proportion <- NA
het <- alltaxa[,c(2,5,8)]
mis <- alltaxa[,c(3,6,9)]
alltaxa$Heterozygous_Proportion <- apply(het,1,mean)
alltaxa$MissRate <- apply(mis,1,mean)
ggplot(alltaxa, aes(x = MissRate )) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  theme_classic()
ggplot(alltaxa, aes(x = Heterozygous_Proportion )) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  theme_classic()

ABmds <- read.table("/Users/guoyafei/Documents/07_105XD/Dlineage.mds",header=T,stringsAsFactors = F)
info <- read.table("/Users/guoyafei/Documents/07_105XD/105.taxa.info.txt",header=T,stringsAsFactors = F)
allmds <- merge(ABmds,info,by.x="FID",by.y="Taxa.vmap3.")

ggplot(allmds, aes(C1, C2, color=as.factor(allmds$Origin_city))) +
  geom_point(size=2) +
  #scale_color_manual(values = c("#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  scale_color_manual(values = brewer.pal(11, "Set3")[c(1:11)])+
  #scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11))+
  theme_classic()+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size = 10),
    axis.title.x=element_text(size = 10),
  )

library("corrplot")
Dibs <- read.table("/Users/guoyafei/Documents/07_105XD/Dlineage.all.ibs.txt", header=T, row.names = 1, stringsAsFactors = F)
ibs <- as.matrix(Dibs)
ibsB <- as.matrix(Bibs)
rownames(ibs) <- names
colnames(ibs) <- names
#A:0.5 B:0.5 D:0.7
pdf("Dibs.pdf", width=100, height=100)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4, addCoef.col = "grey", type = "lower",number.cex=6,number.digits=3,tl.cex=10,cl.cex=12, cl.lim = c(0, 0.7))
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4, type = "lower",tl.cex=5, cl.cex=5,cl.lim = c(0, 0.7))

dev.off()

library(ggplot2)
P2<- ggplot(sub, aes(x=Lineage, y=CS, fill=Fill),na.rm=TRUE) + 
  #geom_violin() + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  scale_fill_manual(values = c("#6000b4", "#1f97ff", "#ffce6e"))+
  theme_bw()+
  ylab("IBS")+xlab("Lineages")+scale_y_continuous(limits=c(0,0.4)) + 
  theme(axis.text.x=element_text(angle=80, hjust = 1, colour="black", family="Times", size=200),
        axis.text.y=element_text(family="Times",size=300,face="plain"),
        axis.title.y=element_text(family="Times",size = 300,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())

pdf(file = "B_IBS_CS.pdf",width=200,height=100) #结果保存
print(P2)
dev.off()

Ald <- read.table("Alineage_LD_sub.ld",header=T,stringsAsFactors = F)
Ald$dis <- abs(Ald$BP_B-Ald$BP_A)
ggplot(Ald, aes(dis/100,R2)) +
  geom_point() +
  #scale_color_manual(values = c("#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  scale_color_manual(values = brewer.pal(11, "Set3")[c(1:11)])+
  #scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11))+
  theme_classic()+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size = 10),
    axis.title.x=element_text(size = 10),
  )


