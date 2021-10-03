#Required-----
getwd()
setwd("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics")
library(ggplot2)
library(RColorBrewer)
library(aplot)
display.brewer.all()
col <- brewer.pal(n = 8, name = "Set2")[c(1,4,6)]
#Taxon_Ploidy----
AABB <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/02_Site_TaxaQC/chr001_AABB_taxaQCfile.txt", header= T, stringsAsFactors = F)
AABBDD <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/02_Site_TaxaQC/chr001_AABBDD_taxaQCfile.txt", header=T, stringsAsFactors = F)
DD <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/02_Site_TaxaQC/chr005_DD_taxaQCfile.txt", header= T, stringsAsFactors = F)
DD$ploidy <- "DD"
AABBDD$ploidy <- "AABBDD"
AABB$ploidy <- "AABB"
#taxa_Heterozygous_Proportion箱线图
pdata <- rbind(DD[,c(2,4)],AABBDD[,c(2,4)],AABB[,c(2,4)])
p <- ggplot(pdata,aes(x=ploidy,y=HeterozygousProportion,fill=ploidy))+geom_boxplot()
p+xlab("") + ylab("Heterozygous Proportion by Taxon") + labs(fill="ploidy") + theme_classic()+ scale_fill_brewer(palette = "Set2")
#taxa_Heterozygous_Proportion柱状图
pdata2<- pdata[which(pdata$MissRate<=0.21),]#过滤缺失率: 1546->1526
p <- ggplot(pdata,aes(x=HeterozygousProportion,fill=ploidy))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~ploidy)
p+xlab("Heterozygous Proportion by Taxon") + ylab("Density") + labs(fill="ploidy") + theme_classic()+ scale_fill_brewer(palette = "Set2")+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
#taxa_missing_rate柱状图
pdata <- rbind(DD[,c(3,4)],AABBDD[,c(3,4)],AABB[,c(3,4)])
pdata2<- pdata[which(pdata$MissRate<=0.21),]#过滤缺失率: 1546->1526
p <- ggplot(pdata2,aes(x=MissRate,fill=ploidy))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~ploidy)
p+xlab("Missing Rate by Taxon") + ylab("Density") + labs(fill="ploidy") + theme_classic()+ scale_fill_brewer(palette = "Set2")+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
#Site_Ploidy----
AABB <- read.table("chr001_AABB_siteQC_50000.txt", header= T, stringsAsFactors = F)
AABBDD <- read.table("chr001_AABBDD_siteQC_50000.txt", header=T, stringsAsFactors = F)
DD <- read.table("chr005_DD_siteQC_50000.txt", header=T, stringsAsFactors = F)
DD$ploidy <- "DD"
AABBDD$ploidy <- "AABBDD"
AABB$ploidy <- "AABB"
#Site_Heterozygous_Proportion柱状图
pdata2<- pdata[which(pdata$MissRate<=0.21),]#过滤缺失率: 1546->1526
p <- ggplot(pdata,aes(x=HeterozygousProportion,fill=ploidy))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~ploidy)
p+xlab("Heterozygous Proportion by Site") + ylab("Density") + labs(fill="ploidy") + theme_classic()+ scale_fill_brewer(palette = "Set2")+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
#site_missing_rate柱状图
pdata <- rbind(DD[,c(4,6)],AABBDD[,c(4,6)],AABB[,c(4,6)])
pdata2<- pdata[which(pdata$MissingRate<=0.21),]#过滤缺失率: 1546->1526
p <- ggplot(pdata,aes(x=MissingRate,fill=ploidy))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~ploidy)
p+xlab("Missing Rate by Site") + ylab("Density") + labs(fill="ploidy") + theme_classic()+ scale_fill_brewer(palette = "Set2")+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
#site_depth二维散点图_Lineage----
AA <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/chr1A_3k_depth.txt", header=F, stringsAsFactors = F)
BB <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/chr1B_3.1k_depth.txt", header=F, stringsAsFactors = F)
DD <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/chr1D_2.3k_depth.txt", header=F, stringsAsFactors = F)
DD$lineage <- "D"
AA$lineage <- "A"
BB$lineage <- "B"
colnames(AA)[1:2] <- c("Depth","SD")
colnames(BB)[1:2] <- c("Depth","SD")
colnames(DD)[1:2] <- c("Depth","SD")
pdata <- rbind(AA[,c(1,2,3)],BB[,c(1,2,3)],DD[,c(1,2,3)])
#quantile(pdata$Depth,c(0.1,0.9))
#10%     90% 
#3.91843 9.36970 
#quantile(pdata$SD,c(0.1,0.9))
#10%     90% 
#3.86971 6.15844 
Ainfo <- pdata[which(pdata$lineage=="A"),]
Binfo <- pdata[which(pdata$lineage=="B"),]
Dinfo <- pdata[which(pdata$lineage=="D"),]

pAd <- ecdf(Ainfo$Depth)
#> pAd(3.8)
#[1] 0.1016667
#> pAd(9.25)
#[1] 0.9023333
pAs <- ecdf(Ainfo$SD)
#> pAs(3.7)
#[1] 0.105
#> pAs(6)
#[1] 0.91

pBd <- ecdf(Binfo$Depth)
#> pBd(3.8)
#[1] 0.1035484
#> pBd(9.2)
#[1] 0.9045161
pBs <- ecdf(Binfo$SD)
#> pBs(4)
#[1] 0.09935484
#> pBs(6.2)
#[1] 0.8990323

pDd <- ecdf(Dinfo$Depth)
#> pDd(4.7)
#[1] 0.103913
#> pDd(9.7)
#[1] 0.9052174
pDs <- ecdf(Dinfo$SD)
#> pDs(4.1)
#[1] 0.09956522
#> pDs(6.4)
#[1] 0.9017391

p1<-ggplot(pdata,aes(x = Depth, y = SD, color=lineage))+
  geom_point(alpha= 0.3)+
  scale_color_brewer(palette = "Set2")+
  #scale_color_manual(values=c("green","blue","grey"))+
  theme_bw()
  #geom_hline(yintercept = 3.86971,lty="dashed")+
  #geom_hline(yintercept = 6.15844,lty="dashed")+
  #geom_vline(xintercept = 3.91843,lty="dashed")+
  #geom_vline(xintercept = 9.36970,lty="dashed")

p2<-ggplot(pdata,aes(Depth))+
  geom_density(fill="grey",alpha=0.5)+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p3<-ggplot(pdata,aes(SD))+
  geom_density(fill="grey",alpha=0.5)+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  coord_flip()
#累计分布函数
p4<-ggplot(pdata, aes(Depth, colour = lineage)) + 
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_color_brewer(palette = "Set2")+
  stat_ecdf()+
  scale_y_reverse()
p5<-ggplot(pdata, aes(SD, colour = lineage)) + 
  scale_color_brewer(palette = "Set2")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  stat_ecdf()+
  coord_flip()+
  scale_y_reverse()
p1%>%
  insert_top(p2,height = 0.5)%>%
  insert_right(p3,0.5)%>%
  insert_bottom(p4,0.5)%>%
  insert_left(p5,0.5)

#Taxon_Lineage----
AA <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/chr1A_taxaQCfile.txt", header= T, stringsAsFactors = F)
BB <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/chr1B_taxaQCfile.txt", header=T, stringsAsFactors = F)
DD <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/chr1D_taxaQCfile.txt", header= T, stringsAsFactors = F)
DD$lineage <- "DD"
BB$lineage <- "BB"
AA$lineage <- "AA"
#taxa_Heterozygous_Proportion箱线图
pdata <- rbind(DD[,c(2,4)],AA[,c(2,4)],BB[,c(2,4)])
p <- ggplot(pdata,aes(x=lineage,y=HeterozygousProportion,fill=lineage))+geom_boxplot()
p+xlab("") + ylab("Heterozygous Proportion by Taxon") + labs(fill="lineage") + theme_classic()+ scale_fill_brewer(palette = "Accent")+theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))
#taxa_Heterozygous_Proportion柱状图
pdata2<- pdata[which(pdata$MissRate<=0.21),]#过滤缺失率: 1546->1526
p <- ggplot(pdata,aes(x=HeterozygousProportion,fill=lineage))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~lineage)
p+xlab("Heterozygous Proportion by Taxon") + ylab("Density") + labs(fill="lineage") + theme_classic()+ scale_fill_brewer(palette = "Accent")+theme(axis.text.x = element_text(size=15, angle = 45, hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size=15),axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))
#taxa_missing_rate柱状图
pdata <- rbind(DD[,c(3,4)],BB[,c(3,4)],AA[,c(3,4)])
pdata2<- pdata[which(pdata$MissRate<=0.21),]#过滤缺失率: 1546->1526
p <- ggplot(pdata2,aes(x=MissRate,fill=lineage))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~lineage)
p+xlab("Missing Rate by Taxon") + ylab("Density") + labs(fill="lineage") + theme_classic()+ scale_fill_brewer(palette = "Accent")+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size=15),axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))
#Site_Lineage----
AA <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/chr1A_siteQCfile_10k.txt", header= T, stringsAsFactors = F)
BB <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/chr1B_siteQCfile_10k.txt", header=T, stringsAsFactors = F)
DD <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/chr1D_siteQCfile_10k.txt", header=T, stringsAsFactors = F)
DD$lineage <- "DD"
BB$lineage <- "BB"
AA$lineage <- "AA"
#Site_Heterozygous_Proportion柱状图
pdata <- rbind(DD[,c(3,6)],BB[,c(3,6)],AA[,c(3,6)])
p <- ggplot(pdata,aes(x=HeterozygousProportion,fill=lineage))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~lineage)
p+xlab("Heterozygous Proportion by Site") + ylab("Density") + labs(fill="lineage") + theme_classic()+ scale_fill_brewer(palette = "Accent")+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size=15),axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))
#site_missing_rate柱状图
pdata <- rbind(DD[,c(4,6)],AA[,c(4,6)],BB[,c(4,6)])
pdata2<- pdata[which(pdata$MissingRate<=0.21),]#过滤缺失率: 1546->1526
p <- ggplot(pdata,aes(x=MissingRate,fill=lineage))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~lineage)
p+xlab("Missing Rate by Site") + ylab("Density") + labs(fill="lineage") + theme_classic()+ scale_fill_brewer(palette = "Accent")+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size=15),axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))
#MAF distribution_Lineage----
pdata <- rbind(AA[,c(5,6)],BB[,c(5,6)],DD[,c(5,6)])
p <- ggplot(pdata,aes(x=Maf,fill=lineage))+geom_histogram(aes(y = ..density..),bins = 25)+ facet_grid(.~lineage)
p+xlab("Maf") + ylab("Density") + labs(fill="lineage") + theme_classic()+ scale_fill_brewer(palette = "Accent")+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size=15),axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))

pdata$class <- NA
pdata[which(pdata$Maf<=0.05),3] <- "Maf<=0.05"
pdata[which(pdata$Maf>0.05),3] <- "Maf>0.05"
#maf 1
p <- ggplot(pdata,aes(x=Maf)) +
  geom_histogram(aes(y = ..density..),bins = 25) + 
  facet_grid(.~linaege)
p+
  xlab("Maf") + 
  ylab("Proportion") + 
  labs(fill="ploidy") + 
  theme_classic() + 
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
#maf 2
p <- ggplot(pdata,aes(x=class,fill=lineage)) +
  geom_histogram(aes(y = ..density..)) +
  geom_bar(aes(fill=factor(cyl)),position="fill")
p <- ggplot(pdata,aes(x=lineage)) +
  geom_bar(aes(fill=factor(class)),position="fill")
p +
  xlab("lineage") + 
  ylab("Proportion") + 
  labs(fill="") + 
  theme_classic() + 
  scale_fill_brewer(palette = "Accent")+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),axis.title.y = element_text(size=15),axis.title.x = element_text(size=15))
