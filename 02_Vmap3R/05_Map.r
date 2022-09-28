#library(ggplot2)
library(ggmap)
library(RColorBrewer)
#library(sp)
#library(maptools)
#library(maps)
#library(psych)
#画样本分布图----
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/13_Plots/09_Map/")
taxa <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/13_Plots/09_Map/795Taxa.txt", header=T, sep="\t",stringsAsFactors = F)
taxa_AB <- taxa[which(taxa$Ploidy=="AABB"),]
taxa_ABD <- taxa[which(taxa$Ploidy=="AABBDD"),]
taxa_D <- taxa[which(taxa$Ploidy=="DD"),]

#VMap3
setwd("/Users/guoyafei/Documents/02_VmapIII/01_表格")
taxa <- read.table("location.txt", header=T, sep="\t",stringsAsFactors = F)

#绘图
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70") 
mp<-ggplot()+mapworld+ylim(-60,90) +theme_classic()
color <- brewer.pal(8, "Dark2")[c(1,2,3,4,6)]
#AABBDD----
mp2<-mp+geom_point(aes(x=taxa$Longitude, y=taxa$Latitude, color=taxa$RDA_Region),size=2,alpha=0.7)+
  scale_size(range=c(1,1))+ 
  scale_color_manual(values = color) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
 #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
 #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) 
  #theme(panel.grid =element_blank()) + 
  #theme(axis.ticks.y = element_blank())

#AABB----
mp2<-mp+geom_point(aes(x=taxa_AB$Longitude, y=taxa_AB$Latitude, color=taxa_AB$TreeValidatedGroupbySubspecies),size=1.5)+scale_size(range=c(1,1))+ theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())

#DD----
mp2<-mp+geom_point(aes(x=taxa_D$Longitude, y=taxa_D$Latitude, color=taxa_D$TreeValidatedGroupbySubspecies),size=1.5)+scale_size(range=c(1,1))+ theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())
#计算各位点的点的数量并标记----
pos_count <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/pos_count.txt",header=T,stringsAsFactors = F)
mp2<-mp+geom_point(aes(x=pos_count$Longitude, y=pos_count$Latitude),size=2)+scale_size(range=c(1,1))+ theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())




#PCA和地理位置(纬度海拔）聚类画图AABBDD----
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="white") 
mp<-ggplot()+mapworld+ylim(-90,90)+theme_classic()
data1 <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/Out/TXT/452_ABD_19bio_纬度海拔clust.txt", header=T, stringsAsFactors = F)
data <- data1[!is.na(data1$elevation),]
mp2<-mp+geom_point(aes(x=data$lon, y=data$lat_prec, shape=as.factor(data$k_4),color=data$PC2),size=1.5)+scale_size(range=c(1,1))+ theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())
mp2

#PCA和地理位置(经纬度）聚类画图AABBDD----
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="white") 
mp<-ggplot()+mapworld+ylim(-90,90)+theme_classic()
data1 <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/Out/TXT/452_ABD_19bio_经纬度海拔clust.txt", header=T, stringsAsFactors = F)
data <- data1[!is.na(data1$elevation),]
mp2<-mp+geom_point(aes(x=data$lon, y=data$lat_prec, color = as.factor(data$k_6)),size=1.5) + scale_size(range=c(1,1))+ theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())
mp2
#PCA和地理位置(纬度海拔及19变量）聚类画图AABBDD----
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="white") 
mp<-ggplot()+mapworld+ylim(-90,90)+theme_classic()

data1 <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/Out/TXT/452_ABD_19bio_纬度海拔19bio_clust.txt", header=T, stringsAsFactors = F)
data <- data1[!is.na(data1$elevation),]
mp2<-mp+geom_point(aes(x=data$lon, y=data$lat_prec, color=as.factor(data$k_6)),size=1.5)+scale_size(range=c(1,1))+ theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())
mp2
