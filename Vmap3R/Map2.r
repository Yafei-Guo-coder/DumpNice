#library(ggplot2)
library(ggmap)
#library(sp)
#library(maptools)
#library(maps)
#library(psych)
taxa <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/taxa_location.txt", header=T, stringsAsFactors = F)
taxa_AB <- taxa[which(taxa$Ploidy=="AABB"),]
taxa_ABD <- taxa[which(taxa$Ploidy=="AABBDD"),]
taxa_D <- taxa[which(taxa$Ploidy=="DD"),]

#绘图
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="white") 
mp<-ggplot()+mapworld+ylim(-90,90)+theme_classic()

#AABBDD----
mp2<-mp+geom_point(aes(x=taxa_ABD$Longitude, y=taxa_ABD$Latitude, color=taxa_ABD$TreeValidatedGroupbySubspecies),size=1.5)+scale_size(range=c(1,1))+ theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
 #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
 #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())

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



