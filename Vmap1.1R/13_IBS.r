setwd("/Users/guoyafei/Documents/个人项目/Project-2-Migration/migration/IBS_spelt/")
#Xinjiang_wheat
data <- read.table("B_xinjiang.txt",header=T,stringsAsFactors = F)
#Macha
data <- read.table("Macha.txt",header=T,stringsAsFactors = F)
#Persian
data <- read.table("Persian.txt",header=T,stringsAsFactors = F)

library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)
library(psych)

#Rivet_wheat
R <- data[which(data$Type=="Rivet_wheat" | data$Type=="Persian" ),]
#Polish_wheat
P <- data[which(data$Type=="Polish_wheat" | data$Type=="Xinjiang"),]
#Khorasan_wheat
K <- data[which(data$Type=="Khorasan_wheat"),]
#Durum
D <- data[which(data$Type=="Wild_emmer" | data$Type=="Macha"),]
D <- data[which(data$Type=="Domesticated_emmer" ),]

#Landrace
L <- data[which(data$Value=="5"),]

P$Value <- as.factor(P$Value)
R$Value <- as.factor(R$Value)
data$Value <- as.factor(data$Value)
D$Value <- as.factor(D$Value)
ggplot(L, aes(Logititude, Latitude))+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(color=mean),size=1.5)+
  #geom_point(aes(color=mean,shape=Value),size=3)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))

#画图文件：yafei@203
#/data2/yafei/Project3/CS_Vmap/A_ibs.txt
#/data2/yafei/Project3/CS_Vmap/B_ibs.txt
#/data2/yafei/Project3/CS_Vmap/D_ibs.txt

#工作目录：yafei@203:/data2/yafei/Project3/CS_Vmap

#热图  
library("corrplot")
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/13_Plots/03_IBS_heatmap")
#A lineage
ibs <- read.table("A_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

#A:0.5 B:0.5 D:0.7
pdf("IBS_A_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.5))
dev.off()

#B lineage
ibs <- read.table("B_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

pdf("IBS_B_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.5))
dev.off()

#D lineage
ibs <- read.table("D_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

pdf("IBS_A_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.7))
dev.off()

setwd("/Users/guoyafei/Documents/Lulab/Project-2-Migration/基本统计/IBS/")
data <- read.table("126Landrace_barley.txt",header=T,stringsAsFactors = F)
library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)
library(psych)
visit.x<-data$Logititude
visit.y<-data$Latitude
mp <- NULL
mapworld <- borders("world",colour = "gray50",fill="white") 
mp <- ggplot()+mapworld+ylim(-60,90)

mp2 <- mp+geom_point(aes(x=data$Logititude,y=data$Latitude,color=data$wild_A),size=2)+scale_size(range=c(1,1))+ scale_colour_gradient(low = "lightblue",high = "darkblue")+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
mp2 <- mp+geom_point(aes(x=data$Logititude,y=data$Latitude,color=data$dome_A),size=2)+scale_size(range=c(1,1))+ scale_colour_gradient(low = "lightblue",high = "darkblue")+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
mp2 <- mp+geom_point(aes(x=data$Logititude,y=data$Latitude,color=data$free_A),size=2)+scale_size(range=c(1,1))+ scale_colour_gradient(low = "lightblue",high = "darkblue")+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

mp2 <- mp+geom_point(aes(x=data$Logititude,y=data$Latitude,color=data$wild_B),size=2)+scale_size(range=c(1,1))+ scale_colour_gradient(low = "lightblue",high = "darkblue")+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
mp2 <- mp+geom_point(aes(x=data$Logititude,y=data$Latitude,color=data$dome_B),size=2)+scale_size(range=c(1,1))+ scale_colour_gradient(low = "lightblue",high = "darkblue")+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
mp2 <- mp+geom_point(aes(x=data$Logititude,y=data$Latitude,color=data$free_B),size=2)+scale_size(range=c(1,1))+ scale_colour_gradient(low = "lightblue",high = "darkblue")+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

mp2 <- mp+geom_point(aes(x=data$Logititude,y=data$Latitude,color=data$Strangulata),size=2)+scale_size(range=c(1,1))+ scale_colour_gradient(low = "lightblue",high = "darkblue")+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

mp3 <- mp2+ guides(fill=guide_legend(title=NULL))
mp3

library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)
library(psych)

setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/02_ABD_IBS/")
#ABD
data <- read.table("ABD_ibs.txt",sep="\t",header=T,stringsAsFactors = F)

AB <- data[which(data$value=="1" | data$value=="2" | data$value=="5" ),]
AB$AB_mean <- (AB$A_mean_ibs_ABD+AB$B_mean_ibs_ABD)/2
AB$type <- as.factor(AB$type)
D <- data[which(data$value=="2"),]
#A & B lineage
free <- AB[which(AB$value=="1"),]
ggplot(AB, aes(Logititude, Latitude,group=value,shape=type,fill=type))+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(color=AB_mean),size=2.5)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))

ggplot(AB, aes(Logititude, Latitude))+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(color=B_mean_ibs_ABD),size=1.5)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))
#D lineage
ggplot(D, aes(Logititude, Latitude))+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(color=D_mean_ibs_ABD),size=1.5)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))

#####wild_emmer, domesticated emmer, freethreshing Tetra, strangulata and 125 Landrace IBS distribution
#服务器工作目录：yafei@203:/data2/yafei/003_project3/Project3/CS_Vmap
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/02_ABD_IBS")
#读取group名的文件
path <- "/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/14_Subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)
#28个
nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_Durum","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_Khorasan_wheat","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_Polish_wheat","AABB_Rivet_wheat","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
#25个:7,11,16,17是freethreshing
#nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_freethreshing","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_freethreshing","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_freethreshing","AABB_freethreshing","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
A <- read.table("Alineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
rownames(A) <- A$Dxy
B <- read.table("Blineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
rownames(B) <- B$Dxy
D <- read.table("Dlineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
rownames(D) <- D$Dxy
#freethreshing Tetra Alineage
Alineage <- A[which(A$Dxy %in% c(data[[7]][,1],data[[11]][,1],data[[16]][,1],data[[17]][,1])),which(colnames(A)%in%data[[12]][,1])]
#freethreshing Tetra Blineage
Blineage <- B[which(B$Dxy %in% c(data[[7]][,1],data[[11]][,1],data[[16]][,1],data[[17]][,1])),which(colnames(B)%in%data[[12]][,1])]
#strangulata Dlineage
Dlineage <- D[which(D$Dxy %in% c(data[[20]][,1])),which(colnames(D)%in%data[[12]][,1])]
Alineage$mean <- apply(Alineage,1,mean)
Blineage$mean <- apply(Blineage,1,mean)
Dlineage$mean <- apply(Dlineage,1,mean)
newD <- Dlineage[-33,]
#wild emmer Alineage
W_A <- A[which(A$Dxy %in% c(data[[26]][,1])),which(colnames(A)%in%data[[12]][,1])]
#wild emmer Blineage
W_B <- B[which(B$Dxy %in% c(data[[26]][,1])),which(colnames(B)%in%data[[12]][,1])]
#domesticated emmer Alineage
D_A <- A[which(A$Dxy %in% c(data[[6]][,1])),which(colnames(A)%in%data[[12]][,1])]
#domesticated emmer Blineage
D_B <- B[which(B$Dxy %in% c(data[[6]][,1])),which(colnames(B)%in%data[[12]][,1])]
W_A$mean <- apply(W_A,1,mean)
W_B$mean <- apply(W_B,1,mean)
D_A$mean <- apply(D_A,1,mean)
D_B$mean <- apply(D_B,1,mean)
cbind(W_A$mean,W_B$mean)
cbind(D_A$mean,D_B$mean)

#plot
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/02_ABD_IBS/")
#ABD
data <- read.table("ABD_ibs.txt",sep="\t",header=T,stringsAsFactors = F)
AB <- data[which(data$value=="1" | data$value=="2" | data$value=="5" ),]
AB$AB_mean <- (AB$A_mean_ibs_ABD+AB$B_mean_ibs_ABD)/2
AB$type <- as.factor(AB$type)
D <- data[which(data$value=="2"),]
#1color
ggplot(AB, aes(Logititude, Latitude,group=value,shape=type,fill=type))+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(color=AB_mean),size=2.5)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))

#3color
free <- AB[which(AB$value=="1"),]
wild <- AB[which(AB$value=="2"),]
dome <- AB[which(AB$value=="5"),]
dome <- dome[-16,]
library(ggpubr)
a <- ggplot()+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(x=free$Logititude, y=free$Latitude,group=free$value,color=free$AB_mean),size=2.5,shape=15)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))
b <- ggplot()+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(x=wild$Logititude, y=wild$Latitude,group=wild$value,color=wild$AB_mean),size=2.5,shape=16)+
  scale_colour_gradient(low = "orange",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))
c <- ggplot()+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(x=dome$Logititude, y=dome$Latitude,group=dome$value,color=dome$AB_mean),size=2.5,shape=17)+
  scale_colour_gradient(low = "blue",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))
ggarrange(a,b,c, 
          labels = c("Freethreshing", "Wild", "Domesticated"),
          ncol = 1, nrow = 3)



#画图文件：yafei@203
#/data2/yafei/Project3/CS_Vmap/A_ibs.txt
#/data2/yafei/Project3/CS_Vmap/B_ibs.txt
#/data2/yafei/Project3/CS_Vmap/D_ibs.txt

#工作目录：yafei@203:/data2/yafei/Project3/CS_Vmap

#热图  
library("corrplot")
#A lineage
ibs <- read.table("A_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

#A:0.5 B:0.5 D:0.7
pdf("IBS_A_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.5))
dev.off()

#B lineage
ibs <- read.table("B_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

pdf("IBS_B_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.5))
dev.off()

#D lineage
ibs <- read.table("D_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

pdf("IBS_A_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.7))
dev.off()


setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/02_ABD_IBS/")
#Xinjiang_wheat
data <- read.table("B_xinjiang.txt",header=T,stringsAsFactors = F)
#Macha
data <- read.table("Macha.txt",header=T,stringsAsFactors = F)
#Persian
data <- read.table("Persian.txt",header=T,stringsAsFactors = F)
#ABD
data <- read.table("ABD_ibs.txt",sep="\t",header=T,stringsAsFactors = F)

library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)
library(psych)

#Rivet_wheat
R <- data[which(data$Type=="Rivet_wheat" | data$Type=="Persian" ),]
#Polish_wheat
P <- data[which(data$Type=="Polish_wheat" | data$Type=="Xinjiang"),]
#Khorasan_wheat
K <- data[which(data$Type=="Khorasan_wheat"),]
#Durum
D <- data[which(data$Type=="Wild_emmer" | data$Type=="Macha"),]
D <- data[which(data$Type=="Domesticated_emmer" ),]

#Landrace
L <- data[which(data$Value=="5"),]

P$Value <- as.factor(P$Value)
R$Value <- as.factor(R$Value)
data$Value <- as.factor(data$Value)
D$Value <- as.factor(D$Value)

A <- data[which(data$value=="1"),]
D <- data[which(data$value=="2"),]
ggplot(D, aes(Logititude, Latitude))+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(color=D_mean_ibs_spelt),size=1.5)+
  #geom_point(aes(color=mean,shape=Value),size=3)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))


