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


