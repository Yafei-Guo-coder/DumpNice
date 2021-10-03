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
