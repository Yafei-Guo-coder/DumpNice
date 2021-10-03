data <- read.table("Down+Vmap.txt",header=T,stringsAsFactors=F)
library(ggmap)
library(ggplot2)
library(sp)
library(maptools)
library(maps)
register_google(key = "AIzaSyAFnTwa7kEJAqIeHU-Kw6qPZiPUU-zXCLI")
visit.x<-as.numeric(caseAll$Logititude)
visit.y<-as.numeric(caseAll$Latitude)
mp<-NULL #定义一个空的地图
mapworld<-borders("world",colour = "gray50",fill="white") #绘制基本地图
mp<-ggplot()+mapworld+ylim(-60,90)
#利用ggplot呈现，同时地图纵坐标范围从-60到90
mp2<-mp+geom_point(aes(x=visit.x,y=visit.y, size=caseAll$count),color="darkorange")+scale_size(range=c(1,1))
#绘制带点的地图，geom_point是在地图上绘制点，x轴为经度信息，y轴为纬度信息，size是将点的大小按照收集的个数确定，color为暗桔色，scale_size是将点变大一些
mp3<-mp2+labs(x="Longitude", y="Latitude")+
  theme(legend.position = "none",panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 25),axis.text.y = element_text(size = 25),panel.border = element_blank()) #将图例去掉
pdf("castAll.pdf",width=20,height=8)
mp3
dev.off()


###画世界地图###
data <- read.table("/Users/guoyafei/Desktop/dataset_accessions_WGS.txt",header=T,stringsAsFactors = F)
library(maps)
data$color  <- NA
#map1
data[which(data$Region == "AF"), 11 ] <- "darkblue"
data[which(data$Region == "AM"), 11 ] <- "brown"
data[which(data$Region == "EA"), 11 ] <- "purple"
data[which(data$Region == "EU"), 11 ] <- "yellow"
data[which(data$Region == "SCA"), 11 ] <- "red"
data[which(data$Region == "WA"), 11 ] <- "green"
pdf(file = "map1.pdf",width=14,height=7)
map(database="world",fill = T,col = "grey",interior = F,border=NA)
points(x = data$Logititude, y = data$Latitude, col =data$color, pch=20,cex=0.5)
legend("right",legend = c("AF","AM","EA","EU","SCA","WA"),col = c("darkblue","brown","purple","yellow","red","green",rgb(0,255/255,255/255),rgb(255/255,0,255/255)),pch=20,cex=1,box.lty=0)
dev.off()
#map2
data[which(data$Genome == "SS"), 12 ] <- "green"
data[which(data$Genome == "AA"), 12 ] <- "brown"
data[which(data$Genome == "DD"), 12 ] <- "purple"
data[which(data$Genome == "AABB"), 12 ] <- "yellow"
data[which(data$Genome == "AABBDD"), 12 ] <- "red"
pdf(file = "map2.pdf",width=14,height=7)
map(database="world",fill = T,col = "grey",interior = F,border=NA)
points(x = data$Logititude, y = data$Latitude, col =data$color, pch=20,cex=0.5)
legend("right",legend = c("SS","AA","DD","AABB","AABBDD"),col = c("green","brown","purple","yellow","red",rgb(0,255/255,255/255),rgb(255/255,0,255/255)),pch=20,cex=1,box.lty=0)
dev.off()

