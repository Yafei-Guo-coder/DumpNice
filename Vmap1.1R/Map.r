library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)
library(psych)

setwd("/Users/guoyafei/Desktop/")
label <- read.table("map2.txt", header=T, stringsAsFactors = F)

label <- label[!is.na(label$Latitude),]
label$Genome <- as.factor(label$Genome) 
label$Genome  = gsub("AA","#97FFFF",label$Genome)
label$Genome  = gsub("SS","#FFD700",label$Genome)
label$Genome  = gsub("DD","#FF6347",label$Genome)
label$Genome  = gsub("AABB","#8470FF",label$Genome)
label$Genome  = gsub("AABBDD","#D8BFD8",label$Genome)
label$color <- NA
label[which(label$Genome=="AA"),5] <- "#97FFFF"#蓝色
label[which(label$Genome=="SS"),5] <- "#FFD700"#黄色
label[which(label$Genome=="DD"),5] <- "#FF6347"#红色
label[which(label$Genome=="AABB"),5] <- "#8470FF"#紫色
label[which(label$Genome=="AABBDD"),5] <- "#D8BFD8"


col <- label$Genome
visit.x <- label$Logititude
visit.y <- label$Latitude

label$color <- NA
label[which(label$Region=="AF"),5] <- "#97FFFF"#蓝色
label[which(label$Region=="AM"),5] <- "#FFD700"#黄色
label[which(label$Region=="EA"),5] <- "#FF6347"#红色
label[which(label$Region=="EU"),5] <- "#8470FF"#紫色
label[which(label$Region=="SCA"),5] <- "#D8BFD8"
label[which(label$Region=="WA"),5] <- "pink"
visit.color <- label$color

mydata<-read.table("map.txt",header=T,sep="\t",fill=T)
visit.x<-mydata$longitude
visit.y<-mydata$latitude
mydata$color <- NA
mydata[which(mydata$type=="GBS"),4] <- "red"
mydata[which(mydata$type=="Vmap"),4] <- "blue"
visit.color <-mydata$color

mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70") 
mp<-ggplot()+mapworld+ylim(-90,90)

#利用ggplot呈现，同时地图纵坐标范围从-60到90
mp2<-mp+geom_point(aes(x=label$Logititude, y=label$Latitude, color=label$Genome ),size=1,)+scale_size(range=c(1,1))+ theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20)) +
  scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
                    breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
                    labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())
#绘制带点的地图，geom_point是在地图上绘制点，x轴为经度信息，y轴为纬度信息，size是将点的大小按照收集的个数确定，color为暗桔色，scale_size是将点变大一些
mp2<-mp+geom_point(aes(x=label$Logititude, y=label$Latitude,color=label$Genome),size=1)+scale_size(range=c(1,1))+ theme_classic()+
  scale_fill_manual(
    values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8","pink"), 
    breaks=c("AF", "AM", "EA", "EU","SCA","WA"),
    labels=c("AF", "AM", "EA", "EU","SCA","WA"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) +
  theme(panel.grid =element_blank()) + 
  theme(axis.ticks.y = element_blank())
#
mp3+theme(legend.position = )

mp3 

library(spatstat)
library(maps)
library(maptools)
library("FactoMineR")
library("ggplot2")
library("factoextra")
library(mapdata)
library(devtools)
library(REmap)
setwd("/Users/guoyafei/Documents/Lulab/Project-1-Tibet/china-province-border-data")

data <- read.table("7_clim.csv",header=T,stringsAsFactors=F,sep=",")
test <- data[, c(5,6,10,14,19,20)]
test.pr<-princomp(test,cor=TRUE)
summary(test.pr,loadings=TRUE)
load <- as.matrix(c(0.414,0.410,0.410,-0.383,-0.417,-0.415),7,1)
test <- as.matrix(test)
load <- as.matrix(load)
loading <- test%*%load
sub <- cbind(data[,2:3],loading)
sub$name <- c("左贡县","加查县","桑日县","朗县","隆子县","察隅县","察雅县")
geo <- sub[,c(1,2,4)]
data = data.frame(country = mapNames('西藏'),value = 50*sample(7)+200,stringsAsFactors = F)
data$value<-as.numeric(data$value)
pointData = data.frame(geo$name)
remapC(data,maptype = '西藏',color = 'skyblue',
       markPointData = pointData,
       markPointTheme = markPointControl(symbol = 'pin',
                                         symbolSize = 5,
                                         effect = F),
       geoData = geo)  #上海火锅店

#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
setwd("/Users/guoyafei/Documents/个人项目/Project-2-Migration/VRN/GBS/GBS_Vmap/")
library(pophelper)
library(ggplot2)
require(gridExtra)
sfiles <- list.files(path="/Users/guoyafei/Documents/个人项目/Project-2-Migration/VRN/GBS/GBS_Vmap/admixture/", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))

p1 <- plotQ(slist[1],returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="all",showindlab=T,showyaxis=T,showticks=T)
#-----------------------------------------------------------------------------------------------
require(reshape)
require (rworldmap)
require(rworldxtra)
data <- read.table("anno.txt",header=T,stringsAsFactors = F)
#聚类
#根据经纬度给样本聚类
library(cluster)
library(factoextra)
dataA <- data[,c(4,5)]
data2 <- dataA[!is.na(dataA$latitude),]
df = scale(data2)

#先求样本之间两两相似性
result <- dist(df, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")

data2$type <- cutree(result_hc, k=50)
lat_mean <- tapply(data2[,1],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,2],data2$type,mean,na.rm = TRUE)
data2$cluster1 <- NA
data2$cluster2 <- NA
for(i in 1:697) {
  for(j in 1:50){
    if(data2[i,3] == j ){
      data2[i,4] <- as.numeric(lat_mean[j])
      data2[i,5] <- as.numeric(lon_mean[j])
    } 
  }
}
write.table(data2,"data2.txt",sep="\t",row.names = F,quote=F)

data <- read.table("anno.txt",header=T,stringsAsFactors = F)
dataA <- data[,c(14,1,12,13,15,16)]
colnames(dataA)[1:6] <- c("type","id", "Latitude", "Longitude", "Sub.population.6", "value")
wheat2 = dataA[!is.na(dataA$Latitude),]

dataB <- data[,7:12]
colnames(dataB)[1:6] <- c("type","Latitude", "Longitude", "id", "Sub.population.6", "value")
wheat2 = dataB[!is.na(dataB$Latitude),]

dataD <- data[,13:18]
colnames(dataD)[1:6] <- c("type","Latitude", "Longitude", "id", "Sub.population.6", "value")
wheat2 = dataD[!is.na(dataD$Latitude),]

dataD <- data[,19:24]
colnames(dataD)[1:6] <- c("type","Latitude", "Longitude", "id", "Sub.population.6", "value")
wheat2 = dataD[!is.na(dataD$Latitude),]

dataD <- data[,25:30]
colnames(dataD)[1:6] <- c("type","Latitude", "Longitude", "id", "Sub.population.6", "value")
wheat2 = dataD[!is.na(dataD$Latitude),]

dataD <- data[,31:36]
colnames(dataD)[1:6] <- c("type","Latitude", "Longitude", "id", "Sub.population.6", "value")
wheat2 = dataD[!is.na(dataD$Latitude),] 

wheat = wheat2[!is.na(wheat2$Sub.population.6),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.6) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
#5
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4","5"),symbolSize=1,
        zColours=c('red','blue',"green","orange","purple"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#4
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4"),symbolSize=1,
        zColours=c('red','blue',"green","orange"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
