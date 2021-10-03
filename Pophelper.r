library(pophelper)
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/06_Structure/Q_matrix2")
library(pophelper)
library(ggplot2)
require(gridExtra)
sfiles <- list.files(path="/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/06_Structure/Q_matrix2", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))

p1 <- plotQ(alignK(sortQ(slist)[1:3]),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11)
grid.arrange(p1$plot[[1]])

p1 <- plotQ(slist[c(5,6)],returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="all",showindlab=F,showyaxis=T,showticks=T)
p2 <- plotQ(slist[c(1,3)],imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="all",sharedindlab=F,showindlab=F,showyaxis=T,showticks=T)
grid.arrange(p1$plot[[1]],p1$plot[[2]],nrow=2)


tbl=read.table("B_free_land.prune.in.7.Q")
barplot(t(as.matrix(tbl)), col=rainbow(10),xlab="Individual #", ylab="Ancestry", border=NA)

require(reshape)
require (rworldmap)
require(rworldxtra)
data <- read.table("structure_ancestor.txt",header=T,stringsAsFactors = F)
#聚类
#根据经纬度给样本聚类
library(cluster)
library(factoextra)
dataA <- data[,c(32,33)]
data2 <- dataA[!is.na(dataA$Latitude.5),]
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
for(i in 1:141) {
  for(j in 1:50){
    if(data2[i,3] == j ){
      data2[i,4] <- as.numeric(lat_mean[j])
      data2[i,5] <- as.numeric(lon_mean[j])
    } 
  }
}



dataA <- data[,1:6]
colnames(dataA)[1:6] <- c("type","Latitude", "Longitude", "id", "Sub.population.6", "value")
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
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"4","5","6"),symbolSize=1,
        zColours=c('red','blue',"green","orange","purple"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#4
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4"),symbolSize=1,
        zColours=c('red','blue',"green","orange"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")

