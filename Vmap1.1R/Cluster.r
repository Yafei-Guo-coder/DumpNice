#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)
require(gridExtra)
setwd("/Users/guoyafei/Documents/个人项目/Project-2-Migration/migration/add_ZNdata/Structure_map/Select/")
sfiles <- list.files(path="/Users/guoyafei/Documents/个人项目/Project-2-Migration/migration/add_ZNdata/Structure_map/Select/Select_Qmatrix", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))

p1 <- plotQ(slist[3:7],returnplot=T,exportplot=F,quiet=T,basesize=11,
            showindlab=T,showyaxis=T,showticks=T)

p1 <- plotQ(slist[3:7],returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="all",showindlab=F,showyaxis=T,showticks=T,sharedindlab=T,clustercol=brewer.pal(6, "Set2"))

p1 <- plotQ(slist[3:10],returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="Cluster1",showindlab=T,showyaxis=T,showticks=T,sharedindlab=T)

grid.arrange(p1$plot[[1]],p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],nrow=5)

grid.arrange(m,p1$plot[[1]],nrow=2)



#聚类
#根据经纬度给样本聚类
data <- read.table("land.txt", header=T, sep="\t", stringsAsFactors = F)
library(cluster)
library(factoextra)
dataA <- data[,c(3,4)]
data2 <- dataA[!is.na(dataA$Latitude),]
df = scale(data2)

#先求样本之间两两相似性
result <- dist(df, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")

data2$type <- cutree(result_hc, k=40)
lat_mean <- tapply(data2[,1],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,2],data2$type,mean,na.rm = TRUE)
data2$cluster1 <- NA
data2$cluster2 <- NA
for(i in 1:126) {
  for(j in 1:40){
    if(data2[i,3] == j ){
      data2[i,4] <- as.numeric(lat_mean[j])
      data2[i,5] <- as.numeric(lon_mean[j])
    } 
  }
}
write.table(data2,"BB_data2.txt",sep="\t",row.names = F,quote=F)

#map
require("RColorBrewer")
setwd("/Users/guoyafei/Documents/个人项目/Project-2-Migration/migration/add_ZNdata/Structure_map/Select/")
data <- read.table("select_land.txt",header=T,sep="\t",stringsAsFactors = F)
dataA <- data[,c(1,7,8,24,25)]
colnames(dataA)[1:5] <- c("type","Latitude", "Longitude", "Sub.population.6", "value")
wheat2 = dataA[!is.na(dataA$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.6),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.6) 
wheat_reshape2 <- as.data.frame(wheat_reshape)

#8
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4","5","6","7","8"),symbolSize=1,
        zColours=brewer.pal(8, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#6
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4","5","6"),symbolSize=1,
        zColours=brewer.pal(6, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#5
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4","5"),symbolSize=1,
        zColours=brewer.pal(5, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#4
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,36),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4"),symbolSize=1,
        zColours=brewer.pal(4, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#3
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=brewer.pal(3, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#2
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2'),symbolSize=1,
        zColours=brewer.pal(2, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")

par(mfcol=c(5,1))
