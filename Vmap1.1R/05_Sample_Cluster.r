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