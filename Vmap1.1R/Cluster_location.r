library(cluster)
library(factoextra)
data <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/215_location.txt",header=T,stringsAsFactors = F, fill=TRUE)
colnames(data) <- c("id","Latitude", "Longitude")
rownames(data) <- data$id

data2 <- data[,c(2,3)]
df = scale(data2)
#先求样本之间两两相似性
result <- dist(df, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
data2$type <- cutree(result_hc, k=40)
lat_mean <- tapply(data2[,1],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,2],data2$type,mean,na.rm = TRUE)
data2$Latitude_40 <- NA
data2$Longitude_40 <- NA
for(i in 1:215) {
  for(j in 1:40){
    if(data2[i,3] == j ){
      data2[i,4] <- as.numeric(lat_mean[j])
      data2[i,5] <- as.numeric(lon_mean[j])
    } 
  }
}
re <- data2[,c(3,4,5)]
write.table(re,"cluster_40.txt",quote=F,row.names = T)

library(maps)
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70") 
mp<-ggplot()+mapworld+ylim(-90,90)
mp_40<-mp+geom_point(aes(x=re$Longitude_40, y=re$Latitude_40))+
  scale_size(range=c(1,1))+ 
  theme_classic()
