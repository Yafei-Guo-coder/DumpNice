library(cluster)
library(factoextra)
#聚类
#根据经纬度给样本聚类
#data <- read.table("land.txt", header=T, sep="\t", stringsAsFactors = F)
#dataA <- data[,c(3,4)]
#data2 <- dataA[!is.na(dataA$Latitude),]

data <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv/V1/225env.txt", header=F,stringsAsFactors = F)
colname <- c("elevation","temp1","temp2","temp3","temp4","temp5","temp6","temp7","temp8","temp9","temp10","temp11","prec1","prec2","prec3","prec4","prec5","prec6","prec7","prec8","Latitude","Logititude")
rownames(data) <- data[,1]
data2 <- data[!is.na(data$V6),-1]
colnames(data2) <- colname
df = scale(data2,center = T,scale = T)
colnames(df) <- colname
#按列进行标准化
#先求样本之间两两相似性
result <- dist(df, method = "euclidean")
#使用指定距离来计算数据矩阵行之间的距离
#euclidean：欧几里得距离
#maximum：最大距离
#manhattan：绝对距离
#canberra：堪培拉距离
#minkowski：闵可夫斯基距离
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
#Ward: 最小方差方法旨在寻找紧凑的球形簇的完整的联动方法找到相似集群。
#有两种不同的算法，"ward.D"（相当于只沃德选择"ward"不执行沃德（1963）聚类准则）& "ward.D2"实现了标准（Murtagh的及Legendre 2014）在集群更新之前对差异进行平方。注意agnes(*, method="ward")对应于hclust(*, "ward.D2").
#median和centroid在不导致单调的距离测量，或者等效产生的树状图可以具有所谓的倒置或颠倒。

data2$type <- cutree(result_hc, k=30)
lat_mean <- tapply(data2[,21],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,22],data2$type,mean,na.rm = TRUE)
data2$cluster1 <- NA
data2$cluster2 <- NA
for(i in 1:224) {
  for(j in 1:30){
    if(data2[i,23] == j ){
      data2[i,24] <- as.numeric(lat_mean[j])
      data2[i,25] <- as.numeric(lon_mean[j])
    } 
  }
}
write.table(data2,"30_cluster.txt",sep="\t",row.names = F,quote=F)

library(maps)
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70") 
mp<-ggplot()+mapworld+ylim(-90,90)
mp_40<-mp+geom_point(aes(x=re$Longitude_40, y=re$Latitude_40))+
  scale_size(range=c(1,1))+ 
  theme_classic()