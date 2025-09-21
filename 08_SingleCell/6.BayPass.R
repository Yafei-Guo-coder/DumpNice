# 样本环境和地理关系聚类以及环境信息的提取
setwd("/Users/yafeiguo/Documents/项目/拟南芥愈伤/taxa")
library(cluster)
library(factoextra)
library(raster)
library(rgdal)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
library(sf)
library(terra)
library(maps)
#library(viridis)

setwd("/Users/yafeiguo/Documents/项目/拟南芥愈伤/taxa")
taxa <- read.table("samples.txt", header=T,sep="\t", stringsAsFactors = F)
lon <- as.numeric(taxa$Long)
lat <- as.numeric(taxa$Lat)
samples <- data.frame(lon,lat)
temp.data <- samples


elev <- rast("/Volumes/My Passport/worldclim/wc2.1_2.5m_elev.tif")
temp01 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif")
temp02 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_2.tif")
temp03 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_3.tif")
temp04 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_4.tif")
temp05 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_5.tif")
temp06 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_6.tif")
temp07 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_7.tif")
temp08 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_8.tif")
temp09 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_9.tif")
temp10 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_10.tif")
temp11 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_11.tif")
prec1 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_12.tif")
prec2 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_13.tif")
prec3 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_14.tif")
prec4 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_15.tif")
prec5 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_16.tif")
prec6 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_17.tif")
prec7 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_18.tif")
prec8 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_19.tif")
srad01 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_01.tif")
srad02 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_02.tif")
srad03 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_03.tif")
srad04 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_04.tif")
srad05 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_05.tif")
srad06 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_06.tif")
srad07 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_07.tif")
srad08 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_08.tif")
srad09 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_09.tif")
srad10 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_10.tif")
srad11 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_11.tif")
srad12 <- raster("/Volumes/My Passport/worldclim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_12.tif")

temp.data$temp01 <- extract(temp01, samples)
temp.data$temp02 <- extract(temp02, samples)
temp.data$temp03 <- extract(temp03, samples)
temp.data$temp04 <- extract(temp04, samples)
temp.data$temp05 <- extract(temp05, samples)
temp.data$temp06 <- extract(temp06, samples)
temp.data$temp07 <- extract(temp07, samples)
temp.data$temp08 <- extract(temp08, samples)
temp.data$temp09 <- extract(temp09, samples)
temp.data$temp10 <- extract(temp10, samples)
temp.data$temp11 <- extract(temp11, samples) 
temp.data$prec1 <- extract(prec1, samples)
temp.data$prec2 <- extract(prec2, samples)
temp.data$prec3 <- extract(prec3, samples)
temp.data$prec4 <- extract(prec4, samples)
temp.data$prec5 <- extract(prec5, samples)
temp.data$prec6 <- extract(prec6, samples)
temp.data$prec7 <- extract(prec7, samples)
temp.data$prec8 <- extract(prec8, samples)
temp.data$srad01 <- extract(srad01, samples)
temp.data$srad02 <- extract(srad02, samples)
temp.data$srad03 <- extract(srad03, samples)
temp.data$srad04 <- extract(srad04, samples)
temp.data$srad05 <- extract(srad05, samples)
temp.data$srad06 <- extract(srad06, samples)
temp.data$srad07 <- extract(srad07, samples)
temp.data$srad08 <- extract(srad08, samples)
temp.data$srad09 <- extract(srad09, samples)
temp.data$srad10 <- extract(srad10, samples)
temp.data$srad11 <- extract(srad11, samples) 
temp.data$srad12 <- extract(srad12, samples) 
temp.data$elev <- extract(elev, samples)
rownames(temp.data) <- taxa[,1]
temp.data$elev <- temp.data[,34][,2]

#write.table(temp.data, "sample_withenv.txt", sep="\t", row.names = T, quote=F)

#样本聚类
library(cluster)
library(factoextra)
#聚类
#根据经纬度给样本聚类
data <- read.table("sample_withenv.txt", header=T,sep="\t", stringsAsFactors = F)

#dataA <- data[,c(3,4)]
data <- data[!is.na(data$lat),]

df = scale(data,center = T,scale = T)
#按列进行标准化
#先求样本之间两两相似性
set.seed(123)
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
data$type <- cutree(result_hc, k=28)
#画图确定
library(gridExtra)
mp <- NULL
mapworld <- borders("world",colour = "gray70",fill="gray70") 
mp <- ggplot() + mapworld + ylim(-90,90)
#这个死活都不对，为什么产的后一个对象会把前一个对象覆盖掉，
#手动画也不对，这里的问题是使用了同样的sub,sub一变，之前对应计算的画图的值也会变。在画图的时候要定义data数据集，否则最终结果会被覆盖
p <- list()
for(i in 1:28) {
  sub_data <- data[which(data$type == i), ]
  p[[i]] <- mp + 
    geom_point(data = sub_data, aes(x = lon, y = lat), size = 0.5) +
    scale_size(range = c(0.5, 0.5)) +
    theme_classic() +
    ggtitle(paste("Type", i))  # 添加标题以区分不同type的图
}


pdf("cluster.pdf",width = 15,height = 10)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]],p[[27]],p[[28]],nrow=5)
dev.off()

#确定cluster 15,16,24,25 去掉（因为样本数量和地理分布的原因）
#cluster 
lat_mean <- tapply(data[,2],data$type,mean,na.rm = TRUE)
lon_mean <- tapply(data[,1],data$type,mean,na.rm = TRUE)
data$cluster1 <- NA
data$cluster2 <- NA

for(i in 1:1131) {
  for(j in 1:28){
    if(data[i,35] == j ){
      data[i,36] <- as.numeric(lat_mean[j])
      data[i,37] <- as.numeric(lon_mean[j])
    } 
  }
}

data$ID <- rownames(data)
out <- data[,c(38,35,1,2,36,37,3:34)]
write.table(out,"28_cluster.txt",sep="\t",row.names = F,quote=F)



#画地图
setwd("/Users/yafeiguo/Documents/项目/拟南芥愈伤/taxa")
data <- read.table("28_cluster.txt", header=T,stringsAsFactors = F)


library(gridExtra)
mp <- NULL
mapworld <- borders("world",colour = "gray80",fill="gray80") 
mp <- ggplot() + mapworld + ylim(-90,90)

p <- mp + 
    geom_point(data = data, aes(x = lon, y = lat), color = "orange", size = 0.5) +
    scale_size(range = c(0.5, 0.5)) +
    theme_classic() +
    ggtitle("test")  # 添加标题以区分不同type的图


pdf("cluster.pdf",width = 15,height = 10)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]],p[[27]],p[[28]],nrow=5)
dev.off()


#画群体等位基因频率和环境变量的相关性
setwd("/Users/yafeiguo/Documents/项目/拟南芥愈伤/")
library(ggplot2)
data <- read.table("特定位点的关联分析.txt", header = T, stringsAsFactors = F,check.names = F)

p1 <- ggplot(data, aes(data$`temp5-7`, data$`4-11557476`)) +
  geom_point( size=2) +
  theme_classic()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  labs(x = 'Temperature Annual Range', y = 'Population allele frequency (4-11557476)')+
  ggtitle("R-squared: 0.38; P value: 0.00076")

p2 <- ggplot(data, aes(data$`temp10-12`, data$`4-11557476`)) +
  geom_point( size=2) +
  theme_classic()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  labs(x = ' Mean Temperature of Warmest Quarter', y = 'Population allele frequency (4-11557476)')+
  ggtitle("R-squared: 0.3; P value: 0.0033")

p3 <- ggplot(data, aes(data$`prec2-15`, data$`4-11555993`)) +
  geom_point( size=2) +
  theme_classic()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  labs(x = 'Precipitation of Wettest Month', y = 'Population allele frequency (4-11555993)')+
  ggtitle("R-squared: 0.28; P value: 0.0047")

p4 <- ggplot(data, aes(data$`prec5-18`, data$`4-11555993`)) +
  geom_point( size=2) +
  theme_classic()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  labs(x = 'Precipitation of Wettest Quarter', y = 'Population allele frequency (4-11555993)')+
  ggtitle("R-squared: 0.31; P value: 0.0026")

p5 <- ggplot(data, aes(data$`prec8-21`, data$`4-11555993`)) +
  geom_point( size=2) +
  theme_classic()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  labs(x = 'Precipitation of Coldest Quarter', y = 'Population allele frequency (4-11555993)')+
  ggtitle("R-squared: 0.44; P value: 0.00023")











