setwd("/Users/guoyafei/Documents/02_VmapIII/07_Baypass/")
library(reshape)
data <- read.table("climate_summary.txt",header=F,stringsAsFactors = F)
colnames(data) <- c("ID","Variable","Value")
data2 <- cast(data,ID~Variable)
write.table(data2,"35_variable.txt",sep="\t",quote=F,row.names = F)
#整体做PC
pc <- as.matrix(scale(data2[, c(2:36)]))
rm1<- cor(pc)
rs1<- eigen(rm1)
#碎石图, 决定取PC的个数(依据Eigenvalue(Principal Component Variance)), 该例中取4个PC
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=1,
     cex.lab=1,
     cex.axis=1,
     lty=1,
     lwd=1,
     xlab = "PC",
     ylab="Eigenvalue")
#提取特征向量（loading）得到PC值，loadings的输出结果为载荷，是主成分对应于原始变量的系数，即Q矩阵。
U <- as.matrix(rs1$vectors)
PC <- pc %*% U
princ_PC <- PC[,1:10]
rownames(princ_PC) <- data2[,1]
#提取样本PC值的另一种方法
#pca_solar <- prcomp(pc, center = TRUE,scale. = TRUE)
#pca_solar$x
#添加样本经纬度和海拔高度
lonlat <- read.table("lonlat.txt",header=T,stringsAsFactors = F)
rownames(lonlat) <- lonlat$ID
all <- cbind(lonlat[,2:4], princ_PC[rownames(lonlat),])
#kmeans聚类
library(cluster)
library(factoextra)
#按列进行标准化
df = scale(all, center = T, scale = T)
#确定应该分几个cluster
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 20, linetype = 2)
wss <- (nrow(df)-1)*sum(apply(df,2,var))
for (i in 2:50) wss[i] <- sum(kmeans(df,centers=i)$withinss)
plot(1:50, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
abline(h=1000)

#kmeans聚类，标准化后
#set.seed(1234)
set.seed(12345)
km <- kmeans(df, 25, iter.max = 10000) #用于画地图
fviz_cluster(km, data = df,
             #palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid",
             star.plot = T, 
             repel = TRUE,
             ggtheme = theme_minimal())
map <- all[,1:3]
map$type <- as.numeric(km$cluster)
write.table(map,"baypass.type_10PC.txt",quote=F,sep="\t")
-------------------------------层次聚类---------
#求样本之间两两相似性，层次聚类
data2<- data2[,c(1:22)]
result <- dist(data2, method = "euclidean")
result2 <- dist(df, method = "euclidean")
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
data2$type <- cutree(result_hc, k=11)

----------------------------------------计算类群环境变量-----------------------------------
lat_mean <- tapply(data2[,21],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,22],data2$type,mean,na.rm = TRUE)
data2$cluster1 <- NA
data2$cluster2 <- NA
for(i in 1:219) {
  for(j in 1:8){
    if(data2[i,23] == j ){
      data2[i,24] <- as.numeric(lat_mean[j])
      data2[i,25] <- as.numeric(lon_mean[j])
    } 
  }
}
write.table(data2,"13_cluster.txt",sep="\t",row.names = T,quote=F)

----------------------------------------地图上展示聚类结果---------------------------------
#删减样本后，还剩471个样本。
data <- read.table("471_baypass_10PC_edited.txt",header=T,stringsAsFactors = F)
data <- read.table("/Users/guoyafei/Desktop/baypass/471_baypass_10PC_edited.txt", header=T, stringsAsFactors = F)
library(maps)
library(ggplot2)
mp<-NULL
mapworld<-borders("world",colour = "#E7EBF1",fill="#E7EBF1") 

p <- list()
for ( i in 1:25) {
  sub <- as.data.frame(data[which(data$type == i),])
  p[[i]]  <- ggplot(sub, aes(x=Longitude, y=Latitude))+
    mapworld+
    ylim(-50,80)+
    geom_point(data=sub,alpha = 0.7,color="#E69F00")+
    scale_size(range=c(1,1))+ 
    facet_grid(. ~ type)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme_bw()
}

EA <- c(2,6,7,13,20)
CA <- c(5,8,11,12,25)
WA <- c(1,3,16,17,18)
SA <- c(4,10,19)
EU <- c(9,15,22,23)
Other <- c(14,21,24)

test <- c(1,16,17)
data$type <- as.factor(data$type)
sub <- data[which(data$type %in% Other),]
ggplot(sub, aes(x=Longitude, y=Latitude, fill=type,color = type))+
  mapworld+
  ylim(-50,80)+
  geom_point(alpha = 0.7)+
  scale_size(range=c(1,1))+ 
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()

library(gridExtra)
pdf("25群体地理分布.pdf",width = 15,height = 10)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],nrow=5)
dev.off()


pop1:WA:1
pop2:EA:3
pop3:IA:2
pop4:SH:5
pop5:IA:2
pop6:EA:3
pop7:EA:3
pop8:IA:2
pop9:EU:4
pop10:SH:5
pop11:IA:2
pop12:IA:2
pop13:EA:3
pop14:AF:6
pop15:EU:4
pop16:WA:1
pop17:WA:1
pop18:WA:1
pop19:SH:5
pop20:EA:3
pop21:SAM:7
pop22:EU:4
pop23:EU:4
pop24:SAM:7
pop25:IA:2

#25个群体的数据分布
library(ggplot2)
setwd("/Users/guoyafei/Desktop/baypass")
data <- read.table("471_baypass_taxa.Info", header=T,stringsAsFactors = F)
M <- aggregate(data, by=list(data$type),FUN=mean)
a <- colnames(M)[c(3,4,5,8:42)]
pdf("相关性.pdf", width=7, height=3.5)
for(i in c(3,4,5,8:42)){
  #data <- read.table("471_baypass_taxa.Info", header=T,stringsAsFactors = F)
  data$type <- factor(data$type, levels=rownames(M[order(M[,i]),]),ordered = TRUE)
  p <- ggplot(data, aes(x=data$type, y = data[,i-1])) +
    geom_boxplot(color="blue",
                 fill="blue",
                 alpha=0.2) +
    ylab(colnames(data[,i-1,drop=F])) +
    xlab("Group")+
    theme_bw()
  print(p)
}
dev.off()

#25个群体的数据分布
library(ggplot2)
setwd("/Users/guoyafei/Desktop/baypass")
data <- read.table("471_baypass_taxa.Info", header=T,stringsAsFactors = F)
M <- aggregate(data, by=list(data$type),FUN=mean)
a <- colnames(M)[c(3,4,5,8:42)]
pdf("相关性.pdf", width=7, height=3.5)
for(i in c(3,4,5,8:42)){
  #data <- read.table("471_baypass_taxa.Info", header=T,stringsAsFactors = F)
  data$type <- factor(data$type, levels=rownames(M[order(M[,i]),]),ordered = TRUE)
  p <- ggplot(data, aes(x=data$type, y = data[,i-1])) +
    geom_boxplot(color="blue",
                 fill="blue",
                 alpha=0.2) +
    ylab(colnames(data[,i-1,drop=F])) +
    xlab("Group")+
    theme_bw()
  print(p)
}
dev.off()

################ PC-相关性 ################
library(corrgram)
setwd("/Users/guoyafei/Desktop/envgwas/baypass")
data <- read.table("baypass_uniq.txt", header = T, stringsAsFactors = F)
sub <- data[,-c(1,2,3,4,5,18)]
summary(lm(data$solar15_PC1~data$solar15_PC2))

corrgram(sub, order = F,lower.panel = NULL,upper.panel=panel.pie,text.panel = panel.txt, gap = 0.1,
         main="Correlogram of environment variables intercorrelations")







