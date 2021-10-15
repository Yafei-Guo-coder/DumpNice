#all_prec_tavg----
#library("FactoMineR")
#library("factoextra")
data <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/PCA_prec_tavg.txt", header=T, stringsAsFactors = F)
#1-2:经纬度
#3-14:降水
#15-26:温度
#27-45:19bio
#46:elevation
#47,48,49:VcfID, Group, Ploidy
dt <- as.matrix(scale(data[!duplicated(data, fromLast=TRUE), 3:26]))
rm1 <- cor(dt, use = "complete.obs")
rs1 <- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:24,sep="")
head(PC)

#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

df<-data.frame(PC, data$Group, taxa$type, data$Ploidy, data$TreeValidatedGroupbySubspecies)
head(df)
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df, aes(x=PC1, y=PC2, color=df$data.Ploidy))+
  stat_ellipse(aes(fill=df$data.Ploidy),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()

p2<-ggplot(data = df,aes(x=PC1, y=PC2, color=df$data.TreeValidatedGroupbySubspecies))+
  stat_ellipse(aes(fill=df$data.TreeValidatedGroupbySubspecies),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
#聚类
#根据经纬度给样本聚类
library(cluster)
library(factoextra)
taxa1 <- data[,c(1:2)]
taxa <- taxa1[!duplicated(taxa1, fromLast=TRUE),]
df2 = scale(taxa)
#先求样本之间两两相似性
result <- dist(df2, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
taxa$type <- cutree(result_hc, k=4)
#lat_mean <- tapply(taxa[,1],taxa$type,mean,na.rm = TRUE)
#lon_mean <- tapply(taxa[,2],taxa$type,mean,na.rm = TRUE)
#taxa$cluster1 <- NA
#taxa$cluster2 <- NA
#for(i in 1:541) {
#  for(j in 1:20){
#    if(taxa[i,3] == j ){
#      taxa[i,4] <- as.numeric(lat_mean[j])
#      taxa[i,5] <- as.numeric(lon_mean[j])
#    } 
# }
#}
env1 <- data[,c(1:26)]
env2 <- env1[!duplicated(env1, fromLast=TRUE),]
#dt <- as.matrix(scale(data[!duplicated(data, fromLast=TRUE), 3:26]))
env <- as.matrix(env2[,3:26])
rm1 <- cor(env, use = "complete.obs")
rs1 <- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- env %*% U
colnames(PC) <- paste("PC",1:24,sep="")
head(PC)

#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

df<-data.frame(PC,  taxa$type)

p2<-ggplot(data = df,aes(x=PC1, y=PC2, color=as.factor(df$taxa.type)))+
  stat_ellipse(aes(fill=as.factor(df$taxa.type)),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p2


#AABBDD_prec_tavg----
data <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/PCA_prec_tavg.txt", header=T, stringsAsFactors = F)
#聚类
#根据经纬度给样本聚类
library(cluster)
library(factoextra)
ABD1 <- data[which(data$Ploidy=="AABBDD"),1:26]
ABD2 <- ABD1[!duplicated(ABD1, fromLast=TRUE),]
ABD_taxa <- ABD2[,1:2]
df2 = scale(ABD_taxa)
#先求样本之间两两相似性
result <- dist(df2, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
ABD_taxa$type <- cutree(result_hc, k=6)

#PCA分析
ABD <- ABD2[,3:26]
dt <- as.matrix(scale(ABD))
#452 pos
rm1 <- cor(dt, use = "complete.obs")
rs1 <- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:24,sep="")
head(PC)

#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

df<-data.frame(PC,ABD_taxa$type)
head(df)
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df, aes(x=PC1, y=PC2, color=as.factor(df$ABD_taxa.type)))+
  stat_ellipse(aes(fill=as.factor(df$ABD_taxa.type)),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p1
#生成地理分类（根据气温和降水）文件
ABD_taxa$k_3 <- cutree(result_hc, k=3)
ABD_taxa$k_4 <- cutree(result_hc, k=4)
ABD_taxa$k_5 <- cutree(result_hc, k=5)
ABD_taxa$k_6 <- cutree(result_hc, k=6)
ABD_taxa$PC1 <- df$PC1
ABD_taxa$PC2 <- df$PC2
write.table(ABD_taxa,"/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/452_ABD_pos_clust.txt", row.names = F,quote = F)


#AABBDD_prec----
data <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/PCA_prec_tavg.txt", header=T, stringsAsFactors = F)
#聚类
#根据经纬度给样本聚类
library(cluster)
library(factoextra)
ABD1 <- data[which(data$Ploidy=="AABBDD"),1:14]
ABD2 <- ABD1[!duplicated(ABD1, fromLast=TRUE),]
ABD_taxa <- ABD2[,1:2]
df2 = scale(ABD_taxa)
#先求样本之间两两相似性
result <- dist(df2, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
ABD_taxa$type <- cutree(result_hc, k=3)

#PCA分析
ABD <- ABD2[,3:14]
dt <- as.matrix(scale(ABD))
#452 pos
rm1 <- cor(dt, use = "complete.obs")
rs1 <- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:12,sep="")
head(PC)

#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

df<-data.frame(PC,ABD_taxa$type)
head(df)
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df, aes(x=PC1, y=PC2, color=as.factor(df$ABD_taxa.type)))+
  stat_ellipse(aes(fill=as.factor(df$ABD_taxa.type)),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p1
#生成地理分类（根据降水）文件
ABD_taxa$k_3 <- cutree(result_hc, k=3)
ABD_taxa$k_4 <- cutree(result_hc, k=4)
ABD_taxa$k_5 <- cutree(result_hc, k=5)
ABD_taxa$k_6 <- cutree(result_hc, k=6)
ABD_taxa$PC1 <- df$PC1
ABD_taxa$PC2 <- df$PC2
write.table(ABD_taxa,"/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/452_ABD_prec_clust.txt", row.names = F,quote = F)


#AABBDD_tavg
#AABBDD_19bio----
data <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/PCA_prec_tavg_19bio.txt", header=T, stringsAsFactors = F)
#1-2:经纬度
#3-14:降水
#15-26:温度
#27-45:19bio
#46:elevation
#47,48,49:VcfID, Group, Ploidy
#聚类
#根据经纬度给样本聚类
library(cluster)
library(factoextra)
ABD1 <- data[which(data$Ploidy=="AABBDD"),1:45]
ABD2 <- ABD1[!duplicated(ABD1, fromLast=TRUE),]
ABD_taxa <- ABD2[,1:2]
df2 = scale(ABD_taxa)
#先求样本之间两两相似性
result <- dist(df2, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
ABD_taxa$type <- cutree(result_hc, k=6)

#PCA分析
ABD <- ABD2[,27:45]
dt <- as.matrix(scale(ABD))
#452 pos
rm1 <- cor(dt, use = "complete.obs")
rs1 <- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:19,sep="")
head(PC)

#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

df<-data.frame(PC,ABD_taxa$type)
head(df)
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df, aes(x=PC1, y=PC2, color=as.factor(df$ABD_taxa.type)))+
  stat_ellipse(aes(fill=as.factor(df$ABD_taxa.type)),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p1
#生成地理分类（根据19bio）文件
ABD_taxa$k_3 <- cutree(result_hc, k=3)
ABD_taxa$k_4 <- cutree(result_hc, k=4)
ABD_taxa$k_5 <- cutree(result_hc, k=5)
ABD_taxa$k_6 <- cutree(result_hc, k=6)
ABD_taxa$PC1 <- df$PC1
ABD_taxa$PC2 <- df$PC2
write.table(ABD_taxa,"/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/452_ABD_19bio_clust.txt", row.names = F,quote = F)


#AABBDD_19bio_纬度海拔聚类----
data <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/PCA_prec_tavg_19bio_eleva.txt", header=T, stringsAsFactors = F)
#1-2:经纬度
#3-14:降水
#15-26:温度
#27-45:19bio
#46:elevation
#47,48,49:VcfID, Group, Ploidy
#聚类
#根据纬度和海拔给样本聚类
library(cluster)
library(factoextra)
ABD1 <- data[which(data$Ploidy=="AABBDD"),1:46]
ABD2 <- ABD1[!duplicated(ABD1, fromLast=TRUE),]
ABD_taxa <- ABD2[,c(2,46)]
df2 = scale(ABD_taxa)
#先求样本之间两两相似性
result <- dist(df2, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
ABD_taxa$type <- cutree(result_hc, k=6)

#PCA分析
ABD <- ABD2[,27:45]
dt <- as.matrix(scale(ABD))
#452 pos
rm1 <- cor(dt, use = "complete.obs")
rs1 <- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:19,sep="")
head(PC)

#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

df<-data.frame(PC,ABD_taxa$type)
head(df)
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df, aes(x=PC1, y=PC2, color=as.factor(df$ABD_taxa.type)))+
  stat_ellipse(aes(fill=as.factor(df$ABD_taxa.type)),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p1
#生成地理分类（根据19bio）文件
ABD_taxa$k_3 <- cutree(result_hc, k=3)
ABD_taxa$k_4 <- cutree(result_hc, k=4)
ABD_taxa$k_5 <- cutree(result_hc, k=5)
ABD_taxa$k_6 <- cutree(result_hc, k=6)
ABD_taxa$PC1 <- df$PC1
ABD_taxa$PC2 <- df$PC2
ABD_taxa$lon <- ABD2$lon

write.table(ABD_taxa,"/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/Out/TXT/452_ABD_19bio_纬度海拔clust.txt", row.names = F,quote = F)

#AABBDD_19bio_纬度海拔以及19bio变量聚类----
data <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/PCA_prec_tavg_19bio_eleva.txt", header=T, stringsAsFactors = F)
#1-2:经纬度
#3-14:降水
#15-26:温度
#27-45:19bio
#46:elevation
#47,48,49:VcfID, Group, Ploidy
#聚类
#根据纬度和海拔给样本聚类
library(cluster)
library(factoextra)
ABD1 <- data[which(data$Ploidy=="AABBDD"),1:46]
ABD2 <- ABD1[!duplicated(ABD1, fromLast=TRUE),]
ABD_taxa <- ABD2[,c(2,27:46)]
df2 = scale(ABD_taxa)
#先求样本之间两两相似性
result <- dist(df2, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
ABD_taxa$type <- cutree(result_hc, k=8)

#PCA分析
ABD <- ABD2[,27:45]
dt <- as.matrix(scale(ABD))
#452 pos
rm1 <- cor(dt, use = "complete.obs")
rs1 <- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:19,sep="")
head(PC)

#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

df<-data.frame(PC,ABD_taxa$type)
head(df)
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df, aes(x=PC1, y=PC2, color=as.factor(df$ABD_taxa.type)))+
  stat_ellipse(aes(fill=as.factor(df$ABD_taxa.type)),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p1
#生成地理分类（根据19bio）文件
ABD_taxa$k_3 <- cutree(result_hc, k=3)
ABD_taxa$k_4 <- cutree(result_hc, k=4)
ABD_taxa$k_5 <- cutree(result_hc, k=5)
ABD_taxa$k_6 <- cutree(result_hc, k=6)
ABD_taxa$k_7 <- cutree(result_hc, k=7)
ABD_taxa$k_8 <- cutree(result_hc, k=8)
ABD_taxa$k_9 <- cutree(result_hc, k=9)
ABD_taxa$k_10 <- cutree(result_hc, k=10)
ABD_taxa$PC1 <- df$PC1
ABD_taxa$PC2 <- df$PC2
ABD_taxa$lon <- ABD2$lon

write.table(ABD_taxa,"/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/Out/TXT/452_ABD_19bio_纬度海拔19bio_clust.txt", row.names = F,quote = F)

#AABBDD_19bio_经纬度聚类----
data <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/PCA_prec_tavg_19bio_eleva.txt", header=T, stringsAsFactors = F)
#1-2:经纬度
#3-14:降水
#15-26:温度
#27-45:19bio
#46:elevation
#47,48,49:VcfID, Group, Ploidy
#聚类
#根据纬度和海拔给样本聚类
library(cluster)
library(factoextra)
ABD1 <- data[which(data$Ploidy=="AABBDD"),1:46]
ABD2 <- ABD1[!duplicated(ABD1, fromLast=TRUE),]
ABD_taxa <- ABD2[,c(1,2,46)]
df2 = scale(ABD_taxa)
#先求样本之间两两相似性
result <- dist(df2, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
ABD_taxa$type <- cutree(result_hc, k=7)

#PCA分析
ABD <- ABD2[,27:45]
dt <- as.matrix(scale(ABD))
#452 pos
rm1 <- cor(dt, use = "complete.obs")
rs1 <- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
(Standard_deviation <- sqrt(val))
#计算方差贡献率和累积贡献率；
(Proportion_of_Variance <- val/sum(val))
(Cumulative_Proportion <- cumsum(Proportion_of_Variance))
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:19,sep="")
head(PC)

#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

df<-data.frame(PC,ABD_taxa$type)
head(df)
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df, aes(x=PC1, y=PC2, color=as.factor(df$ABD_taxa.type)))+
  stat_ellipse(aes(fill=as.factor(df$ABD_taxa.type)),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p1
#生成地理分类（根据19bio）文件
ABD_taxa$k_3 <- cutree(result_hc, k=3)
ABD_taxa$k_4 <- cutree(result_hc, k=4)
ABD_taxa$k_5 <- cutree(result_hc, k=5)
ABD_taxa$k_6 <- cutree(result_hc, k=6)
ABD_taxa$PC1 <- df$PC1
ABD_taxa$PC2 <- df$PC2
ABD_taxa$lon <- ABD2$lon

write.table(ABD_taxa,"/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/06_Environment/Out/TXT/452_ABD_19bio_经纬度海拔clust.txt", row.names = F,quote = F)
