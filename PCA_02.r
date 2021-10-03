setwd("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Test_chr001/03_Climate/")
library("FactoMineR")
library("ggplot2")
library("factoextra")
data <- read.csv("19bio.txt",fill=TRUE,na.strings = "",header = T, sep="\t", stringsAsFactors = F)
rownames(data) <- data$sample
anno <- read.table("anno.txt", header=T, stringsAsFactors = F)
anno2 <- anno[match(data$sample,anno$VcfId),]
data$ploidy <- anno2$Ploidy
data$Region <- anno2$ContinentAbbreviation

dt_all <- as.data.frame(data[!duplicated(data, fromLast=TRUE), 2:24])
dt <- as.matrix(scale(data[!duplicated(data, fromLast=TRUE), 4:22]))

rm1<-cor(dt)
rm1
rs1<-eigen(rm1)
rs1
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

#将iris数据集的第5列数据合并进来；
df<-data.frame(PC, dt_all$ploidy, dt_all$Region)
head(df)


library(ggplot2)
#提取主成分的方差贡献率，生成坐标轴标题；
xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
#绘制散点图并添加置信椭圆；
p1<-ggplot(data = df,aes(x=PC1,y=PC2,color=df$dt_all.ploidy))+
  stat_ellipse(aes(fill=df$dt_all.ploidy),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p2<-ggplot(data = df,aes(x=PC1,y=PC2,color=df$dt_all.Region))+
  stat_ellipse(aes(fill=df$dt_all.Region),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p1
p2

#------
dt_all$lat_bin <- NA
dt_all[which(dt_all$lat >= -45 & dt_all$lat < -35),24] <- "[-45,-35)"
dt_all[which(dt_all$lat >= -35 & dt_all$lat < -25),24] <- "[-35,-25)"
dt_all[which(dt_all$lat >= -25 & dt_all$lat < -15),24] <- "[-25,-15)"
dt_all[which(dt_all$lat >= -15 & dt_all$lat < -5),24] <- "[-15,-5)"
dt_all[which(dt_all$lat >= -5 & dt_all$lat < 5),24] <- "[-5,5)"
dt_all[which(dt_all$lat >= 5 & dt_all$lat < 15),24] <- "[5,15)"
dt_all[which(dt_all$lat >= 15 & dt_all$lat < 25),24] <- "[15,25)"
dt_all[which(dt_all$lat >= 25 & dt_all$lat < 35),24] <- "[25,35)"
dt_all[which(dt_all$lat >= 35 & dt_all$lat < 45),24] <- "[35,45)"
dt_all[which(dt_all$lat >= 45 & dt_all$lat < 55),24] <- "[45,55)"
dt_all[which(dt_all$lat >= 55 & dt_all$lat < 65),24] <- "[55,65)"
dt_all[which(dt_all$lat >= 65 & dt_all$lat < 75),24] <- "[65,75)"
df$lat_bin <- dt_all$lat_bin
p2<-ggplot(data = df,aes(x=PC1,y=PC2,color=df$lat_bin))+
  stat_ellipse(aes(fill=df$dt_all.Region),
               type ="norm", geom ="polygon",alpha=0.2,color=NA)+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  guides(fill=F)+
  theme_classic()
p1
