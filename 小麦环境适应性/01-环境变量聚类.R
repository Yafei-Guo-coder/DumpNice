library(corrgram)
data <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)
sub <- data[,c(3,6:40)]
sub2 <- sub[,c(1:9,11:14,22:36)]
sub3 <-  sub2[!duplicated(sub2, fromLast=TRUE), ] 

##################################################环境变量的PCA分析######################################################
#方法1（无效，由于存在相关性方向，无法很好的区分）
dt = t(scale(sub3, center = T, scale = T))
rm1 <- cor(dt)
rs1<- eigen(rm1)
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
     cex=1,
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
colnames(PC) <- paste("PC",1:269,sep="")
PC <- as.data.frame(PC)
library(ggplot2)

#提取主成分的方差贡献率，生成坐标轴标题；
pc1<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
pc2<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
pc3<-paste0("PC3(",round(Proportion_of_Variance[3]*100,2),"%)")

ggplot(data = PC,aes(x=PC1,y=PC2))+
  geom_point()+
  labs(x=pc1,y=pc2,color="")+
  geom_text(aes(label = rownames(PC)), size = 5)+
  guides(fill=F)+
  theme_classic()
ggplot(data = PC,aes(x=PC2,y=PC3))+
  geom_point()+
  labs(x=pc2,y=pc3,color="")+
  geom_text(aes(label = rownames(PC)), size = 5)+
  guides(fill=F)+
  theme_classic()
#另一种画图方式
plot(PC$PC1,PC$PC2,col='white')
text(PC$PC1,PC$PC2, rownames(PC), cex=1.5)
title("Eigenvector plot of baseball data")
arrows(0, 0, PC$PC1,PC$PC2, cex=0.5, col="red", length=0.1)

#方法2(有效，不使用相关性方向去掉)
df = scale(sub3, center = T, scale = T)
rm1 <- cor(df)
rm2 <- abs(rm1)
rs1<- eigen(rm2)
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
     cex=1,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
(U<-as.matrix(rs1$vectors))
#进行矩阵乘法，获得PC score；
PC <- rm1 %*% U
PC <- rm2 %*% U
colnames(PC) <- paste("PC",1:29,sep="")
PC <- as.data.frame(PC)
library(ggplot2)

#提取主成分的方差贡献率，生成坐标轴标题；
pc1<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
pc2<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
pc3<-paste0("PC3(",round(Proportion_of_Variance[3]*100,2),"%)")

ggplot(data = PC,aes(x=PC1,y=PC2))+
  geom_point()+
  labs(x=pc1,y=pc2,color="")+
  geom_text(aes(label = rownames(PC)), size = 5)+
  guides(fill=F)+
  theme_classic()

ggplot(data = PC,aes(x=PC2,y=PC3))+
  geom_point()+
  labs(x=pc2,y=pc3,color="")+
  geom_text(aes(label = rownames(PC)), size = 5)+
  guides(fill=F)+
  theme_classic()

#另一种画图方式
plot(PC$PC1,PC$PC2,col='white')
text(PC$PC1,PC$PC2, rownames(PC), cex=1)
title("Eigenvector plot of baseball data")
arrows(0, 0, PC$PC1,PC$PC2, cex=0.5, col="red", length=0.1)

###############################################环境变量的层次聚类分析####################################################
#方法1（无效，基于欧几里得距离，把相关性很强但方向相反的变量距离酸的最远）
library(factoextra)
df = t(scale(sub3, center = T, scale = T))
result <- dist(df,method = "euclidean")
#这个方法无法根据相关性计算距离，另一个amap包里的Dist方法可以，但是没有安装成功。
#按列进行标准化
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
plot(result_hc, hang=-1, cex=.8, main="Average Linkage Clustering")

#方法2（有效，基于相关性算距离，相关性越高，不管方向，距离越近）
df = scale(sub3, center = T, scale = T)
cor_matrix <- cor(df)
cluster_result <- hclust(as.dist(1 - abs(cor_matrix)), method = "ward.D2")
plot(cluster_result, hang=-1, cex=.8, main="Average Linkage Clustering")

################################################环境变量的相关性分析####################################################
library(corrgram)
#方法1（对应于层次聚类的方法一，无效）
vars2 <- c("soil8","prec7","temp8","prec3","prec6","soil13",
           "prec1","prec2","prec5","soil10","prec8","soil12","soil9",
           "temp3","temp9","temp1","temp11","temp6","temp10","temp5","prec4","solar2",
           "soil11","temp4","temp7","Elevation","temp2","solar1","solar3")

#方法2（对应于层次聚类的方法二，有效）
vars2 <- c("temp9","temp1","temp11","temp6","Elevation","temp10","temp5",
           "soil13","prec1","prec2","prec5","prec8","soil11","soil12","soil9","soil10","temp3","temp4","temp7",
           "soil8","solar2","prec7","temp8","temp2","solar1","solar3","prec4","prec3","prec6")
#corrplot
corrgram(sub3[,vars2], order = F,lower.panel = panel.shade,upper.panel=panel.pie,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")

#################################################环境变量的热图分析######################################################
#无法说明问题
library(pheatmap)
pheatmap(scale(sub3),cluster_rows = F,clustering_distance_cols  = "correlation", border_color = "white",cutree_cols  = 3)

##############################################三种环境变量的主成分分析###################################################
#library(FactoMineR)
#library(factoextra)
library(psych)
library(ggplot2)

EA <- c(2,6,7,13,20)
CA <- c(5,8,11,12,25)
WA <- c(1,3,16,17,18)
SA <- c(4,10,19)
EU <- c(9,15,22,23)
AM <- c(21,24)
AF <- c(14)
color <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#999999","#000000")

type1_var <- c("solar1","solar2","solar3","prec3","prec4","prec6","prec7","temp2","temp8")
type2_var <- c("temp1","temp5","temp6","temp9","temp10","temp11","Elevation")
type3_var <- c("prec1","prec2","prec5","prec8","temp3","temp4","temp7","soil9","soil10","soil11","soil12","soil13")

#type1 <- sub2[,which(colnames(sub2) %in% type1_var)]
#type2 <- sub2[,which(colnames(sub2) %in% type2_var)]
#type3 <- sub2[,which(colnames(sub2) %in% type3_var)]
#PC_type1 <- principal(type1, nfactors=3,rotate="none")
#PC_type2 <- principal(type2, nfactors=3,rotate="none")
#PC_type3 <- principal(type3, nfactors=3,rotate="none")

dt = scale(sub2, center = T, scale = T)
type1 <- dt[,which(colnames(dt) %in% type1_var)]
type2 <- dt[,which(colnames(dt) %in% type2_var)]
type3 <- dt[,which(colnames(dt) %in% type3_var)]

rm1 <- cor(type3)
rs1<- eigen(rm1)
val <- rs1$values
Proportion_of_Variance <- val/sum(val)
pc1<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
pc2<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
pc3<-paste0("PC3(",round(Proportion_of_Variance[3]*100,2),"%)")

U<-as.matrix(rs1$vectors)
PC <- type3 %*% U
colnames(PC) <- paste("PC",1:dim(PC)[2],sep="")
PC <- as.data.frame(PC)
PC$type <- as.factor(data$type)
PC$region <- NA
PC[which(PC$type %in% EA), dim(PC)[2]] <- "EA"
PC[which(PC$type %in% EU), dim(PC)[2]] <- "EU"
PC[which(PC$type %in% CA), dim(PC)[2]] <- "CA"
PC[which(PC$type %in% SA), dim(PC)[2]] <- "SA"
PC[which(PC$type %in% WA), dim(PC)[2]] <- "WA"
PC[which(PC$type %in% AM), dim(PC)[2]] <- "AM"
PC[which(PC$type %in% AF), dim(PC)[2]] <- "AF"
M <- aggregate(PC, by=list(PC$type),FUN=mean)
PC$type <- factor(PC$type, levels=rownames(M[order(M$PC1),]),ordered = TRUE)

ggplot(PC, aes(x=type,y=PC1,color=region,fill=region)) +
  geom_boxplot(alpha=0.8) +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  xlab("Group")+
  theme_classic()

###############################################环境距离和遗传距离的相关性############################################

#方式一：(在计算环境距离/地理距离和遗传距离的相关性时使用），根据环境PC选变量的问题：环境PC最相关的变量，在基因组上相应的位点很少，并不能很好的反应环境适应性。除了prec是PC2，其余都是PC1.
#  temp11(0.98);temp1(0.96);temp6(0.96)(temp1鉴定的又多又好)
#  solar3(0.99);solar1(0.94)(solar1鉴定的又多又好)
#  prec1(0.98);prec5(0.85);prec6(0. 85);prec3(0.83)
#    更新成prec3 prec4 prec6上面prec1;prec5鉴定的位点都很少
#  soil9(0.94);soil12(0.94);soil13(0.93);soil11(-0.83);soil10(soil12鉴定的又多又好)
#    更新成soil12(0.94)上面soil9，soil13，soil11鉴定位点情况都不好,位点也比较少
#方式二：（计算群体结构，定义最终适应性位点的时候使用）根据基因组贡献程度以及环境之间的相关性(小于0.3)，选择环境变量(V4)
#  temp8;temp9;temp4,temp5,temp2,elevation,temp7,temp3
#  solar1;solar2
#  prec7;prec8;prec4,prec3,prec6（剔除了prec8）
#  soil8;soil7;soil10;soil9;soil6,soil12;soil3(因为生物学意义，去掉soil6;soil8;soil7)

#环境变量的筛选
library(ggplot2)
setwd("/Users/guoyafei/Desktop/3type/fst")
data <- read.table("471_baypass_taxa.Info", header=T,row.names = 1,stringsAsFactors = F)
env <- data[,c(4:3,6:40)]

library(tidyr)
library(reshape2)
library ("geosphere")

#地理位置距离
geo_sub <- data[,c("Longitude","Latitude","type")]
M <- aggregate(geo_sub, by=list(geo_sub$type),FUN=mean)
row.names(M) <- paste("pop",M$type,sep="")
geo_sub <- M[,colnames(M)[2:3],drop=F]
result3 <- as.data.frame(as.matrix(distm(geo_sub, fun = distGeo)))
row.names(result3) <- paste("pop",1:25,sep="")
colnames(result3) <- paste("pop",1:25,sep="")
result3$id <- rownames(result3)
c<- melt(result3,id="id")
c$seq <- paste(c$id,c$variable,sep="-")
geo_order <- c[order(c$seq),]

#方式一：用多环境合并计算环境的距离（多环境根据PC最相关的环境来确定，
#all
#all <- c(4,6,8,10,11,15,16,17,18,26,27,29,30,32,37)
#old
#all <- c(4,8,9,11,17,27,29,30,32,37)
#now
all <- c(1,4,6,7,10,13,15,16,17,18,26,34,35,8,9,11,12,27,28,29,33,3,30,31,32,36,37)

#old
#prec <- c(4,8,9,11)
#solar <- c(4,27,29)
#temp <- c(4,30,32,37)

#now
#prec <- c(4,6,7,10,13,15,16,17,18,26,34,35,38)
solar <- c(1,24)
temp <- c(1,28)
prec <- c(1,14)

seq <- list(all,prec,solar,temp)

name <- c("all","prec","solar","temp")
#file <- paste(name,c("PC1V3_fst.txt","PC1V2_fst.txt","PC1_fst.txt","PC1_fst.txt"),sep="")
#file <- c("allPC1V3_fst.txt","precPC1V2_fst.txt","type1-R-V5_fst.txt","type2-R-V5_fst.txt")

file <- c("allPC1V3_fst.txt","type3-R-V5_fst.txt","type1-R-V5_fst.txt","type2-R-V5_fst.txt")
file_out <- paste(name,".pdf",sep="")
for ( i in c(1:4)) {
  #环境距离
  env2 <- env[,seq[[i]]]
  M <- aggregate(env2, by=list(env2$type),FUN=mean)
  row.names(M) <- paste("pop",M$type,sep="")
  env_sub <- M[,colnames(env2)[-1],drop=F]
  #df3 = scale(geo_sub,center = T,scale = T)
  #df3 <- geo_sub 
  #colnames(df3) <- colnames(geo_sub)
  result3 <- as.data.frame(as.matrix(dist(env_sub, method = "euclidean")))
  row.names(result3) <- paste("pop",1:25,sep="")
  colnames(result3) <- paste("pop",1:25,sep="")
  result3$id <- rownames(result3)
  c<- melt(result3,id="id")
  c$seq <- paste(c$id,c$variable,sep="-")
  env_order <- c[order(c$seq),]
  
  #遗传距离
  #PC1
  data <- read.table(file[i], header=T,stringsAsFactors = F)
  data2 <- data[,c(1,3,2,4,5)]
  colnames(data2) <- colnames(data)
  #meanfst
  all <- rbind(data,data2)[,c(2,3,4)]
  duijiao <- read.table("对角线.txt",header=T,stringsAsFactors = F)
  all <- rbind(all,duijiao)
  all[which(all$value < 0),3] <- 0
  all$seq <- paste(all$id1,all$id2,sep="-")
  fst_order <- all[order(all$seq),]
  out <- cbind(geo_order,fst_order)
  colnames(out) <- c("id1","id2","geo","seq1","id3","id4","fst","seq2")
  out$type <- name[i]
  a <- cor.test(out$geo,out$fst)
  
  out2 <- cbind(env_order,fst_order)
  colnames(out2) <- c("id1","id2","env","seq1","id3","id4","fst","seq2")
  out2$type <- name[i]
  a2 <- cor.test(out2$env,out2$fst)
  
  #遗传距离
  #随机变异
  data <- read.table("shuf_fst.txt", header=T,stringsAsFactors = F)
  data2 <- data[,c(1,3,2,4,5)]
  colnames(data2) <- colnames(data)
  #meanfst
  all <- rbind(data,data2)[,c(2,3,4)]
  duijiao <- read.table("对角线.txt",header=T,stringsAsFactors = F)
  all <- rbind(all,duijiao)
  all[which(all$value < 0),3] <- 0
  all$seq <- paste(all$id1,all$id2,sep="-")
  shuf_order <- all[order(all$seq),]
  shuf_out <- cbind(geo_order,shuf_order)
  colnames(shuf_out) <- c("id1","id2","geo","seq1","id3","id4","fst","seq2")
  shuf_out$type <- "shuf"
  b <- cor.test(shuf_out$geo,shuf_out$fst)
  
  shuf_out2 <- cbind(env_order,shuf_order)
  colnames(shuf_out2) <- c("id1","id2","env","seq1","id3","id4","fst","seq2")
  shuf_out2$type <- "shuf"
  b2 <- cor.test(shuf_out2$env,shuf_out2$fst)
  
  #plot:geo_fst
  pdf(file_out[i],height=3,width=4)
  til <- paste(a$p.value,as.numeric(a$estimate),b$p.value,as.numeric(b$estimate),sep=",")
  p <- rbind(out,shuf_out)
  pp <- ggplot(p, aes(x=geo/100000, y=fst/(1-fst),color=type)) +
    geom_point(size=0.5,alpha=0.3)+
    theme_bw()+
    xlab("Geographic distance(100km)")+
    ggtitle(til)+
    geom_smooth(method = "lm", se=TRUE, 
                formula = y ~ x) +
    scale_color_manual(values = c("#E69F00","#56B4E9"))
  
  #plot:env_fst
  q <- rbind(out2,shuf_out2)
  til2 <- paste(a2$p.value,as.numeric(a2$estimate),b2$p.value,as.numeric(b2$estimate),sep=",")
  qq <- ggplot(q, aes(x=env, y=fst/(1-fst),color=type)) +
    geom_point(size=0.5,alpha=0.3)+
    theme_bw()+
    xlab("Environmental distance")+
    ggtitle(til2)+
    geom_smooth(method = "lm", se=TRUE, 
                formula = y ~ x) +
    scale_color_manual(values = c("#E69F00","#56B4E9"))
  print(pp)
  print(qq)
  dev.off()
}

#方式二：用环境PC1计算环境的距离(这种方法效果不好)
setwd("/Users/guoyafei/Desktop/GF/等位基因频率变化")
PC1 <- read.table("/Users/guoyafei/Desktop/baypass/群体变量.txt",header=T,stringsAsFactors = F)
PC <- PC1[,c(1,4,5,2,3)]

name <- c("all","prec","soil","solar","temp")
file <- paste(name,"-R-V4_fst.txt",sep="")
file_out <- paste(name,"-all-R-PC2.pdf",sep="")
#所有环境做PC,结果存储在*-all-R-PC*.pdf
PC_file <- paste(name,"-PC.txt",sep="")
#只有筛选过（有一定的遗传贡献以及R<0.3)的环境变量做PC）效果稍微好一点，结果存储在*-sub-R-PC*.pdf
#PC_file <- paste(name,"-sub-PC.txt",sep="")
for ( i in c(2:5)) {
  #环境距离
  #env2 <- PC[,i,drop=F]
  env <- read.table(PC_file[i],header=T,stringsAsFactors = F)
  env2 <- env[,2,drop=F]
  result3 <- as.data.frame(as.matrix(dist(env2, method = "euclidean")))
  row.names(result3) <- paste("pop",1:25,sep="")
  colnames(result3) <- paste("pop",1:25,sep="")
  result3$id <- rownames(result3)
  c<- melt(result3,id="id")
  c$seq <- paste(c$id,c$variable,sep="-")
  env_order <- c[order(c$seq),]
  
  #遗传距离
  #PC1
  data <- read.table(file[i], header=T,stringsAsFactors = F)
  data2 <- data[,c(1,3,2,4,5)]
  colnames(data2) <- colnames(data)
  #meanfst
  all <- rbind(data,data2)[,c(2,3,4)]
  duijiao <- read.table("对角线.txt",header=T,stringsAsFactors = F)
  all <- rbind(all,duijiao)
  all[which(all$value < 0),3] <- 0
  all$seq <- paste(all$id1,all$id2,sep="-")
  fst_order <- all[order(all$seq),]
  out <- cbind(geo_order,fst_order)
  colnames(out) <- c("id1","id2","geo","seq1","id3","id4","fst","seq2")
  out$type <- name[i]
  a <- cor.test(out$geo,out$fst)
  
  out2 <- cbind(env_order,fst_order)
  colnames(out2) <- c("id1","id2","env","seq1","id3","id4","fst","seq2")
  out2$type <- name[i]
  a2 <- cor.test(out2$env,out2$fst)
  
  #遗传距离
  #随机变异
  data <- read.table("shuf_fst.txt", header=T,stringsAsFactors = F)
  data2 <- data[,c(1,3,2,4,5)]
  colnames(data2) <- colnames(data)
  #meanfst
  all <- rbind(data,data2)[,c(2,3,4)]
  duijiao <- read.table("对角线.txt",header=T,stringsAsFactors = F)
  all <- rbind(all,duijiao)
  all[which(all$value < 0),3] <- 0
  all$seq <- paste(all$id1,all$id2,sep="-")
  shuf_order <- all[order(all$seq),]
  shuf_out <- cbind(geo_order,shuf_order)
  colnames(shuf_out) <- c("id1","id2","geo","seq1","id3","id4","fst","seq2")
  shuf_out$type <- "shuf"
  b <- cor.test(shuf_out$geo,shuf_out$fst)
  
  shuf_out2 <- cbind(env_order,shuf_order)
  colnames(shuf_out2) <- c("id1","id2","env","seq1","id3","id4","fst","seq2")
  shuf_out2$type <- "shuf"
  b2 <- cor.test(shuf_out2$env,shuf_out2$fst)
  
  #plot:geo_fst
  pdf(file_out[i],height=3,width=4)
  til <- paste(a$p.value,as.numeric(a$estimate),b$p.value,as.numeric(b$estimate),sep=",")
  p <- rbind(out,shuf_out)
  pp <- ggplot(p, aes(x=geo/100000, y=fst/(1-fst),color=type)) +
    geom_point(size=0.5,alpha=0.3)+
    theme_bw()+
    xlab("Geographic distance(100km)")+
    ggtitle(til)+
    geom_smooth(method = "lm", se=TRUE, 
                formula = y ~ x) +
    scale_color_manual(values = c("#E69F00","#56B4E9"))
  
  #plot:env_fst
  q <- rbind(out2,shuf_out2)
  til2 <- paste(a2$p.value,as.numeric(a2$estimate),b2$p.value,as.numeric(b2$estimate),sep=",")
  qq <- ggplot(q, aes(x=env, y=fst/(1-fst),color=type)) +
    geom_point(size=0.5,alpha=0.3)+
    theme_bw()+
    xlab("Environmental distance")+
    ggtitle(til2)+
    geom_smooth(method = "lm", se=TRUE, 
                formula = y ~ x) +
    scale_color_manual(values = c("#E69F00","#56B4E9"))
  print(pp)
  print(qq)
  dev.off()
}

###########################################曼哈顿图(lfmm & baypass)#################################################
library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Desktop/环境适应性位点/曼哈顿/")

input <- c("baypass.solar1_21chr.txt","lfmm.solar1_21chr.txt")
output <- c("baypass.solar1_21chr.pdf","lfmm.solar1_21chr.pdf")

for (i in 1){
  gwasResults2 <- read.table(input[i], header=T, stringsAsFactors = F)
  gwasResults <- gwasResults2[sample(nrow(gwasResults2), 200000), ]
  don <- gwasResults %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) 
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  p <- ggplot(don, aes(x=BPcum, y=P)) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    scale_color_manual(values = rep(c("grey","skyblue","#E69F00"), 7)) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +   
    #Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=15),
      axis.text.y=element_text(size=15),
      axis.title.y=element_text(size = 15),
      axis.title.x=element_text(size = 15),
    )
  #scale_y_continuous(limits = c(0,7))+
  #geom_point(data=point,aes(x=BPcum,y=-log10(P)),color="red")
  pdf(output[i],height = 2.5,width = 15)
  print(p)
  dev.off()
}

#########################################################环境变量的三角图###################################################
#首先对环境变量取PC1的值



















