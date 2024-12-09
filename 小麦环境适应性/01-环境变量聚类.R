library(corrgram)
setwd("/Users/guoyafei/Desktop/3type/location")
data <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info2", header=T, row.names = 1, stringsAsFactors = F)
sub <- data[,c(3,6:13,27:40)]
#sub2 <- sub[,c(1:9,11:14,22:36)]
sub3 <-  sub[!duplicated(sub, fromLast=TRUE), ] 

####################################### 环境变量的PCA分析 ######################################################
#方法1（无效，由于存在相关性方向，无法很好的区分）-------
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

#方法2(有效，不使用相关性方向，去掉)-----------
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
colnames(PC) <- paste("PC",1:28,sep="")
PC <- as.data.frame(PC)
library(ggplot2)
library(scatterplot3d)
my_color = c("#66C2A5FF", "#FC8D62FF", "#8DA0CBFF")
colors = my_color[as.numeric(PC$type)]

type1 <- c("solar1","solar2","solar3","prec3","prec4","prec6","prec7","temp2","temp8")
type2 <- c("temp1","temp5","temp6","temp9","temp10","temp11","Elevation")
type3 <- c("prec1","prec2","prec5","prec8","temp3","temp4","temp7","soil9","soil10","soil11","soil12","soil13")

PC$type <- "NA"

PC[which(PC$name %in% type1), 30] <- "1"
PC[which(PC$name %in% type2), 30] <- "2"
PC[which(PC$name %in% type3), 30] <- "3"

p1 = scatterplot3d(PC[,1:3], color = colors, main="Environmental Variables", pch = 16)
legend(p1$xyz.convert(8.5, 2.5, 5), legend = levels(iris$Species),
       col = my_color, pch = 16)
text(p1$xyz.convert(PC[,1:3]),labels = PC[,29],
     cex= 0.7, col = "black",adj=c(1,-1),font=2)


#提取主成分的方差贡献率，生成坐标轴标题；
pc1 <- paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
pc2 <- paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
pc3 <- paste0("PC3(",round(Proportion_of_Variance[3]*100,2),"%)")

ggplot(data = PC,aes(x=PC1,y=PC2,color=type))+
  geom_point()+
  labs(x=pc1,y=pc2,color="")+
  geom_text(aes(label = rownames(PC)), size = 5)+
  guides(fill=F)+
  theme_classic()

ggplot(data = PC,aes(x=PC2,y=PC3,color=type))+
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

####################################### 环境变量的层次聚类分析 ####################################################
#方法1（无效，基于欧几里得距离，把相关性很强但方向相反的变量距离算的最远）
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
cluster_result <- hclust(as.dist(1 - abs(cor_matrix)), method = "")

plot(cluster_result, hang=-1, cex=.8, main="Average Linkage Clustering")

####################################### 环境变量的相关性分析 ####################################################
library(corrgram)
#28个环境变量
vars2 <- c("solar1","prec7","solar2","temp9","temp8","temp5","prec8","temp4","prec4","temp10",
           "solar3","temp1","temp7","prec3","temp6","temp2","prec6","temp11","soil11","soil13",
           "soil10","soil12","prec5","soil9","prec2","prec1")

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

####################################### 环境变量的热图分析 ######################################################
#无法说明问题
library(pheatmap)
pheatmap(scale(sub3),cluster_rows = F,clustering_distance_cols  = "correlation", border_color = "white",cutree_cols  = 3)
                        
####################################### 三种环境变量的主成分分析 ###################################################
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

rm1 <- cor(type1)
rs1<- eigen(rm1)
val <- rs1$values
Proportion_of_Variance <- val/sum(val)
pc1<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
pc2<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
pc3<-paste0("PC3(",round(Proportion_of_Variance[3]*100,2),"%)")

U<-as.matrix(rs1$vectors)
PC <- type1 %*% U
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
type <- M[,c(1:5)]
write.table(type,"solar-PC.txt",quote=F,row.names = F,sep="\t")

ggplot(PC, aes(x=type,y=PC1,color=region,fill=region)) +
  geom_boxplot(alpha=0.8) +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  xlab("Group")+
  theme_classic()

####################################### 环境距离和遗传距离的相关性 ############################################
#方式一：(在计算环境距离/地理距离和遗传距离的相关性时使用），根据环境PC选变量的问题：环境PC最相关的变量，在基因组上相应的位点很少，并不能很好的反应环境适应性。除了prec是PC2，其余都是PC1.
#  temp11(0.98);temp1(0.96);temp6(0.96)(temp1鉴定的又多又好)
#  solar3(0.99);solar1(0.94)(solar1鉴定的又多又好)
#  prec1(0.98);prec5(0.85);prec6(0. 85);prec3(0.83)
#    更新成prec3 prec4 prec6上面prec1;prec5鉴定的位点都很少
#  soil9(0.94);soil12(0.94);soil13(0.93);soil11(-0.83);soil10(soil12鉴定的又多又好)
#    更新成soil12(0.94)上面soil9，soil13，soil11鉴定位点情况都不好,位点也比较少
#环境变量的筛选
library(ggplot2)
library(tidyr)
library(reshape2)
library (geosphere)
setwd("/Users/guoyafei/Desktop/3type/fst")
data <- read.table("471_baypass_taxa.Info", header=T,row.names = 1,stringsAsFactors = F)
env <- data[,c(4:3,6:40)]

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
#取和环境变量PC1最相关的变量
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

####################################### 曼哈顿图(lfmm & baypass)#################################################
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

####################################### 群体对环境变量距离PC1的三元图###################################################
#首先对环境变量取PC1的值
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggtern)
#library (geosphere)
setwd("/Users/guoyafei/Desktop/3type/三元图")
data <- read.table("/Users/guoyafei/Desktop/3type/fst/471_baypass_taxa.Info", header=T,row.names = 1,stringsAsFactors = F)
env <- data[,c(4:3,6:40)]
sub <- data[,c(3,6:40)]
type1_var <- c("solar1","solar2","solar3","prec3","prec4","prec6","prec7","temp2","temp8")
type2_var <- c("temp1","temp5","temp6","temp9","temp10","temp11","Elevation")
type3_var <- c("prec1","prec2","prec5","prec8","temp3","temp4","temp7","soil9","soil10","soil11","soil12","soil13")

#pop pair distance header
dt = scale(sub, center = T, scale = T)
var <- list(type1_var <- type1_var, type2_var <- type2_var, type3_var <- type3_var)
a <- rep(paste("pop",1:25,sep=""), 25)
b <- paste("pop",rep(1:25,each=25),sep="")
c<-paste(b,a,sep="-")
all <- as.data.frame(cbind(b,a,c))
all <- all[order(all$c),]

#环境变量PCA
for(i in c(1:3)){
  type <- dt[,which(colnames(dt) %in% var[[i]])]
  rm1 <- cor(type)
  rs1<- eigen(rm1)
  #提取结果中的特征值，即各主成分的方差；
  val <- rs1$values
  #计算方差贡献率和累积贡献率；
  Proportion_of_Variance <- val/sum(val)
  #提取主成分的方差贡献率，生成坐标轴标题；
  pc1<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
  pc2<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
  pc3<-paste0("PC3(",round(Proportion_of_Variance[3]*100,2),"%)")
  #提取结果中的特征向量(也称为Loadings,载荷矩阵)；
  U<-as.matrix(rs1$vectors)
  #进行矩阵乘法，获得PC score；
  PC <- type %*% U
  colnames(PC) <- paste("PC",1:dim(PC)[2],sep="")
  PC <- as.data.frame(PC)
  PC$type <- data$type
  M <- aggregate(PC, by=list(PC$type),FUN=mean)
  
  #pop pair distance PC1 distribution
  env_sub <- M[,c(2,3),drop=F]
  result3 <- as.data.frame(as.matrix(dist(env_sub, method = "euclidean")))
  row.names(result3) <- paste("pop",1:25,sep="")
  colnames(result3) <- paste("pop",1:25,sep="")
  result3$id <- rownames(result3)
  c<- melt(result3,id="id")
  c$seq <- paste(c$id,c$variable,sep="-")
  type_order <- c[order(c$seq),]
  all <- cbind(all,as.numeric(type_order[,3]))

}

#画图测试
ggplot(data = M,aes(x=PC1,y=PC2,color=as.factor(type)))+
  geom_point()+
  labs(x=pc1,y=pc2,color="")+
  guides(fill=F)+
  theme_classic()

#pop pair distance PC1 distribution
all_dist <- as.data.frame(all)
colnames(all_dist) <- c("pop1","pop2","pop","type1","type2","type3")
all_dist$pop1 <- as.character(all_dist$pop1)
all_dist$pop2 <- as.character(all_dist$pop2)
all_dist$pop <- as.character(all_dist$pop)
all_dist$type1 <- as.numeric(all_dist$type1)
all_dist$type2 <- as.numeric(all_dist$type2)
all_dist$type3 <- as.numeric(all_dist$type3)
sub_dist <- all_dist[which(all_dist$pop1 != all_dist$pop2 & all_dist$pop1 != "pop21" & all_dist$pop1 != "pop24"& all_dist$pop2 != "pop21" & all_dist$pop2 != "pop24"),]
rownames(sub_dist) <- sub_dist$pop
a <- sub_dist[,c(4:6)]
b <- a[!duplicated(a, fromLast=TRUE),]
sub_dist <- b
sub_dist$scale_type1 <- (sub_dist$type1-min(sub_dist$type1))/(max(sub_dist$type1)-min(sub_dist$type1))
sub_dist$scale_type2 <- (sub_dist$type2-min(sub_dist$type2))/(max(sub_dist$type2)-min(sub_dist$type2))
sub_dist$scale_type3 <- (sub_dist$type3-min(sub_dist$type3))/(max(sub_dist$type3)-min(sub_dist$type3))
test <- sub_dist[,c(4:6)]
test$all <- test$scale_type1 + test$scale_type2 + test$scale_type3
test$prop1 <- test$scale_type1/test$all
test$prop2 <- test$scale_type2/test$all
test$prop3 <- test$scale_type3/test$all
test$id <- c(1:dim(test)[1])
test$color <- "type7"
test[which(test$prop1 > 0.6),9] <- "type1"
test[which(test$prop2 > 0.6),9] <- "type2"
test[which(test$prop3 > 0.6),9] <- "type3"
test[which(test$prop1 < 0.2 & test$color == "type7"),9] <- "type4"
test[which(test$prop2 < 0.2 & test$color == "type7"),9] <- "type5"
test[which(test$prop3 < 0.2 & test$color == "type7"),9] <- "type6"

#检查群体对距离的分布
for(i in c(1:3)){
  out <- paste("type",i,"poppair.PC1.pdf",sep="")
  pdf(out)
  p <- ggplot(test, aes(test[,i]))+
    theme_bw()+
    geom_histogram(binwidth = 0.05,alpha = 0.6)
  print(p)
  dev.off()
}  

#简单测试
o <- test[,c(5,8)]
o$type <- "type1"
colnames(o) <- c("prop","id","type")
p <- test[,c(6,8)]
p$type <- "type2"
colnames(p) <- c("prop","id","type")
q <- test[,c(7,8)]
q$type <- "type3"
colnames(q) <- c("prop","id","type")
all <- rbind(o,p,q)
ggplot(all, aes( x = id, weight = prop, fill = type))+
  geom_bar(position = "dodge")

pdf("三元图.pdf")
ggtern(data=test,aes(scale_type1,scale_type2,scale_type3,color= color)) +  
  #stat_density_tern(aes(fill=..level.., alpha=..level..), geom='polygon') +
  #scale_fill_gradient2(high = "blue") + 
  scale_color_manual(values = c("#0072B2", "#D55E00","#CC79A7","#56B4E9","#009E73","#F0E442","#999999"))+
  geom_point() +
  theme_showarrows() +
  ggtitle("My Favorite Color") +
  xlab("Solar") + 
  ylab("Temperature") +
  zlab("Precipitation") +
  guides(color = "none", fill = "none", alpha = "none")
dev.off()

#分类图
color <- c("#0072B2", "#D55E00","#CC79A7","#56B4E9","#009E73","#F0E442","#999999")
var <- paste("type",1:7,sep="")
for(i in c(1:7)){
  type <- test[which(test$color == var[i]),c(5:8)]
  colnames(type)[1:3] <- c("Solar","Temperature","Precipitation")
  melt_type <- melt(type,id="id")
  out <- paste(var[i],"2","pdf",sep=".")
  pdf(out)
  p <- ggplot(melt_type, aes(x=as.factor(variable), y = value, group=as.factor(variable)))+
    geom_boxplot(fill=color[i])+
    theme_classic()+  
    theme(
      legend.position="none",
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=15),
      axis.text.y=element_text(size=15),
      axis.title.y=element_text(size = 18),
      axis.title.x=element_text(size = 18),
    )+
    ylim(0,1)+
    ylab("Relative Distance")+
    xlab("Environment Variables")
  print(p)
  dev.off()
}

#具体的群体变量对
a <- rownames(test[which(test$color == "type1"),])
b <- rownames(test[which(test$color == "type2"),])
c <- rownames(test[which(test$color == "type3"),])
d <- rownames(test[which(test$color == "type4"),])
e <- rownames(test[which(test$color == "type5"),])
f <- rownames(test[which(test$color == "type6"),])
g <- rownames(test[which(test$color == "type7"),])


####################################### 单个群体的环境变量PC1的三元图###################################################
#首先对环境变量取PC1的值
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggtern)
#library (geosphere)
setwd("/Users/guoyafei/Desktop/3type/三元图")
data <- read.table("/Users/guoyafei/Desktop/3type/fst/471_baypass_taxa.Info", header=T,row.names = 1,stringsAsFactors = F)
env <- data[,c(4:3,6:40)]
sub <- data[,c(3,6:40)]
type1_var <- c("solar1","solar2","solar3","prec3","prec4","prec6","prec7","temp2","temp8")
type2_var <- c("temp1","temp5","temp6","temp9","temp10","temp11","Elevation")
type3_var <- c("prec1","prec2","prec5","prec8","temp3","temp4","temp7","soil9","soil10","soil11","soil12","soil13")

#pop single header
dt = scale(sub, center = T, scale = T)
var <- list(type1_var <- type1_var, type2_var <- type2_var, type3_var <- type3_var)
all <- paste("pop",1:25,sep="")

#环境变量PCA
for(i in c(1:3)){
  type <- dt[,which(colnames(dt) %in% var[[i]])]
  #####获得PC score 方法一，与方法二结果一样，loading是指观测变量与主成分的相关系数
  #rm1 <- cor(type)
  #rs1<- eigen(rm1)
  #####提取结果中的特征值，即各主成分的方差；
  #val <- rs1$values
  #####计算方差贡献率和累积贡献率；
  #Proportion_of_Variance <- val/sum(val)
  #####提取主成分的方差贡献率，生成坐标轴标题；
  #pc1<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
  #pc2<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
  #pc3<-paste0("PC3(",round(Proportion_of_Variance[3]*100,2),"%)")
  #####提取结果中的特征向量；
  #U<-as.matrix(rs1$vectors)
  #进行矩阵乘法，获得主成分得分；
  #PC <- type %*% U
  
  #####获得PC score 方法二，与方法二结果一样
  PC_type <- principal(type, nfactors=5,rotate="none")
  #获得PC score；计算主成分得分
  PC <- PC_type$scores
  colnames(PC) <- paste("PC",1:dim(PC)[2],sep="")
  PC <- as.data.frame(PC)
  PC$type <- data$type
  M <- aggregate(PC, by=list(PC$type),FUN=mean)
  #pop single PC1 distribution
  env_sub <- M[,c(2),drop=F]
  row.names(env_sub) <- paste("pop",1:25,sep="")
  all <- cbind(all,as.numeric(env_sub[,1]))
}

#获取主成分与原始变量的线性关系
library(pheatmap)
out_var <- paste("type",1:7,sep="")
for(i in c(1:3)){
  type <- dt[,which(colnames(dt) %in% var[[i]])]
  PC_type <- principal(type, nfactors=4,rotate="none")
  #计算主成分得分系数
  corr <- PC_type$weights
  out <- paste(out_var[i],"corr","pdf",sep=".")
  pdf(out)
  p <- pheatmap(corr[,1,drop=F], show_rownames=TRUE,cluster_col = F,cluster_row = FALSE)
  print(p)
  dev.off()
}

#画图测试
ggplot(data = M,aes(x=PC1,y=PC2,color=as.factor(type)))+
  geom_point(size=10)+
  #labs(x=pc1,y=pc2,color="")+
  guides(fill=F)+
  theme_classic()

#pop single PC1 distribution
all_dist <- as.data.frame(all)
colnames(all_dist) <- c("pop","type1","type2","type3")
all_dist$pop <- as.character(all_dist$pop)
all_dist$type1 <- as.numeric(all_dist$type1)
all_dist$type2 <- as.numeric(all_dist$type2)
all_dist$type3 <- as.numeric(all_dist$type3)
sub_dist <- all_dist[which(all_dist$pop != "pop24" & all_dist$pop != "pop21"),]

rownames(sub_dist) <- sub_dist$pop
a <- sub_dist[,c(2:4)]
b <- a[!duplicated(a, fromLast=TRUE),]

sub_dist <- b

sub_dist$scale_type1 <- (sub_dist$type1-min(sub_dist$type1))/(max(sub_dist$type1)-min(sub_dist$type1))
sub_dist$scale_type2 <- (sub_dist$type2-min(sub_dist$type2))/(max(sub_dist$type2)-min(sub_dist$type2))
sub_dist$scale_type3 <- (sub_dist$type3-min(sub_dist$type3))/(max(sub_dist$type3)-min(sub_dist$type3))

test <- sub_dist[,c(4:6)]
test$all <- test$scale_type1 + test$scale_type2 + test$scale_type3
test$prop1 <- test$scale_type1/test$all
test$prop2 <- test$scale_type2/test$all
test$prop3 <- test$scale_type3/test$all
test$id <- c(1:dim(test)[1])

#检查单群体PC1的分布
for(i in c(1:3)){
  out <- paste("type",i,"single.PC1.pdf",sep="")
  pdf(out)
  p <- ggplot(test, aes(test[,i]))+
    theme_bw()+
    geom_histogram(binwidth = 0.05,alpha = 0.6)
  print(p)
  dev.off()
}



test$color <- "type7"
test[which(test$prop1 > 0.6),9] <- "type1"
test[which(test$prop2 > 0.6),9] <- "type2"
test[which(test$prop3 > 0.6),9] <- "type3"
test[which(test$prop1 < 0.2 & test$color == "type7"),9] <- "type4"
test[which(test$prop2 < 0.2 & test$color == "type7"),9] <- "type5"
test[which(test$prop3 < 0.2 & test$color == "type7"),9] <- "type6"

pdf("三元图2.single.pdf")
ggtern(data=test,aes(scale_type1,scale_type2,scale_type3,color= color)) +  
  #stat_density_tern(aes(fill=..level.., alpha=..level..), geom='polygon') +
  #scale_fill_gradient2(high = "blue") + 
  scale_color_manual(values = c("#0072B2", "#D55E00","#CC79A7","#56B4E9","#009E73","#F0E442","#999999"))+
  geom_point() +
  theme_showarrows() +
  ggtitle("My Favorite Color") +
  xlab("Solar") + 
  ylab("Temperature") +
  zlab("Precipitation") +
  guides(color = "none", fill = "none", alpha = "none")
dev.off()

#分类图
color <- c("#0072B2", "#D55E00","#CC79A7","#56B4E9","#009E73","#F0E442","#999999")
var <- paste("type",1:7,sep="")
for(i in c(1:7)){
  type <- test[which(test$color == var[i]),c(5:8)]
  colnames(type)[1:3] <- c("Solar","Temperature","Precipitation")
  melt_type <- melt(type,id="id")
  out <- paste(var[i],"single2","pdf",sep=".")
  pdf(out)
  p <- ggplot(melt_type, aes(x=as.factor(variable), y = value, group=as.factor(variable)))+
    geom_boxplot(fill=color[i])+
    theme_classic()+
    theme(
      legend.position="none",
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=15),
      axis.text.y=element_text(size=15),
      axis.title.y=element_text(size = 18),
      axis.title.x=element_text(size = 18),
    )+
    ylim(0,1)+
    ylab("Relative Distance")+
    xlab("Environment Variables")
  print(p)
  dev.off()
}

#具体的群体变量对
a <- rownames(test[which(test$color == "type1"),])
b <- rownames(test[which(test$color == "type2"),])
c <- rownames(test[which(test$color == "type3"),])
d <- rownames(test[which(test$color == "type4"),])
e <- rownames(test[which(test$color == "type5"),])
f <- rownames(test[which(test$color == "type6"),])
g <- rownames(test[which(test$color == "type7"),])

########################################################## 2023.11.20 ###################################################
#################### 重新计算环境变量的相关性，去掉相关性很高的变量以及跟其他变量相关性都不高的变量######################
library(corrgram)
setwd("/Users/guoyafei/Desktop/3type/location")
data <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)


sub <- data[,c(3,6:40)] #36个环境变量
sub2 <- sub[,c(1,4,6,8:9,12,14,22:27,29,34:36)] #17个环境变量

sub2 <- data[,c(3,6,7,8,10,11,12,16,17,18,26,27,29,30,31,32,33,36,37,40)]#20个环境变量
sub3 <-  sub2[!duplicated(sub2, fromLast=TRUE), ] 
sub3 <-  sub[!duplicated(sub, fromLast=TRUE), ] 

#去掉美洲和非洲的样本
sub <- data[which(data$type != 14 & data$type != 21 & data$type != 24),c(3,6:40)]
#PCA分析
#方法2(有效，不使用相关性方向，因为在与环境变量做关联时，无关环境是梯度上升还是梯度下降)
df = scale(sub3, center = T, scale = T)
rm1 <- cor(df)
rm2 <- abs(rm1)
rs1<- eigen(rm2)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
Standard_deviation <- sqrt(val)
#计算方差贡献率和累积贡献率；
Proportion_of_Variance <- val/sum(val)
Cumulative_Proportion <- cumsum(Proportion_of_Variance)
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
U <- as.matrix(rs1$vectors)
#进行矩阵乘法，获得PC score；
PC <- rm1 %*% U
PC <- rm2 %*% U
colnames(PC) <- paste("PC",1:20,sep="")
PC <- as.data.frame(PC)

library(ggplot2)
library(scatterplot3d)
my_color = c("#66C2A5FF", "#FC8D62FF", "#8DA0CBFF")

#type1 <- c("solar1","solar3","prec7","temp2","temp9")
#type2 <- c("temp1","temp10","solar2","temp8","Elevation")
#type3 <- c("prec5","prec8","prec3","temp7","soil9","soil11","soil13")
type1 <- c("temp1","temp5","temp6","temp9","temp10","temp11","Elevation")
type2 <- c("solar1","solar3","temp2","prec3","prec6","prec7")
type3 <- c("prec1","prec2","prec5","soil9","soil11","soil12","soil13")

PC$type <- "NA"

PC$name <- rownames(PC)
PC[which(PC$name %in% type1), 21] <- "1"
PC[which(PC$name %in% type2), 21] <- "2"
PC[which(PC$name %in% type3), 21] <- "3"

colors = my_color[as.numeric(PC$type)]
p1 = scatterplot3d(PC[,1:3], color = colors, main="Environmental Variables", pch = 16)
p1 = scatterplot3d(PC[,1:3], main="Environmental Variables", pch = 16)
legend(p1$xyz.convert(8.5, 2.5, 5), legend = levels(iris$Species),
       col = my_color, pch = 16)
text(p1$xyz.convert(PC[,1:3]),labels = PC[,22],
     cex= 0.7, col = "black",adj=c(1,-1),font=2)

#提取主成分的方差贡献率，生成坐标轴标题；
pc1 <- paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
pc2 <- paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")
pc3 <- paste0("PC3(",round(Proportion_of_Variance[3]*100,2),"%)")

ggplot(data = PC,aes(x=PC1,y=PC2,color=type))+
  geom_point()+
  labs(x=pc1,y=pc2,color="")+
  geom_text(aes(label = rownames(PC)), size = 5)+
  guides(fill=F)+
  theme_classic()+
  theme(legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))


ggplot(data = PC,aes(x=PC2,y=PC3,color=type))+
  geom_point()+
  labs(x=pc2,y=pc3,color="")+
  geom_text(aes(label = rownames(PC)), size = 5)+
  guides(fill=F)+
  theme_classic()+
  theme(legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))

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

library(FactoMineR)
library(factoextra)
library(psych)
PC2 <- principal(rm2, nfactors=5,rotate="none")
#PC_solar <- principal(solar, nfactors=3,rotate="none")
#PC_prec <- principal(prec, nfactors=3,rotate="none")
#PC_soil <- principal(soil, nfactors=3,rotate="none")

#另一种画图方式
plot(PC$PC1,PC$PC2,col='white')
text(PC$PC1,PC$PC2, rownames(PC), cex=1)
title("Eigenvector plot of baseball data")
arrows(0, 0, PC$PC1,PC$PC2, cex=0.5, col="red", length=0.1)

#############################################################环境变量的热图分析##############################################################
#无法说明问题
library(pheatmap)
pheatmap(scale(sub3),cluster_rows = F,clustering_distance_cols  = "correlation", border_color = "white",cutree_cols  = 3)

pheatmap(rm2, cluster_cols = T,cluster_rows = T, clustering_method = "ward",border_color = "white")
pheatmap(rm2, cluster_cols = T,cluster_rows = T, clustering_method = "ward.D2",border_color = "white")
pheatmap(rm2, cluster_cols = T,cluster_rows = T, clustering_method = "single",border_color = "white")
pheatmap(rm2, cluster_cols = T,cluster_rows = T, clustering_method = "complete",border_color = "white")
pheatmap(rm2, cluster_cols = T,cluster_rows = T, clustering_method = "average",border_color = "white")
pheatmap(rm2, cluster_cols = T,cluster_rows = T, clustering_method = "mcquitty",border_color = "white")
pheatmap(rm2, cluster_cols = T,cluster_rows = T, clustering_method = "median",border_color = "white")
pheatmap(rm2, cluster_cols = T,cluster_rows = T, clustering_method = "centroid",border_color = "white")

############################################################# 环境变量在群体的箱线图分布 ##############################################################

#环境变量的筛选
library(ggplot2)
setwd("/Users/guoyafei/Desktop/baypass")
data <- read.table("471_baypass_taxa.Info", header=T,row.names = 1,stringsAsFactors = F)
env <- data[,c(1:4,6:40)]

#PCA
library(FactoMineR)
library(factoextra)
library(psych)

#population regions
EA <- c(2,6,7,13,20)
CA <- c(5,8,11,12,25)
WA <- c(1,3,16,17,18)
SA <- c(4,10,19)
EU <- c(9,15,22,23)
AM <- c(21,24)
AF <- c(14)

env$region <- NA
env[which(env$type %in% EA),dim(env)[2]] <- "EA"
env[which(env$type %in% WA),dim(env)[2]] <- "WA"
env[which(env$type %in% SA),dim(env)[2]] <- "SA"
env[which(env$type %in% EU),dim(env)[2]] <- "EU"
env[which(env$type %in% CA),dim(env)[2]] <- "CA"
env[which(env$type %in% AM),dim(env)[2]] <- "AM"
env[which(env$type %in% AF),dim(env)[2]] <- "AF"

env$type <- as.factor(env$type)
for(i in c(1:3,5:39)){
  color <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#999999","#000000")
  #乌拉尔图: "#999999"（灰色）, 野生一粒:"#E69F00"(橘黄色), 栽培一粒:"#56B4E9"(明兰色), AABBDD:"#009E73"(橄榄绿), AABB:"#F0E442"（明黄）,"#0072B2"(天蓝), EA："#D55E00"(橘红色), SCA："#CC79A7"(皮粉色))
  outfile <- paste(colnames(env)[i],".pdf",sep="")
  pdf(outfile, height=3, width=7)
  p <- ggplot(env, aes(x=type,y=env[,i],color=region,fill=region)) +
    geom_boxplot(alpha=0.8) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    #scale_colour_discrete(breaks = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#999999","#000000"), labels = c('EU','WA','CA','SA','EA','AF','AM'))+
    ylab(colnames(env)[i]) +
    xlab("Group") +
    theme_bw()
  print(p)
  dev.off()
}

############################################################# 环境变量的含义 ##########################################################
library(raster)
library(rgdal)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
library(viridis) 
data <- raster("/Users/guoyafei/Downloads/share/spatial03/worldclim/cmip6/7_fut/10m/BCC-CSM2-MR/ssp585/wc2.1_10m_bioc_BCC-CSM2-MR_ssp585_2021-2040.tif")

taxa <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/225.taxa.txt", sep="\t",header=T, stringsAsFactors = F)
lon <- as.numeric(taxa$Logititude)
lat <- as.numeric(taxa$Latitude)
samples <- data.frame(lon,lat)
temp.data <- samples
extract(data, samples)
temp.data$temp1 <- extract(temp1, samples)
rownames(temp.data) <- taxa[,1]
all <- temp.data[,-c(1,2)]
write.table(all, "select_bio2.txt", sep="\t", row.names = T,quote=F)

colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))
oregon <- readOGR('.', 'Oregon_10N')
levelplot(temp6_2018, 
          margin=FALSE,# suppress marginal graphics
          ylim=c(-50,60),
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=-5:5, font=4)      # legend ticks and labels 
          ),    
          par.settings=list(
            axis.line=list(col='transparent') # suppress axes and legend outline
          ),
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=colr,                   # colour ramp
          at=seq(-10, 50)) 
#layer(sp.polygons(oregon, lwd=3))           # add oregon SPDF with latticeExtra::layer

########################################################### 环境变量的PC在地图上的分布 ##########################################################
setwd("/Users/guoyafei/Desktop/3type/分布")
library(raster)
library(rgdal)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
library(viridis) 
library(FactoMineR)
library(factoextra)
library(psych)
library(scatterplot3d)
library(lattice)

type1 <- c("Latitude","Longitude","solar2","temp1","temp8","temp10","Elevation")
type2 <- c("Latitude", "Longitude","solar1","solar3","temp2","temp9","prec7")
type3 <- c("Latitude", "Longitude","soil9","soil11","soil13","prec3","prec5","prec8","temp7")

all <- c("Latitude","Longitude","solar2","temp1","temp8","temp10","Elevation","solar1","solar3","temp2","temp9","prec7","soil9","soil11","soil13","prec3","prec5","prec8","temp7")
sub1 <- data[, which(colnames(data) %in% type1)]
sub2 <- data[, which(colnames(data) %in% type2)]
sub3 <- data[, which(colnames(data) %in% type3)]

name1 <-  sub1[!duplicated(sub1, fromLast=TRUE), c(3:7) ] 
name2 <-  sub2[!duplicated(sub2, fromLast=TRUE), c(3:7) ] 
name3 <-  sub3[!duplicated(sub3, fromLast=TRUE), c(3:7) ] 

pca1 <- prcomp(name3,scale. = TRUE)
df1 <- pca1$x
df1 <- as.data.frame(df1)
write.table(df1, "type3.PC.txt", quote=F, row.names = T)

########################################################### 合并tif环境变量 ##########################################################
#做个测试 66:/data2/yafei/polygenic/worldclim
#export LD_LIBRARY_PATH=/data1/home/yafei/008_Software/anaconda3/pkgs/aspera-connect-3.9.6-0/opt/aspera/connect/lib:$LD_LIBRARY_PATH
library(raster)
library(sp)
library(sf)

library(raster)
library(ncdf4)
library(rasterVis)
library(lattice)

#conda activate worldclim
#crop & mask 函数可以对tif文件进行剪切
pre_5 = raster("/data2/yafei/polygenic/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_05.tif")
pre_6 = raster("/data2/yafei/polygenic/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_06.tif")
pre_7 = raster("/data2/yafei/polygenic/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_07.tif")
pre_8 = raster("/data2/yafei/polygenic/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_08.tif")
pre_9 = raster("/data2/yafei/polygenic/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_09.tif")

bio12 <- raster("/data2/yafei/polygenic/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")
#北美地区小麦主要收获季节是6-8月
color_palette <- colorRampPalette(c("#e66101","#fdb863","#ffffbf","#b2abd2","#5e3c99","#c2a5cf","#bababa","#a6dba0","#92c5de","#d7dadb","#abd9e9","#2c7bb6"))
color_palette <- colorRampPalette(c("#e66101","#fdb863","#ffffbf","#b2abd2","#5e3c99","#c2a5cf","#bababa","#abd9e9","#2c7bb6"))

m1 <- addLayer(pre_6, pre_7,pre_8)
r <- stackApply(m1, indices = 1, fun = sum)
shape <- read_sf("/data2/yafei/polygenic/worldclim/NA_CEC_Eco_Level2/NA_CEC_Eco_Level2.shp")
#查看矢量/栅格数据的CRS是否一致
#crs(shp)
#crs(r)
#转换投影
shp_proj <- st_transform(shape, crs=crs(r))
#掩膜
out <- crop(r, c(-175,-50,10,85))
r_mask <- mask(out, shp_proj)
pdf("NA-6-8.pdf")
plot(r_mask, col=color_palette(100), legend=TRUE)
dev.off()

#欧洲地区小麦主要收获季节是7-9月
m1 <- addLayer(pre_7, pre_8, pre_9)
r <- stackApply(m1, indices = 1, fun = sum)
shape <- read_sf("/data2/yafei/polygenic/worldclim/NUTS_RG_10M_2021_3035.shp/NUTS_RG_10M_2021_3035.shp")
shp_proj <- st_transform(shape, crs=crs(r))
out <- crop(r, c(-30,50,30,75))
r_mask <- mask(out, shp_proj)
pdf("EU-7-9.pdf")
plot(r_mask, col=color_palette(100), legend=TRUE)
dev.off()

#中国主要是5月-7月
m1 <- addLayer(pre_5, pre_6,pre_7)
r <- stackApply(m1, indices = 1, fun = sum)
shape <- read_sf("/data2/yafei/polygenic/worldclim/CHA.shp/country.shp")
shp_proj <- st_transform(shape, crs=crs(r))
out <- crop(r, c(60,150,10,60))
r_mask <- mask(out, shp_proj)
pdf("CH-5-7.pdf")
plot(r_mask, col=color_palette(100), legend=TRUE)
dev.off()

#library(terra)
#计算季节降水量，3–5月为春季，6–8月为夏季，9–11月为秋季，12月至翌年2月为冬季
#SIndices = sort(rep(1:4,3))
#library(gimms)
#Pre_CHA = stack("pre_CHA.tif")
#SPre = monthlyComposite(Pre_CHA, fun=sum, indices=SIndices, cores=4)

########################################### 确定landrace样本地图分布, kmeans聚类，划分欧亚大陆6个区域样本 #############################################
library(gdata)
library(ggmap)
library(RColorBrewer)
library(cluster)
library(factoextra)
setwd("/Users/guoyafei/Desktop/样本分群")
#绘图:全部的
mp <- NULL
mapworld <- borders("world",colour = "gray90",fill="gray90") 
mp <- ggplot() + mapworld + ylim(-60,90)  + theme_classic()
color <- brewer.pal(8, "Dark2")
color2 <- brewer.pal(8, "Set1")
color <- c(color,color2)

data <- read.xls("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx",sheet=3)
data <- as.data.frame(data)
row.names(data) <- data$ID
sub2 <- data[which(data$region != "Sth_AM" & data$region != "Nth_AM" & data$region != "AF" ),]
sub <- data[which(data$region != "Sth_AM" & data$region != "Nth_AM" & data$region != "AF" ),c(2,3,89:96,263:264)]


#data <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_523landrace_withBioclimate.txt", header=T, stringsAsFactors = F)
#row.names(data) <- data$ID
#data <- na.omit(data[which(data$region != "Sth_AM" & data$region != "Nth_AM" & data$region != "AF" ),])

#根据82个环境变量进行PCA，只需要做一次
clim <- data[,c(4:85)]
#根据环境变量给样本做PCA
dt <- as.matrix(scale(clim))
rm1<-cor(dt)
rs1<- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values

#换算成标准差(Standard deviation);
Standard_deviation <- sqrt(val)
#计算方差贡献率和累积贡献率；
Proportion_of_Variance <- val/sum(val)
Cumulative_Proportion <- cumsum(Proportion_of_Variance)
#碎石图绘制;
#par(mar=c(6,6,2,2))
#plot(rs1$values,type="b",cex=2,cex.lab=2,cex.axis=2,lty=2,lwd=2,xlab = "PC",ylab="Eigenvalue (Principal Component Variance)")
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
U<-as.matrix(rs1$vectors)
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:82,sep="")
PC <- as.data.frame(PC[,1:10])
all <- cbind(data[,c(2:3)], PC[,1:10])
  
#kmeans聚类看一下，手动调整
df = scale(sub)
#聚类数量 vs. 总体平方和
fviz_nbclust(df, kmeans, method = "wss")
#聚类数量 vs. 差距统计
gap_stat <- clusGap(df,FUN = kmeans,nstart = 25,K.max = 25,B = 50)
fviz_gap_stat(gap_stat)

set.seed(1)
km <- kmeans(df, centers = 14, nstart = 25)
fviz_cluster(km, data = df)
aggregate(USArrests, by=list(cluster=km$cluster), mean)
final_data <- cbind(rownames(sub),sub, cluster = km$cluster)

final_data$cluster <- as.factor(final_data$cluster)

mp+geom_point(data =final_data, aes(x=lon, y=lat,color=cluster),size=1)+
  scale_size(range=c(1,1)) +
  #scale_color_manual(values = color) +
  theme(panel.border = element_blank())


data <- read.table("82clim.txt", header=T,stringsAsFactors = F)
row.names(data) <- data$rownames.all.
PC$ID <- row.names(PC)
data$ID <- rownames(data)
all <- merge(PC,data,by="ID",all.x=T)
all$cluster <- as.factor(all$cluster)
mp+geom_point(data =sub, aes(x=lon, y=lat,color=cluster),size=1)+
  scale_size(range=c(1,1)) + 
  scale_color_manual(values = color) +
  theme(panel.border = element_blank())
sub <- all[which(all$cluster != "0"),]
ggplot()+geom_point(data =all, aes(x=PC1, y=PC2,color=cluster),size=1)+
  scale_size(range=c(1,1)) + 
  scale_color_manual(values = color) +
  theme(panel.border = element_blank())

#写进文档里，并手动根据kmeans聚类调整后保存为VMap3_479landracePC_6regions.txt，这个文件的PC是通过70个非土壤的环境变量做的
#write.table(final_data, "PC_cluster5.txt", quote=F, sep="\t", row.names=F)
#data <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_479landracePC_6regions.txt",header=T,stringsAsFactors = F)
#final_data$manual <- data[rownames(final_data),11]
#重新写进文档里，并手动调整后保存为VMap3_476landracePC_6regions.txt，这个文件的PC是通过82个环境变量做的
#write.table(final_data, "VMap3_476landracePC_6regions.txt",quote=F, sep="\t",row.names = F)

#根据PC结果画图
data <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_476landracePC_6regions.txt", header=T,stringsAsFactors = F)
data <- final_data[which(final_data$manual != "0"),]
data$manual <- as.factor(data$manual)
mp+geom_point(data =data, aes(x=lon, y=lat,color=manual),size=1)+
  scale_size(range=c(1,1)) + 
  theme(panel.border = element_blank())+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
  

#根据经纬度和环境PC进行PCA
data <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_398landrace_82clim_PC10.txt", header=T,stringsAsFactors = F)
sub <- data[,c(6:87)]
dt <- as.matrix(scale(sub))
rm1<-cor(dt)
rs1<- eigen(rm1)
#提取结果中的特征值，即各主成分的方差；
val <- rs1$values
#换算成标准差(Standard deviation);
Standard_deviation <- sqrt(val)
#计算方差贡献率和累积贡献率；
Proportion_of_Variance <- val/sum(val)
Cumulative_Proportion <- cumsum(Proportion_of_Variance)
#提取结果中的特征向量(也称为Loadings,载荷矩阵)；
U<-as.matrix(rs1$vectors)
#进行矩阵乘法，获得PC score；
PC <- dt %*% U
colnames(PC) <- paste("PC",1:82,sep="")
PC <- as.data.frame(PC)

PC$region.manual. <- as.factor(data$region.manual.)
PC$lon <- data$lon
PC$lat <- data$lat

mp+geom_point(data =PC, aes(x=lon, y=lat,color=region.manual.),size=2)+
  scale_size(range=c(1,1)) + 
  theme(panel.border = element_blank())+
  theme_classic()+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
#PC1 (32.37%)；PC2 (24.12%)
ggplot()+geom_point(data =PC, aes(x=PC1, y=PC2,color=region.manual.),size=3)+
  scale_size(range=c(1,1)) + 
  xlab("PC1 (32.37%)") +
  ylab("PC2 (24.12%)") +
  theme(panel.border = element_blank())+
  theme_classic()+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

###################################################### 环境PC在地图上的映射 ###########################################
data <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_398landrace_82clim_PC10.txt", header=T,stringsAsFactors = F)
data$region.manual. <- as.factor(data$region.manual.)
mp <- NULL
mapworld <- borders("world",colour = "gray90",fill="gray90") 
mp <- ggplot() + mapworld + ylim(-60,90)  + theme_classic()
mp+geom_point(data =data, aes(x=lon, y=lat,color=PC1),size=1.5)+
  scale_size(range=c(1,1)) + 
  theme(panel.border = element_blank())+
  theme_classic()+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

###################################################### 环境变量RDA分析 ################################################
setwd("/Users/guoyafei/Desktop/RDA")
data <- read.table("out2.txt", header= T, stringsAsFactors = F)
type <- c(rep("solar",time=15), rep("temp", time=12), rep("prec",time=24),rep("temp", time=11), rep("prec",time=8), rep("soil",time=12),rep("allPC",time=10), rep("solarPC",time=10), rep("tempPC",time=10), rep("precPC",time=10), rep("soilPC",time=12), rep("STR",time=10))
group <- c(rep("solar",time=15), rep("temp", time=12), rep("prec",time=24),rep("temp", time=11), rep("prec",time=8), rep("soil",time=12),rep("allPC",time=10), rep("solar",time=10), rep("temp",time=10), rep("prec",time=10), rep("soil",time=12), rep("STR",time=10))
data$type <- type
data$group <- group
solar <- data[which(data$group == "soil"),]
ggplot(solar, aes(x=a, y=b))+
  geom_col()+
  #facet_grid(type~.)+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90, hjust=1))




####################################### 样本根据环境变量进行聚类画图 #################
library(ggmap)
library(RColorBrewer)
library(cluster)
library(factoextra)
library(pheatmap)
library(corrgram)
library(GGally)
setwd("/Users/guoyafei/Desktop/climate")
data <- read.table("VMap3_sheet3_manual.txt", header =T, stringsAsFactors = F)

#经纬度,82个环境变量
clim <- data[,c(2,3,8:34,59:69,35:58,70:89)]

#经纬度,前三个环境
clim <- data[,c(2,3,172,173,174,187,188,189,210,211,212,242,243,244)]
p <- scale(clim)
p[43,81] <- 10
pheatmap(p,cluster_rows = F,clustering_distance_cols  = "correlation", border_color = "white",cutree_cols  = 3)
corrgram(clim, order = T,lower.panel = panel.shade, upper.panel = panel.shade,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")
ggcorr(clim, method = c("everything", "pearson"),nbreaks = 5,hjust = 1, size = 4, color = "grey50",layout.exp = 1)

data <- data[which(data$baypass_cluster_V1 !="0"),]
clim <- data[,c(2,3,8:89)]

#kmeans聚类看一下，手动调整
df = scale(clim)
#聚类数量 vs. 总体平方和
fviz_nbclust(df, kmeans, method = "wss")
#聚类数量 vs. 差距统计
gap_stat <- clusGap(df,FUN = kmeans,nstart = 25,K.max = 25,B = 50)
fviz_gap_stat(gap_stat)
set.seed(120)
km <- kmeans(df, centers = 14, nstart = 25)
fviz_cluster(km, data = df,repel = FALSE,geom = c("point"))
final_data <- cbind(rownames(clim),clim, cluster = km$cluster)
final_data$cluster <- as.factor(final_data$cluster)
pdf("plot17_divide.pdf",height = 14,width=20)
mp <- NULL
mapworld <- borders("world",colour = "gray90",fill="gray90") 
mp <- ggplot() + mapworld + ylim(-60,90)  + theme_classic()
mp+geom_point(data =final_data, aes(x=lon, y=lat,color = cluster),size=1)+
  scale_size(range=c(1,1)) + 
  theme(panel.border = element_blank())+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

######################################## 画不同群体的环境变量的分布 ##############
library(ggmap)
library(RColorBrewer)
library(cluster)
library(factoextra)
library(pheatmap)
library(corrgram)
library(GGally)
library(reshape2)
setwd("/Users/guoyafei/Desktop/climate")
data <- read.table("VMap3_sheet3_manual.txt", header =T, stringsAsFactors = F)
data <- data[which(data$baypass_cluster_V1 !="0"),]
#经纬度,82个环境变量
clim <- data[,c(7,2,3,8:34,59:69,35:58,70:89)]
clim2 <- scale(data[,c(2,3,8:34,59:69,35:58,70:89)])
clim3 <- as.data.frame(cbind(clim[,1],clim2))
colnames(clim3)[1] <- "baypass_cluster_V1"

#经纬度,前三个环境
clim3 <- data[,c(7,2,3,172,173,174,187,188,189,210,211,212,242,243,244)]

M <- aggregate(clim3, by=list(clim3$baypass_cluster_V1),FUN=mean)[,-1]
#按照纬度排序
a <- M[order(M$lat),]
a$baypass_cluster_V1 <- factor(a$baypass_cluster_V1,levels =a$baypass_cluster_V1)
cats <- melt(a,id="baypass_cluster_V1")
#按照经度排序
a <- M[order(M$lon ),]
a$baypass_cluster_V1 <- factor(a$baypass_cluster_V1,levels =a$baypass_cluster_V1)
cats <- melt(a,id="baypass_cluster_V1")

solar <- cats[which(cats$variable == "solar15_PC1" | cats$variable == "solar15_PC2"| cats$variable == "solar15_PC3" ),]
ggplot(solar, aes(x=baypass_cluster_V1, y=value, group= variable,shape= variable,color=variable)) + 
  geom_line() +geom_point(size=4) +theme_bw()+xlab("population")+
  theme(plot.title = element_text(color="red", size=15, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))

temp <- cats[which(cats$variable == "temp23_PC1" | cats$variable == "temp23_PC2"| cats$variable == "temp23_PC3" ),]
ggplot(temp, aes(x=baypass_cluster_V1, y=value, group= variable,shape= variable,color=variable)) + 
  geom_line() +geom_point(size=4) +theme_bw()+xlab("population")+
  theme(plot.title = element_text(color="red", size=15, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))

prec <- cats[which(cats$variable == "prec32_PC1" | cats$variable == "prec32_PC2"| cats$variable == "prec32_PC3"),]
ggplot(prec, aes(x=baypass_cluster_V1, y=value, group= variable,shape= variable,color=variable)) + 
  geom_line() +geom_point(size=4) +theme_bw()+xlab("population")+
  theme(plot.title = element_text(color="red", size=15, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))

soil <- cats[which(cats$variable == "soil12_PC1" | cats$variable == "soil12_PC2"| cats$variable == "soil12_PC3"),]
ggplot(soil, aes(x=baypass_cluster_V1, y=value, group= variable,shape= variable,color=variable)) + 
  geom_line() +geom_point(size=4) +theme_bw()+xlab("population")+
  theme(plot.title = element_text(color="red", size=15, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))

lat <- cats[which(cats$variable == "lat"),]
ggplot(lat, aes(x=baypass_cluster_V1, y=value, group= variable,shape= variable,color=variable)) + 
  geom_line() +geom_point(size=4) +theme_bw()+xlab("population")+
  theme(plot.title = element_text(color="red", size=15, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))

lon <- cats[which(cats$variable == "lon"),]
ggplot(lon, aes(x=baypass_cluster_V1, y=value, group= variable,shape= variable,color=variable)) + 
  geom_line() +geom_point(size=4) +theme_bw()+xlab("population")+
  theme(plot.title = element_text(color="red", size=15, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))
