library(corrgram)
data <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)
sub <- data[,c(3,6:40)]
sub2 <- sub[,c(1:9,11:14,22:36)]
sub3 <-  sub2[!duplicated(sub2, fromLast=TRUE), ] 
pheatmap(scale(sub3),cluster_rows = F,clustering_distance_cols  = "correlation",border_color = "white",cutree_cols  = 6)

df = scale(t(sub3),center = T,scale = T)
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

sub2 <- sub[,c(1,4,5,8,11,13,16,24,29,34,35,36)]
corrgram(sub3, order = F,lower.panel = panel.shade,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")

a <- hclust(sub3, method = "complete")

library(usdm)
vif(sub3)
v1 <- vifcor(sub3,th=0.3)

v1
#环境变量的筛选
library(ggplot2)
setwd("/Users/guoyafei/Desktop/baypass")
data <- read.table("471_baypass_taxa.Info", header=T,row.names = 1,stringsAsFactors = F)
env <- data[,c(3:4,6:40)]
#M <- aggregate(env, by=list(env$type),FUN=mean)
#solar <- M[,25:27]
#temp <- M[,c(2,28:38)]
#prec <- M[,4:11]
#soil <- M[,c(12:24)]

solar <- env[,24:26]
temp <- env[,c(1,27:37)]
prec <- env[,3:10]
#soil <- env[,c(11:23)]#提取所有13个变量作为土壤的变量
soil <- env[,c(12:15,23)]#只提取9-13作为土壤的相关变量
all <- list(solar,temp,prec,soil)

#PCA
library(FactoMineR)
library(factoextra)
library(psych)

#查看loading
#PC_temp <- principal(temp, nfactors=3,rotate="none")
#PC_solar <- principal(solar, nfactors=3,rotate="none")
#PC_prec <- principal(prec, nfactors=3,rotate="none")
#PC_soil <- principal(soil, nfactors=3,rotate="none")

#个体相关性：
#temp1(0.99);temp6(0.92);temp10(0.90);temp11(0.90)
#solar1(0.99);solar3(0.94)
#prec1(0.97);prec5(0.84);prec6(0.83);prec3(0.81);prec2(0.80)
#soil9(0.94);soil12(0.93)

#群体相关性：
#PC1:temp11(0.98);temp1(0.96);temp6(0.96) 
#solar3(0.99);solar1(0.94)
#prec1(0.98);prec5(0.85);prec6(0.85);prec3(0.83);prec2(0.80)
#soil9(0.94);soil12(0.94);soil13(0.93);soil11(-0.83)
#population regions
EA <- c(2,6,7,13,20)
CA <- c(5,8,11,12,25)
WA <- c(1,3,16,17,18)
SA <- c(4,10,19)
EU <- c(9,15,22,23)
AM <- c(21,24)
AF <- c(14)

out <- c("solar","temp","prec","soil")
for(i in c(1:4)){
  sub <- all[[i]]
  dt <- as.matrix(scale(sub))
  rm1<- cor(dt)
  rs1<- eigen(rm1)
  val <- rs1$values
  U <- as.matrix(rs1$vectors)
  PC <- dt %*% U
  PC <- as.data.frame(PC)
  colnames(PC)[1:3] <- paste("PC",1:3,sep="")
  PC$type <- data$type
  M <- aggregate(PC, by=list(PC$type),FUN=mean)
  PC$type <- factor(PC$type, levels=rownames(M[order(M$PC1),]),ordered = TRUE)
  PC$region <- NA
  PC[which(PC$type %in% EA),dim(PC)[2]] <- "EA"
  PC[which(PC$type %in% WA),dim(PC)[2]] <- "WA"
  PC[which(PC$type %in% SA),dim(PC)[2]] <- "SA"
  PC[which(PC$type %in% EU),dim(PC)[2]] <- "EU"
  PC[which(PC$type %in% CA),dim(PC)[2]] <- "CA"
  PC[which(PC$type %in% AM),dim(PC)[2]] <- "AM"
  PC[which(PC$type %in% AF),dim(PC)[2]] <- "AF"
  color <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#999999","#000000")
  #乌拉尔图: "#999999"（灰色）, 野生一粒:"#E69F00"(橘黄色), 栽培一粒:"#56B4E9"(明兰色), AABBDD:"#009E73"(橄榄绿), AABB:"#F0E442"（明黄）,"#0072B2"(天蓝), EA："#D55E00"(橘红色), SCA："#CC79A7"(皮粉色))
  outfile <- paste(out[i],".pdf",sep="")
  pdf(outfile, height=3, width=7)
  p <- ggplot(PC, aes(x=type,y=PC1,color=region,fill=region)) +
    geom_boxplot(alpha=0.8) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    #scale_colour_discrete(breaks = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#999999","#000000"), labels = c('EU','WA','CA','SA','EA','AF','AM'))+
    ylab(out[i]) +
    xlab("Group") +
    theme_bw()
  print(p)
  dev.off()
}

#write.table(PC,"/Users/guoyafei/Desktop/GF/等位基因频率变化/solar-PC.txt",quote=F,row.names = F)
#cor.test(temp$temp1,temp$PC1)

M <- aggregate(sub, by=list(sub$type),FUN=mean)
sub$type <- factor(sub$type, levels=rownames(M[order(M$PC1),]),ordered = TRUE)
#soil
sub$type <- factor(sub$type, levels=c("10", "19", "21", "6",  "4",  "15", "14", "23", "2",  "9",  "16", "18", "5",  "17", "8",  "1",  "25", "22", "11", "3",  "13", "20", "7","12"),ordered = TRUE)
ggplot(sub2, aes(x=sub2$type, y = sub2$PC1)) +
    geom_boxplot(color="blue",
                 fill="blue",
                 alpha=0.2) +
    ylab("soil$PC1") +
    xlab("Group")+
    theme_bw()

Proportion_of_Variance <- val/sum(val)
Cumulative_Proportion <- cumsum(Proportion_of_Variance)
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


df<-data.frame(PC, dt_all$ploidy, dt_all$Region)
library(ggplot2)

xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

p1<-ggplot(data = PC,aes(x=PC1,y=PC2))+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  theme_classic()

