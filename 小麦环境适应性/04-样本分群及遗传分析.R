################################################## 样本环境PC：整体环境和单一环境 #######################################
all <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_523landrace_82clim_PC10.txt", check.names = F,header=T,stringsAsFactors = F)
row.names(all) <- all$ID
sub <- all[,c(7:21)]
#15个：光：[1] 0.4583169 0.8507076 0.9366573 0.9696977 0.9872101
sub <- all[,c(22:33,58:68)]
#23个：温：[1] 0.6423842 0.8660568 0.9289074 0.9804116 0.9875202
sub <- all[,c(34:57,69:76)]
#32个：水：[1] 0.4580800 0.6615094 0.8217867 0.8817675 0.9231409 
sub <- all[,c(77:88)]
#12个：土：[1] 0.3659615 0.5028298 0.6232972 0.7033572 0.7694534 0.8240760 0.8770495 0.9191473
sub <- all[,c(7:88)]
#82个：[1] 0.3083311 0.5390865 0.6820894 0.7801551 0.8068087 0.8301100 0.8515767 0.8680436 0.8831071 0.8971792
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
colnames(PC) <- paste("12soil_PC",1:12,sep="")
PC <- as.data.frame(PC)
PC$ID <- all$ID

data <- merge(all,PC,by="ID")
data$manual_region <- as.factor(data$manual_region)

merge <- data
all2 <- merge(merge,PC,by="ID")
all3 <- merge(all2,PC, by="ID")
all4 <- merge(all3, PC, by="ID")
all5 <- merge(all4,PC,by="ID")
#write.table(all5,"VMap3_523landrace_climPC.txt", quote=F, row.names = F,sep="\t")

data <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_523landrace_climPC.txt",check.names = F, header=T, stringsAsFactors = F)
data$manual_region <- as.factor(data$manual_region)
sub <- data[which(data$manual_region !="0"),]

ggplot()+geom_point(data =sub, aes(x=`12soil_PC2`, y=`12soil_PC3`,color=manual_region),size=3)+
  scale_size(range=c(1,1)) + 
  xlab("PC1 (32.37%)") +
  ylab("PC2 (24.12%)") +
  theme(panel.border = element_blank())+
  theme_classic()+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

data2$seq <- as.factor(data2$seq)
ggplot(data2, aes(x=seq, y=value, group = variable,color = variable)) +
  geom_point(size=2)+
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("PC") + ylab("Cumulative Proportion") +
  geom_hline(yintercept=0.9,color='#E69F00',linetype = "dashed")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

################################################## GF准备文件 #######################################
#kmeans区分样本区域
library(readxl)
library(ggmap)
library(RColorBrewer)
library(cluster)
library(factoextra)
data <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_523landrace_climPC.txt",check.names = F, header = T, stringsAsFactors = F)
rownames(data) <- data$ID
sub <- data[,c(2,3,89:98)]

#kmeans聚类看一下，手动调整
df = scale(sub)
#聚类数量 vs. 总体平方和
fviz_nbclust(df, kmeans, method = "wss")
#聚类数量 vs. 差距统计
gap_stat <- clusGap(df,FUN = kmeans,nstart = 25,K.max =25,B = 50)
fviz_gap_stat(gap_stat)

set.seed(1)
km <- kmeans(df, centers = 14, nstart = 25)
fviz_cluster(km, data = df)
#aggregate(USArrests, by=list(cluster=km$cluster), mean)
final_data <- cbind(data[,c(2:3)], cluster = km$cluster)

final_data$cluster <- as.factor(final_data$cluster)
mp <- NULL
mapworld <- borders("world",colour = "gray90",fill="gray90") 
mp <- ggplot() + mapworld + ylim(-60,90)  + theme_classic()
mp+geom_point(data =final_data, aes(x=lon, y=lat,color=cluster),size=1)+
  scale_size(range=c(1,1)) + 
  #scale_color_manual(values = color) +
  theme(panel.border = element_blank())

sub <- final_data[which(final_data$cluster == "15"),]
mp+geom_point(data =sub, aes(x=lon, y=lat,color=cluster),size=1)+
  scale_size(range=c(1,1)) + 
  #scale_color_manual(values = color) +
  theme(panel.border = element_blank())

colnames(data[,c(2,3,6,89:93,171:175,186:190,209:213,241:245)])

################################################## 样本群体结构分布 ############################
setwd("/Users/guoyafei/Desktop/样本分群/群体结构")
library(ggplot2)
library(gdata)
all <- read.xls("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx",sheet=3,na.strings=c("NA","#DIV/0!"))
all <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3_382landrace_82clim_PC10.txt", check.names = F,header=T,stringsAsFactors = F)
row.names(all) <- all$ID
sub <- all[!is.na(all$manual_conti),]

sub$manual_conti <- as.factor(sub$manual_conti)
ggplot(sub, aes(x=str_PC1_AB, y=str_PC2_AB, group = manual_conti, color = manual_conti)) +
  geom_point(size=2)+
  theme_classic()+
  xlab("PC1") + ylab("PC2") +
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

mp <- NULL
mapworld <- borders("world",colour = "gray90",fill="gray90") 
mp <- ggplot() + mapworld + ylim(-60,90)  + theme_classic()
mp+geom_point(data =sub, aes(x=lon, y=lat,color=manual_conti),size=1.5)+
  scale_size(range=c(1,1)) + 
  theme(panel.border = element_blank())+
  theme_classic()+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

################################################## RDA分析鉴定适应性位点 ################
library(ggplot2)
library(robustbase)
library(robust)
library(qvalue)
library(vegan)
library(readxl)
library(ggpubr)
library(data.table)

rdadapt<-function(rda,K)
{
  loadings<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

all <- read.xls("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx",sheet=3)
all <- as.data.frame(all)
row.names(all) <- all$ID

env <- all[,c(186:188)]
prec <- all[,c(209:211)]
solar <- all[,c(171:173)]
soil <- all[,c(241:243)]
bio <- all[,c(89:91)]

data <- read.table("~/Desktop/RDA/ABD_30k.012.geno", header=F,stringsAsFactors=F)
rownames(data) <- all$ID
genoA <- data[,3:30002]

geno <- genoA[,sample(c(1:30000),size=5000)]

rda_all <- rda(geno~., env, scale = FALSE)
res_rdadapt <- rdadapt(rda_all, 3)
row.names(res_rdadapt) <- colnames(geno)

p1<- ggplot() +
  geom_line(aes(x=c(1:length(rda_all$CCA$eig)), y=as.vector(rda_all$CCA$eig)), linetype="dotted", size = 1.5, color="darkgrey") +
  geom_point(aes(x=c(1:length(rda_all$CCA$eig)), y=as.vector(rda_all$CCA$eig)), size = 3, color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) +
  ylab("Inertia") +
  theme_bw()

#which(qvalue(res_rdadapt[,1], fdr = 0.05)$signif)
p2<- ggplot() +
  geom_point(aes(x=rda_all$CCA$v[,1], y=rda_all$CCA$v[,2]), col = "gray86") +
  geom_point(aes(x=rda_all$CCA$v[which(res_rdadapt[,2] < 0.05),1], y=rda_all$CCA$v[which(res_rdadapt[,2] < 0.05),2]), col = "orange") +
  geom_segment(aes(xend=rda_all$CCA$biplot[,1]/10, yend=rda_all$CCA$biplot[,2]/10, x=0, y=0), colour="black", size=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.2*rda_all$CCA$biplot[,1]/10, y=1.2*rda_all$CCA$biplot[,2]/10, label = colnames(env))) +
  xlab("RDA 1") + ylab("RDA 2") +
  theme_bw() +
  theme(legend.position="none")

p3<- ggplot() +
  geom_point(aes(x=c(1:length(res_rdadapt[,1])), y=-log10(res_rdadapt[,1])), col = "gray83") +
  geom_point(aes(x=c(1:length(res_rdadapt[,1]))[which(res_rdadapt[,2] < 0.1)], y=-log10(res_rdadapt[which(res_rdadapt[,2] < 0.1),1])), col = "orange") +
  xlab("SNPs") + ylab("-log10(p.values)") +
  theme_bw()

#筛选显著位点
rda_sig_snp<- subset(res_rdadapt, res_rdadapt$q.values < 0.05)
write.table(rda_sig_snp,"rda_sig_snp.txt",row.names = T)
pdf("test2.pdf")
ggarrange(ggarrange(p1, p2, ncol = 2, labels = c("A", "B")), ggarrange(p3,ncol = 1, labels = c("C")), nrow = 2)
dev.off()

