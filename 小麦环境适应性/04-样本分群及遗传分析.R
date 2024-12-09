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

#kmeans聚类看一下
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


################################################## IBD & IBE & PC ##################################
library(gdata)
library (geosphere)
library(reshape)
library(ggplot2)
setwd("/Users/guoyafei/Desktop/baypass/fst")
data <- read.table("/Users/guoyafei/Desktop/baypass/fst/sixRegion_fst.out.txt", header = T, stringsAsFactors = F)
all <- read.xls("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx",sheet=3,na.strings=c("NA","#DIV/0!"))
row.names(all) <- all$ID

#地理距离
sub <- all[!is.na(all$manual_conti),]
M <- aggregate(sub, by=list(sub$manual_region), FUN=mean)
sixRegion <- M[,c(1,3,4,90,91,92)]
geo_sub <- sixRegion[,2:3,drop=F]
result3 <- as.data.frame(as.matrix(distm(geo_sub, fun = distGeo)))
row.names(result3) <- paste(1:6,sep="")
colnames(result3) <- paste(1:6,sep="")
result3$id <- rownames(result3)
c <- melt(result3,id="id")[,c(2,1,3)]
d <- c[which(as.numeric(c$variable) < as.numeric(c$id)),]
adap <- cbind(data[which(data$type == "adap"),],d[,3])
neutral <- cbind(data[which(data$type == "neutral"),],d[,3])
alldata <- rbind(adap,neutral)
colnames(alldata)[6] <- "Dis"
ggplot(alldata, aes(Dis,weightedfst,color=type)) +
  geom_point( size=3)+
  theme_classic()
  
#环境距离
geo_sub <- sixRegion[,4:6,drop=F]
result3 <- as.data.frame(as.matrix(dist(geo_sub, method = "euclidean")))
row.names(result3) <- paste(1:6,sep="")
colnames(result3) <- paste(1:6,sep="")
result3$id <- rownames(result3)
c <- melt(result3,id="id")[,c(2,1,3)]
d <- c[which(as.numeric(c$variable) < as.numeric(c$id)),]
adap <- cbind(data[which(data$type == "adap"),],d[,3])
neutral <- cbind(data[which(data$type == "neutral"),],d[,3])
alldata <- rbind(adap,neutral)
colnames(alldata)[6] <- "Dis"
ggplot(alldata, aes(Dis,weightedfst,color=type)) +
  geom_point( size=3)+
  theme_classic()
  
#baypass
data <- read.table("/Users/guoyafei/Desktop/baypass/fst/baypass_fst.out.txt", header = T, stringsAsFactors = F)

sub <- all[which(all$baypass_cluster != "0"),]
M <- aggregate(sub, by=list(sub$baypass_cluster), FUN=mean)
sixRegion <- M[,c(267,3,4,172,173)]
geo_sub <- sixRegion[,2:3,drop=F]
result3 <- as.data.frame(as.matrix(distm(geo_sub, fun = distGeo)))
row.names(result3) <- paste(1:14,sep="")
colnames(result3) <- paste(1:14,sep="")
result3$id <- rownames(result3)
c <- melt(result3,id="id")[,c(2,1,3)]
d <- c[which(as.numeric(c$variable) < as.numeric(c$id)),]
adap <- cbind(data[which(data$type == "adap"),],d[,3])
neutral <- cbind(data[which(data$type == "neutral"),],d[,3])
alldata <- rbind(adap,neutral)
colnames(alldata)[6] <- "Dis"
ggplot(alldata, aes(Dis,weightedfst,color=type)) +
  geom_point( size=3)+
  theme_classic()

#环境距离
geo_sub <- sixRegion[,4:5,drop=F]
result3 <- as.data.frame(as.matrix(dist(geo_sub, method = "euclidean")))
row.names(result3) <- paste(1:14,sep="")
colnames(result3) <- paste(1:14,sep="")
result3$id <- rownames(result3)
c <- melt(result3,id="id")[,c(2,1,3)]
d <- c[which(as.numeric(c$variable) < as.numeric(c$id)),]
adap <- cbind(data[which(data$type == "adap"),],d[,3])
neutral <- cbind(data[which(data$type == "neutral"),],d[,3])
alldata <- rbind(adap,neutral)
colnames(alldata)[6] <- "Dis"
ggplot(alldata, aes(Dis,weightedfst,color=type)) +
  geom_point( size=3)+
  theme_classic()

PC <- read.table("neutral_ABD_pca.eigenvec", header=F,stringsAsFactors = F)
rownames(PC) <- PC$V1
PCdata <- cbind(sub[,c(1,5,6)], PC[sub$ID,3:7])
ggplot(PCdata, aes(V3,V4,color=manual_conti)) +
  geom_point( size=3)+
  theme_classic()

################################################## 经纬度和什么环境变量相关 ##################################

library(ggplot2)
setwd("/Users/guoyafei/Desktop/baypass")
data <- read.table("latlon.txt",header=T, row.names = 1, stringsAsFactors = F)
lat <- data[,1,drop=F]
lon <- data[,2,drop=F]
lat$type <- "lat"
lon$type <- "lon"
colnames(lat) <- c("value","type")
colnames(lon) <- c("value","type")
lat$climate <- row.names(lat)
lon$climate <- row.names(lon)
all <- rbind(lat, lon)

all$climate <- factor(all$climate,levels = c("solar","temp","perc","soil"))

ggplot(data=all, mapping=aes(x = type, y=value, fill=climate))+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  scale_fill_manual(values = c("#D55E00","#CC79A7","#56B4E9","#009E73","#F0E442","#0072B2", "#999999"))

################################################## 基因的选择次数 ##################################
setwd("/Users/guoyafei/Desktop/baypass")
data <- read.table("selectTime_500b_top05.txt", header=T, stringsAsFactors = F)
sub <- data[which(data$select_time != "1"),]

sub$select_time <- factor(sub$select_time,levels = c("2","3","4",">4"))

sub <- data
sub$select_time <- factor(sub$select_time,levels = c("2","3","4"))
sub$type <- factor(sub$type,levels = c("solar","temp","prec","soil"))
ggplot(data=sub, mapping=aes(x = gene_num, y=type, fill=select_time))+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))

################################################## 基因的适应性强度分类 ##################################
setwd("/Users/guoyafei/Desktop/baypass")
data <- read.table("class.txt", header=F, stringsAsFactors = F)
library(ggplot2)
sub <- data[which(data$V3 != "all" & data$V3 != "other"),]
sub$fraction <- "NA"
sub[which(sub$V2 == "all"),4] <- sub[which(sub$V2 == "all"),1] / sum(sub[which(sub$V2 == "all"),1])
sub[which(sub$V2 == "lat"),4] <- sub[which(sub$V2 == "lat"),1] / sum(sub[which(sub$V2 == "lat"),1])
sub[which(sub$V2 == "lon"),4] <- sub[which(sub$V2 == "lon"),1] / sum(sub[which(sub$V2 == "lon"),1])
sub[which(sub$V2 == "type1"),4] <- sub[which(sub$V2 == "type1"),1] / sum(sub[which(sub$V2 == "type1"),1])
sub[which(sub$V2 == "type2"),4] <- sub[which(sub$V2 == "type2"),1] / sum(sub[which(sub$V2 == "type2"),1])
sub[which(sub$V2 == "type3"),4] <- sub[which(sub$V2 == "type3"),1] / sum(sub[which(sub$V2 == "type3"),1])
sub[which(sub$V2 == "type4"),4] <- sub[which(sub$V2 == "type4"),1] / sum(sub[which(sub$V2 == "type4"),1])

all <- sub[which(sub$V2 == "lat"),]

all$ymax <- cumsum(all$fraction)
all$ymin <- c(0, head(all$ymax, n=-1))
all$labelPosition <- (all$ymax + all$ymin) / 2
all$label <- paste0(all$V3, "\n value: ", all$V1)
ggplot(all, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=V3)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")


sub$select_time <- factor(sub$select_time,levels = c("2","3","4"))
sub$type <- factor(sub$type,levels = c("solar","temp","prec","soil"))
ggplot(data=sub, mapping=aes(x = gene_num, y=type, fill=select_time))+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))

################################################## 环境变量解释经纬度相关位点的遗传方差情况 ##################################
setwd("~/Desktop/RDA")
data <- read.table("RDA_lat_lon_inter.txt",header=T,stringsAsFactors = F)
data <- data[which(data$climate != "all"),]
data$type <- factor(data$type,levels = c("lat","lon"))
data$climate <- factor(data$climate,levels = c("solar","temp","prec","soil"))
ggplot(data=data, mapping=aes(x = type, y=RDA, fill=climate))+
  #geom_bar(position="dodge")+
  geom_bar(position="dodge", stat="identity")+
  theme_classic()+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))

################################################## 经纬度效应值分布 ##################################
setwd("/Users/guoyafei/Desktop/baypass")
lat1 <- read.table("lat_type.bed", header=F, stringsAsFactors = F)
lat <- lat1[sample(c(1:dim(lat1)[1]),size=5000),]
pdf("lat_all_BF.pdf", height=4,width=8)
ggplot(data=lat1, mapping=aes(x = V4, y=V5,color=V8))+
  facet_grid(.~V8)+
  #geom_bar(position="dodge")+
  geom_point(alpha=0.5)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))
dev.off()

pdf("lon_all_BF.pdf", height=4,width=8)
lon1 <- read.table("lon_type.bed", header=F, stringsAsFactors = F)
lon <- lon1[sample(c(1:dim(lon1)[1]),size=5000),]
ggplot(data=lon1, mapping=aes(x = V4, y=V5, color = V8))+
  #geom_bar(position="dodge")+
  facet_grid(.~V8)+
  geom_point(alpha=0.5)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))
dev.off()

strong1 <- lat[which(lat$V4 >= 10 & lat$V4 < 15),]
pdf("lat_strong1_BF.pdf", height=4,width=8)
ggplot(data=strong1, mapping=aes(x = V4, y=V5, color = V8))+
  #geom_bar(position="dodge")+
  facet_grid(.~V8)+
  geom_point(alpha=0.5)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))
dev.off()

strong1 <- lon[which(lon$V4 >= 10 & lon$V4 < 15),]
pdf("lon_strong1_BF.pdf", height=4,width=8)
ggplot(data=strong1, mapping=aes(x = V4, y=V5, color = V8))+
  #geom_bar(position="dodge")+
  facet_grid(.~V8)+
  geom_point(alpha=0.5)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))
dev.off()

strong1 <- lat[which(lat$V4 >= 15 & lat$V4 < 20),]
pdf("lat_strong2_BF.pdf", height=4,width=8)
ggplot(data=strong1, mapping=aes(x = V4, y=V5, color = V8))+
  #geom_bar(position="dodge")+
  facet_grid(.~V8)+
  geom_point(alpha=0.5)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))
dev.off()

strong1 <- lon[which(lon$V4 >= 15 & lon$V4 < 20),]
pdf("lon_strong2_BF.pdf", height=4,width=8)
ggplot(data=strong1, mapping=aes(x = V4, y=V5, color = V8))+
  #geom_bar(position="dodge")+
  facet_grid(.~V8)+
  geom_point(alpha=0.5)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))
dev.off()

strong1 <- lat[which(lat$V4 >= 20),]
pdf("lat_strong3_BF.pdf", height=4,width=8)
ggplot(data=strong1, mapping=aes(x = V4, y=V5, color = V8))+
  #geom_bar(position="dodge")+
  facet_grid(.~V8)+
  geom_point(alpha=0.5)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))
dev.off()

strong1 <- lon[which(lon$V4 >= 20),]
pdf("lon_strong3_BF.pdf", height=4,width=8)
ggplot(data=strong1, mapping=aes(x = V4, y=V5, color = V8))+
  #geom_bar(position="dodge")+
  facet_grid(.~V8)+
  geom_point(alpha=0.5)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))
dev.off()


################################################## 经纬度效应值分布 + envGWAS ######################################
library(ggExtra)
setwd("/Users/guoyafei/Desktop/envgwas/fst")
lat1 <- read.table("lon_type_fst_allele2.bed", header=F, stringsAsFactors = F)

sub <- lat1[,14:28]
max_values <- apply(sub, 1, max,na.rm = TRUE)
min_values <- apply(sub, 1, min,na.rm = TRUE)
mean_values <- apply(sub, 1, mean,na.rm = TRUE)
median_values <- apply(sub, 1, median,na.rm = TRUE)
all <- cbind(lat1, max_values, min_values, mean_values, median_values)

strong1 <- all[sample(c(1:dim(all)[1]),size=8000),]
strong1 <- all[which(all$V7 > 3),]

strong1 <- all[which(all$V9 < 0.05),]
strong1 <- strong1[sample(c(1:dim(strong1)[1]),size=8000),]

type1 <- strong1[which(strong1$V10 == "type1"),]
type2 <- strong1[which(strong1$V10 == "type2"),]
type3 <- strong1[which(strong1$V10 == "type3"),]
type4 <- strong1[which(strong1$V10 == "type4"),]

strong1$x <- abs(strong1$V34)
#strong1$x <- abs(strong1$V5)
strong1$y <- abs(strong1$V8)

p1 <- ggplot(data=strong1, mapping=aes(x = x, y = y))+
  geom_point(alpha=0.1)+
  theme_classic()+
  #geom_density_2d(color = "blue", size = 1) + # 将散点图转换为密度图
  geom_smooth(method = "lm", color = "#a50f15",se = TRUE,  fill = "lightgray") +
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("Baypass Effect Size")+
  ylab("envGWAS Effect Size")+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))

ggMarginal(p1, type = "histogram", fill = "#a50f15", color = "black")
ggMarginal(p1, type = "histogram", fill = "#54278f", color = "black")

data <- read.table("~/Desktop/test.txt3",header=F,stringsAsFactors = F)
ggplot(data=data, mapping=aes(x = V1, y = V2))+
  geom_point(alpha=0.1)+
  theme_classic()+
  xlab("Allele 1 Effect Size")+
  ylab("Allele 2 Effect Size")+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 25))

################################################## 经纬度效应值分布 + envGWAS + fst ######################################
setwd("/Users/guoyafei/Desktop/envgwas/fst")
lat1 <- read.table("lat_type_fst.bed", header=F, stringsAsFactors = F)

sub <- lat1[,14:28]
max_values <- apply(sub, 1, max,na.rm = TRUE)
min_values <- apply(sub, 1, min,na.rm = TRUE)
mean_values <- apply(sub, 1, mean,na.rm = TRUE)
median_values <- apply(sub, 1, median,na.rm = TRUE)
all <- as.data.frame(cbind(lat1, max_values, min_values, mean_values, median_values))
all <- all[!is.na(all),]

strong1 <- all[sample(c(1:dim(all)[1]),size=8000),]
#strong1 <- lat1[which(lat1$V4 >= 10 & lat1$V4 < 15),]
#strong1 <- lat1[which(lat1$V4 >= 15 & lat1$V4 < 20),]
#strong1 <- lat1[which(lat1$V4 >= 20),]
strong1 <- all[which(all$V7 > 3),]

strong1 <- all[which(all$V9 < 0.05),]
strong1 <- strong1[sample(c(1:dim(strong1)[1]),size=8000),]

type1 <- strong1[which(strong1$V10 == "type1"),]
#BF_beta: lat: Adjusted R-squared:  0.3492   p-value: < 2.2e-16
#BF_beta: lon: Adjusted R-squared:  0.3284  
#envGWAS_beta: lat: Adjusted R-squared:  0.1749  
#envGWAS_beta: lon: Adjusted R-squared:  0.1208 
type2 <- strong1[which(strong1$V10 == "type2"),]
#BF_beta: lat: Adjusted R-squared:  0.3793  p-value: < 2.2e-16
#BF_beta: lon: Adjusted R-squared:  0.1068 
#envGWAS_beta: lat: Adjusted R-squared:  0.2014  
#envGWAS_beta: lon: Adjusted R-squared:  0.01333 
type3 <- strong1[which(strong1$V10 == "type3"),]
#BF_beta: lat: Adjusted R-squared:  0.2805 
#BF_beta: lon: Adjusted R-squared:  0.4471 
#envGWAS_beta: lat: Adjusted R-squared:  0.2067
#envGWAS_beta: lon: Adjusted R-squared:  0.1164 
type4 <- strong1[which(strong1$V10 == "type4"),]
#BF_beta: lat: Adjusted R-squared:  0.3528 
#BF_beta: lon: Adjusted R-squared:  0.2037 
#envGWAS_beta: lat: Adjusted R-squared:  0.403 
#envGWAS_beta: lon: Adjusted R-squared:  0.09374 

pdf("envGWAS_pop5-6.pdf", width=5, height=5)
ggplot(data=strong1, mapping=aes(x = abs(V8), y=V28))+
  #geom_bar(position="dodge")+
  #facet_grid(.~V10)+
  geom_point(alpha=0.3)+
  theme_classic()+
  #xlim(0,15)+
  #ylim(0,0.06)+
  xlab("envGWAS Effect Size")+
  ylab("Fst")+
  #geom_vline(xintercept=3,color='#E69F00',linetype = "dashed")+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  #geom_smooth(se = TRUE, method = "lm", formula = y ~ x)+
  scale_fill_manual(values = c("#dadaeb","#9e9ac8","#6a51a3","#3f007d"))+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 25))
dev.off()


