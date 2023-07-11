setwd("/Users/guoyafei/Desktop/环境适应性位点")
#整体fst分布箱线图----
library(ggplot2)
library(tidyverse)
library(viridis)
name <- c("shuf","all","prec","soil","solar","temp")
file <- paste(name,c("_fst.txt","PC1V3_fst.txt","PC1V2_fst.txt","PC1V3_fst.txt","PC1_fst.txt","PC1_fst.txt"),sep="")
all <- data.frame(type="", id1="", id2="", meanfst="", weightedfst="", stringsAsFactors=FALSE)
for ( i in c(1:6)) {
  data <- read.table(file[i], header=T,stringsAsFactors = F)
  #meanfst
  data[which(data$meanfst < 0),4] <- 0
  all <- rbind(all,data)
}
sub <- as.data.frame(all[-1,])
sub$meanfst <- as.numeric(sub$meanfst)
sub$type <- factor(sub$type, levels=c("shuf","all","prec","temp","soil","solar"),ordered = TRUE)
ggplot(sub, aes(x=type, y = meanfst)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  #theme_ipsum() +
  ylab("meanfst") +
  xlab("Group")+
  theme_bw()

shuf <- sub[which(sub$type=="shuf"),4]
all <- sub[which(sub$type=="all"),4]
temp <- sub[which(sub$type=="temp"),4]
soil <- sub[which(sub$type=="soil"),4]
solar <- sub[which(sub$type=="solar"),4]
prec <- sub[which(sub$type=="prec"),4]

t.test(shuf,all)
t.test(shuf,temp)
#0.003
t.test(shuf,prec)
#0.014
t.test(shuf,soil)
t.test(shuf,solar)
t.test(all,prec)

#群体间fst分布热图----
library(corrgram)
library(reshape2)
library(tidyr)
setwd("/Users/guoyafei/Desktop/GF/等位基因频率变化/")
library("corrplot")
#data <- read.table("both_solar1.fst.txt", header=T,stringsAsFactors = F)
data <- read.table("shuf_fst.txt", header=T,stringsAsFactors = F)
data <- read.table("solarPC1_fst.txt", header=T,stringsAsFactors = F)
data2 <- data[,c(1,3,2,4,5)]
colnames(data2) <- colnames(data)
#meanfst
all <- rbind(data,data2)[,c(2,3,4)]
duijiao <- read.table("对角线.txt",header=T,stringsAsFactors = F)
all <- rbind(all,duijiao)
all[which(all$value < 0),3] <- 0
colnames(all) <- c("id1","variable","value")
cats <- cast(all,id1~variable,value.var="value",fun.aggregate = mean)

cats <- as.matrix(cats)
order <- cats[row.names(present[order(present$bio_solar1),]),row.names(present[order(present$bio_solar1),])]
#pdf("weight_A_fst.pdf",width=30,height=30)
corrplot(order,method = "color",col.lim = c(0, 0.5),type = 'upper',tl.col="black",tl.srt = 45,addrect=1,addCoef.col = "grey",number.cex=0.5,number.digits=2,tl.cex=1,cl.cex=1,cl.lim = c(0, 1))
#dev.off()

#等位基因频率变化----
setwd("/Users/guoyafei/Desktop/GF/等位基因频率变化/")
gfData <- read.table("GF.solar1_36bio.txt", row.names = "pop",header=T,stringsAsFactors = F)

candidate <- gfData[,grep("both",names(gfData))]
reference <- gfData [,grep("ref",names(gfData))]
baypass <- gfData [,grep("baypass",names(gfData))]
lfmm <- gfData [,grep("lfmm",names(gfData))]
present <- gfData[,c(1,2,grep("bio",names(gfData)))]
#bioclimatic <- paste("bio_",1:36,sep = "") 
bioclimatic <- colnames(gfData)[3:38]

candidate$solar1 <- present[,25]

p <- vector()
e <- vector()
for (i in c(1:1261)) {
  a <- cor.test(candidate[,i], candidate[,1263])
  p <- append(p,a$p.value)
  e <- append(e,a$estimate)
}

ggplot(candidate, aes(y = candidate$both21.263538003,x = solar1))+
  geom_point()

ggplot(candidate, aes(y = candidate$both21.263538003,x = solar1)) +
  geom_point(size=0.5,alpha=0.3,color="#999999")+
  theme_bw()+
  geom_smooth(method = "lm", se=TRUE, 
              color="#0072B2", formula = y ~ x) 

#环境距离和遗传距离的相关性-----
#群体相关性
setwd("/Users/guoyafei/Desktop/GF/等位基因频率变化")

方式一：根据环境PC选变量的问题：环境PC最相关的变量，在基因组上相应的位点很少，并不能很好的反应环境适应性。除了prec是PC2，其余都是PC1.
  temp11(0.98);temp1(0.96);temp6(0.96)(temp1鉴定的又多又好)
  solar3(0.99);solar1(0.94)(solar1鉴定的又多又好)
  prec1(0.98);prec5(0.85);prec6(0.85);prec3(0.83)
    更新成prec3 prec4 prec6上面prec1;prec5鉴定的位点都很少
  soil9(0.94);soil12(0.94);soil13(0.93);soil11(-0.83);soil10(soil12鉴定的又多又好)
    更新成soil12(0.94)上面soil9，soil13，soil11鉴定位点情况都不好,位点也比较少
方式二：根据基因组贡献程度以及环境之间的相关性(小于0.3)，选择环境变量(V4)
  temp8;temp9;temp4,temp5,temp2,elevation,temp7,temp3
  solar1;solar2
  prec7;prec8;prec4,prec3,prec6
  soil8;soil7;soil10;soil9;soil6,soil12;soil3

#环境变量的筛选
library(ggplot2)
setwd("/Users/guoyafei/Desktop/GF/等位基因频率变化")

#环境PCA
data <- read.table("471_baypass_taxa.Info", header=T,row.names = 1,stringsAsFactors = F)
env <- data[,c(3:4,6:40)]
M <- aggregate(env, by=list(env$type),FUN=mean)
#solar <- M[,25:27]
solar <- M[,25:26]
#temp <- M[,c(2,28:38)]
temp <- M[,c(2,31:34,36:38)]
#prec <- M[,4:11]
prec <- M[,c(6,7,9,10,11)]
#soil <- M[,c(13:16,24)]
soil <- M[,c(12:13,15:18,21:24)]

#查看loading
PC_temp <- principal(temp, nfactors=3,rotate="none")
PC_solar <- principal(solar, nfactors=2,rotate="none")
PC_prec <- principal(prec, nfactors=3,rotate="none")
PC_soil <- principal(soil, nfactors=3,rotate="none")

sub <- prec
dt <- as.matrix(scale(sub))
rm1<- cor(dt)
rs1<- eigen(rm1)
val <- rs1$values
U <- as.matrix(rs1$vectors)
PC <- dt %*% U
PC <- as.data.frame(PC)
colnames(PC)[1:2] <- paste("PC",1:2,sep="")
write.table(PC,"/Users/guoyafei/Desktop/GF/等位基因频率变化/prec-sub-PC.txt",quote=F,row.names = F)


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
#更新后：
all <- c(4,8,9,11,17,27,29,30,32,37)
#prec
#prec <- c(4,6,8,10,11)
#8,9,11 
prec <- c(4,8,9,11)
#soil
#soil <- c(4,15,16,17,18,26)
#去掉10
#soil <- c(4,16,17,18,26)
#15,17,19
soil <- c(4,17)
#solar
solar <- c(4,27,29)
#temp
temp <- c(4,30,32,37)
seq <- list(all,prec,soil,solar,temp)

env <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)
name <- c("all","prec","soil","solar","temp")
file <- paste(name,c("PC1V3_fst.txt","PC1V2_fst.txt","PC1V3_fst.txt","PC1_fst.txt","PC1_fst.txt"),sep="")
file_out <- paste(name,".pdf",sep="")
for ( i in c(1:5)) {
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


#适应性位点在全基因组上的物理分布
#install.packages('RIdeogram')
require(RIdeogram)
setwd("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/13_Plots/01_Density")
#gene_density <- read.table("ArinaLrFor_LTR_1.txt", header=T,stringsAsFactors = F)
gene_density2 <- read.table("gene_density.txt", header=T, stringsAsFactors = F)
gene_density <- read.table("21-1M_VMap3_SnpDensity.txt", header=T, stringsAsFactors = F)
wheat_karyotype <- read.table("wheat_karyotype.txt", header=T, stringsAsFactors = F)
ideogram(karyotype = wheat_karyotype, overlaid = gene_density2)
convertSVG("chromosome.svg", device = "pdf")

library(CMplot)
setwd("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据")
mydata<-read.table("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据/fastcall2_001_pos.txt",header=TRUE,sep="\t")
head(mydata)
# snp chr pos
# snp1_1  1 2041
# snp1_2  1 2062
# snp1_3  1 2190
CMplot(mydata,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300) 
mydata<-read.table("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据/fastcall_001_pos.txt",header=TRUE,sep="\t")
head(mydata)
# snp         chr       pos
# snp1_1    1        2041
# snp1_2    1        2062
# snp1_3    1        2190
CMplot(mydata,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300) 


#lineage density
A <- read.table("Alineage.txt",header=F,stringsAsFactors = F)
ggplot(A, mapping = aes(x = V4)) +
  geom_density( alpha = 0.5, color = "#999999",fill="#999999") +
  #geom_density(fill = "blue",alpha = 0.3,color="blue")+
  theme_bw()


B <- read.table("Blineage.txt",header=F,stringsAsFactors = F)
ggplot(B, mapping = aes(x = V4)) +
  geom_density( alpha = 0.5, color = "#377EB8",fill="#377EB8") +
  #geom_density(fill = "blue",alpha = 0.3,color="blue")+
  theme_bw()

D <- read.table("Dlineage.txt",header=F,stringsAsFactors = F)
ggplot(D, mapping = aes(x = V4)) +
  geom_density( alpha = 0.5, color = "#FF7F00",fill="#FF7F00") +
  #geom_density(fill = "blue",alpha = 0.3,color="blue")+
  theme_bw()

ggplot(D, mapping = aes(x = V4),color = "blue",fill="blue") +
  geom_line(stat = "density")
