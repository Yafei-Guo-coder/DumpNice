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

#方式一：(在计算环境距离/地理距离和遗传距离的相关性时使用），根据环境PC选变量的问题：环境PC最相关的变量，在基因组上相应的位点很少，并不能很好的反应环境适应性。除了prec是PC2，其余都是PC1.
#  temp11(0.98);temp1(0.96);temp6(0.96)(temp1鉴定的又多又好)
#  solar3(0.99);solar1(0.94)(solar1鉴定的又多又好)
#  prec1(0.98);prec5(0.85);prec6(0.85);prec3(0.83)
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
#更新后
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
library(CMplot)
setwd("/Users/guoyafei/Desktop/环境适应性位点/位置分布/fst")
mydata <- read.table("adap_fst.pos.21chr.withend.txt", header=F, stringsAsFactors = F)
colnames(mydata) <- c("chr","pos","snp")
data <- mydata[,c(3,1,2)]
head(data)
# snp         chr       pos
# snp1_1    1        2041
CMplot(data,plot.type="d",bin.size=1e2,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300)

setwd("/Users/guoyafei/Desktop/环境适应性位点/位置分布/all")
mydata <- read.table("soil_all.pos.21chr.withend.txt", header=F, stringsAsFactors = F)
colnames(mydata) <- c("chr","pos","snp")
data <- mydata[,c(3,1,2)]
head(data)
# snp         chr       pos
# snp1_1    1        2041
CMplot(data,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300)

setwd("/Users/guoyafei/Desktop/环境适应性位点/位置分布/fst_V2")
mydata <- read.table("soil_fst_V2.pos.21chr.withend.txt", header=F, stringsAsFactors = F)
colnames(mydata) <- c("chr","pos","snp")
data <- mydata[,c(3,1,2)]
head(data)
# snp         chr       pos
# snp1_1    1        2041
CMplot(data,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="pdf",memo="snp_density",dpi=300)

#曼哈顿图(lfmm & baypass)
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

#受选择位点的鉴定(RAiSD)
#方法一:去除着丝粒区后，取top千分之一作为阈值（小麦白粉病）
library(UpSetR)
#temp
temp <- list(
  CA <- c("pop8","pop12","pop11","pop5","pop25"),
  EU <- c("pop22","pop9","pop15","pop23"),
  SA <- c("pop4","pop10","pop19"),
  EA <- c("pop2","pop7","pop13","pop20","pop6"),
  WA <- c("pop3","pop17","pop1","pop16","pop18")
)
#prec
prec <- list(
  CA <- c("pop25","pop8","pop11","pop5","pop12"),
  EU <- c("pop23","pop9","pop15","pop22"),
  SA <- c("pop10","pop4","pop19"),
  EA <- c("pop2","pop7","pop6","pop13","pop20"),
  WA <- c("pop18","pop1","pop16","pop17","pop3")
)
#soil
soil <- list(
  CA <- c("pop25","pop5","pop8","pop11","pop12"),
  EU <- c("pop23","pop9","pop15","pop22"),
  SA <- c("pop10","pop4","pop19"),
  EA <- c("pop6","pop2","pop7","pop13","pop20"),
  WA <- c("pop1","pop17","pop16","pop18","pop3")
)
#solar
solar <- list(
  CA <- c("pop5","pop12","pop8","pop25","pop11"),
  EU <- c("pop15","pop9","pop23","pop22"),
  SA <- c("pop10","pop4","pop19"),
  EA <- c("pop2","pop7","pop6","pop20","pop13"),
  WA <- c("pop18","pop16","pop3","pop17","pop1")
)
all <- list(
  temp <- temp,
  prec <- prec,
  soil <- soil,
  solar <- solar
)

setwd("/Users/guoyafei/Desktop/选择/选择1")
var1 <- c("temp","prec","soil","solar")
var2 <- c("CA","EU","SA","EA","WA")
for(i in c(1:4)){
  for(j in c(1:5)){
    tmp <- all[[i]][[j]]
    infile <- paste(tmp, var1[i],"pos.bed", sep=".")
    a <- list()
    for(m in c(1:length(infile))){
      data <- read.table(infile[m], header=F,stringsAsFactors = F)
      a[[m]] <- paste(data[,1],data[,2],sep="-")
    }
    names(a) <- tmp
    outfile <- paste(var1[i],"_",var2[j],".pdf",sep="")
    pdf(outfile,height=6,width=10)
    p <- upset(fromList(a), keep.order = TRUE,sets = tmp, text.scale = c(2),point.size = 2.5, line.size = 1.5)
    print(p)
    dev.off()
  }
}

#画位点的分布：
#TraesCS7A02G037700 
#19:16835876-16837161
library(ggplot2)
library(gridExtra)
setwd("/Users/guoyafei/Desktop/选择/选择1/位点/")
file <- paste("pop",1:25,".bed3",sep="")
p <- list() 
for(i in c(2,6,7,13,20)){
  data <- read.table(file[i], header=F,stringsAsFactors = F)
  data$pos <- data$V2 - 13504417
  p[[i]] <- ggplot(data,aes(pos,V4))+
    theme_bw()+
    geom_vline(xintercept = 16835876-13504417-100000, color = 'red', size = 0.5) + 
    geom_vline(xintercept = 16836550-13504417, color = 'red', size = 0.5) + 
    geom_hline(yintercept = 33, color = 'gray', size = 0.5, linetype="dashed") +
    xlim(1,3500000)+
    ylab("mu")+
    ggtitle(paste("pop",i,sep=""))+
    geom_point(alpha=0.6,size=0.5)
}

pdf("EA.pdf", width = 3, height = 9)
grid.arrange(p[[2]],p[[6]],p[[7]],p[[13]],p[[20]],nrow=5)
dev.off()

pdf("pop.pdf", width = 15, height = 9)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],nrow=5)
dev.off()

#画lassip的影扫描和软扫描
#37:16835876-16837161
#37:16435876-17237161（基因上下游400k）
setwd("/Users/guoyafei/Desktop/选择/选择1/硬软扫描")
library(ggplot2)
library(gridExtra)
data <- read.table("/Users/guoyafei/Desktop/选择/选择1/硬软扫描/pop1.lassip.bed", header=F,stringsAsFactors = F)
file <- paste("pop",1:25,".lassip.bed",sep="")
p <- list() 
for(i in c(1:25)){
  data <- read.table(file[i], header=F,stringsAsFactors = F)
  data$pos <- data$V5 - 16450306
  point <- data[which(data$V9 == max(data$V9)), 10]
  
  point2 <- data[which(data$V9 == max(data$V9)), 8]
  print(point2)
  p[[i]] <- ggplot(data,aes(pos,V9))+
    theme_bw()+
    geom_vline(xintercept = 16835876-16450306-100000, color = 'red', size = 0.5) + 
    geom_vline(xintercept = 16836550-16450306-100000, color = 'red', size = 0.5) +
    geom_vline(xintercept = point, color = 'blue', size = 0.5) +
    ylab("T")+
    ggtitle(paste("pop",i,sep=""))+
    geom_line()
  p[[i+25]] <- ggplot(data,aes(pos,V8))+
    theme_bw()+
    geom_vline(xintercept = 16835876-16450306-100000, color = 'red', size = 0.5) + 
    geom_vline(xintercept = 16836550-16450306+100000, color = 'red', size = 0.5) + 
    geom_vline(xintercept = point, color = 'blue', size = 0.5) +
    ylab("m")+
    ggtitle(paste("pop",i,sep=""))+
    geom_line()
}
pdf("lassip.pdf", width = 15, height = 18)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[26]],p[[27]],p[[28]],p[[29]],p[[30]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[31]],p[[32]],p[[33]],p[[34]],p[[35]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[36]],p[[37]],p[[38]],p[[39]],p[[40]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],p[[41]],p[[42]],p[[43]],p[[44]],p[[45]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[46]],p[[47]],p[[48]],p[[49]],p[[50]],nrow=10)
dev.off()
