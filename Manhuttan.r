library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR")
data1 <- read.table("South_all_change.txt",sep="\t",header=T,stringsAsFactors = F) 
#data2 <- data1[-which(data1$p == "NaN"),]
data1$SNP <- c(1:330025)
#gwasR <- data2[,c(2,3,4,7)]
gwasR <- data1[,c(7,1,2,5)]
#gwasR$P <- 10^(-gwasR$MeanY)
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
highsnp <- c(213040,213071,213159,213222,213259,213516,213558,213599,213647,213664,213731,213737,213748,213837,213841,213876,213882,213902,213951,213970,213979,213993,214021,214034,214054,214061,214073,214079,214115,214149,214166,214173,214206,214211,329310,329469,329542,329652,329989)
manhattan(gwasR, annotateTop = T, highlight = highsnp, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"), suggestiveline=FALSE,genomewideline=F,logp=F, ylim=c(-2,25))
#manhattan(gwasR, main="Manhattan plot", ylim=c(0, 10), cex=0.6, cex.axis=0.9, col = c("blue","orange"), suggestiveline = F, genomewideline = F, chrlabs = c(1:21))
#qq(gwasR$P, main="Q-Q plot of GWAS p-value", xlim=c(0,7), ylim=c(0,12), pch=18, col = "blue4", cex=1.5, las=1)
#South_A
abline(h = 3.849, col="red", lwd=2, lty=2)
#South_B
abline(h = 3.73, col="red", lwd=2, lty=2)
#South_D
abline(h = 3.58, col="red", lwd=2, lty=2)
#North_A
abline(h = 3.76, col="red", lwd=2, lty=2)
#North_B
abline(h = 3.98, col="red", lwd=2, lty=2)
#North_D
abline(h = 3.61, col="red", lwd=2, lty=2)

data1 <- read.table("North_all_change.txt",sep="\t",header=T,stringsAsFactors = F) 
#data2 <- data1[-which(data1$p == "NaN"),]
data1$SNP <- c(1:332899)
#gwasR <- data2[,c(2,3,4,7)]
gwasR <- data1[,c(7,1,2,5)]
#gwasR$P <- 10^(-gwasR$MeanY)
colnames(gwasR) <- c("SNP", "CHR", "BP","P")

highsnp <- c(99120,99232,99317,99324,99381,99530,99629,99658,99668,99693,99720,99752,99767,99808,99844,99847,99860,99888,99900,99917,99919,99922,99926,99941,99956)
manhattan(gwasR, annotateTop = T, highlight = highsnp, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"), suggestiveline=FALSE,genomewideline=F,logp=F, ylim=c(-2,25))

#20210902 FuGWAS
setwd("/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/07_气孔导度数据/20210928/")
data <- read.table("stoma.mlm.txt",header=F,stringsAsFactors = F)
colnames(data) <- c("SNP", "CHR", "BP","P")
sub <- data[order(data$P),]
sub$logP <- -log10(sub$P)
sub2 <- sub[which(sub$logP>2),1:4]
#data <- data[-1,]
png("height1.png")
#manhattan(sub2, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"), suggestiveline=FALSE,genomewideline=F,logp=T, ylim=c(2,10))
manhattan(data, col = c("#b2df8a","#b2df8a","#a6cee3","#a6cee3","#fdbf6f","#fdbf6f"), highlight = highsnp,suggestiveline=FALSE,genomewideline=F,logp=T, ylim=c(1,7))
dev.off()

png("stoma1.png")
#manhattan(sub2, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"), suggestiveline=FALSE,genomewideline=F,logp=T, ylim=c(2,10))
manhattan(data, highlight = highsnp, col = c("#b2df8a","#b2df8a","#a6cee3","#a6cee3","#fdbf6f","#fdbf6f"), suggestiveline=FALSE,genomewideline=F,logp=T, ylim=c(1.5,7))
dev.off()

highsnp <- c("21-30861571","23-18781242","23-22706233")
manhattan(gwasR, annotateTop = T, highlight = highsnp, col = c("#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#33a02c","#b2df8a","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#1f78b4","#a6cee3","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f","#ff7f00","#fdbf6f"), suggestiveline=FALSE,genomewideline=F,logp=F, ylim=c(-2,25))


setwd("/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/07_气孔导度数据/20210928/")
data <- read.table("height_all.mlm.txt",header=F,stringsAsFactors = F)
colnames(data) <- c("SNP", "CHR", "BP","P")
qq_dat <- data.frame(obs=-log10(sort(data$P,decreasing=FALSE)),
                     exp=-log10( ppoints(length(data$P))))
 pd_qq <- ggplot(data=qq_dat,aes(exp,obs))+
  geom_point(alpha=0.7,color="#7F7F7FFF")+
  geom_abline(color="#D62728FF")+
  xlab("Expected -log10(P-value)")+
  ylab("Observed -log10(P-value)")+
  scale_x_continuous(limits = c(0,7))+
  scale_y_continuous(limits = c(0,7))+
  theme(
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    #axis.line = element_line(size=0.8,color="black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())

png("height_qq.png")
pd_qq
dev.off()

sub3 <- sub[which(sub$logP>5),]

data <- read.table("/Users/guoyafei/Desktop/相关性.txt", header=T, stringsAsFactors = F)
sub <- data[!is.na(data$株高.郑老师),]
sub2 <- sub[!is.na(sub$气孔导度),]
library("ggpubr")
ggscatter(data, x = "株高.郑老师", y = "气孔导度", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Plant height", ylab = "Stomatal conductance")
ggscatter(data, x = "ZX株高", y = "气孔导度", color = "orgCty",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tiller", ylab = "Stomatal conductance")
ggplot(data=data, aes(x=ZX株高, y=气孔导度 ))+
  geom_point(size=3,aes(x=ZX株高, y=气孔导度,color=orgCty ))+
  stat_smooth(method="lm")+
  xlab("plant height")+
  ylab("Stomatal conductance")+
  theme_classic()


library(gggplot2)
ggplot(data,aes(x = ZX株高,color=orgCty)) +
  geom_histogram(aes(x = ZX株高),stat="bin",binwidth=1, boundary = 0)+
  theme_classic()+
  xlab("plant height")



geom_text(aes(label=as.character(round(..density..,2))),stat="bin",binwidth=0.01,boundary = 0,vjust=-0.5)


