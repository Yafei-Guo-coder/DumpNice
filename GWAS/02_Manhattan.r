library(qqman)
library(tidyverse)
#WA Top5 nlr----
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Manhattan/gene/")
data1 <- read.table("WA_gene_7chr.txt",sep="\t",header=F,stringsAsFactors = F) 
data1$SNP <- c(1:76)
gwasR <- data1[,c(6,1,2,4)]
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
# 1)计算chr长度
#chr_len <- gwasR %>% 
#  group_by(CHR) %>% 
#  summarise(chr_len=max(BP))
chr_len <- tibble(
   CHR=1:21,
   chr_len= c(594102056,689851870,495453186,780798557,801256715,651852609,750843639,830829764,615552423,744588157,673617499,509857067,709773743,713149757,566080677,618079260,720988478,473592718,736706236,750620385,638686055)
   )
# 2)计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 
# 3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(gwasR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)

#查看转化后的数据
head(Snp_pos,2)
#1) 准备X轴标签位置--在每条chr的中间
#X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
start <- chr_pos$total
stop <- c(594102056,1283953926,1779407112,2560205669,3361462384,4013314993,4764158632,5594988396,6210540819,6955128976,7628746475,8138603542,8848377285,9561527042,10127607719,10745686979,11466675457,11940268175,12676974411,13427594796,14066280851)
X_axis <- tibble(
   CHR = 1:21,
   center = (start+stop)/2)

#2）绘制“改良版”Manhattan图
p1 <- ggplot(Snp_pos, aes(x=BPcum, y=P)) +
  #设置点的大小，透明度
  geom_point( aes(color=as.factor(CHR)), size=1) +
  #设置颜色
  scale_color_manual(values = rep(c( "orange"), 21 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  #scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  #geom_hline(yintercept = c(6, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), size= 1.2, linetype = c("dotted", "twodash")) + 
  #设置主题
  theme_bw() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size = 10),
    axis.title.x=element_text(size = 10),
  )+
  ylab('XP-CLR ratio')+xlab('Chromosome Position')


#EU Top5 nlr----
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Manhattan/gene/")
data1 <- read.table("EU_gene_7chr.txt",sep="\t",header=F,stringsAsFactors = F) 
data1$SNP <- c(1:61)
gwasR <- data1[,c(6,1,2,4)]
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
# 1)计算chr长度
#chr_len <- gwasR %>% 
#  group_by(CHR) %>% 
#  summarise(chr_len=max(BP))
chr_len <- tibble(
  CHR=1:21,
  chr_len= c(594102056,689851870,495453186,780798557,801256715,651852609,750843639,830829764,615552423,744588157,673617499,509857067,709773743,713149757,566080677,618079260,720988478,473592718,736706236,750620385,638686055)
)
# 2)计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 
# 3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(gwasR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)

#查看转化后的数据
head(Snp_pos,2)
#1) 准备X轴标签位置--在每条chr的中间
#X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
start <- chr_pos$total
stop <- c(594102056,1283953926,1779407112,2560205669,3361462384,4013314993,4764158632,5594988396,6210540819,6955128976,7628746475,8138603542,8848377285,9561527042,10127607719,10745686979,11466675457,11940268175,12676974411,13427594796,14066280851)
X_axis <- tibble(
  CHR = 1:21,
  center = (start+stop)/2)

#2）绘制“改良版”Manhattan图
p2 <- ggplot(Snp_pos, aes(x=BPcum, y=P)) +
  #设置点的大小，透明度
  geom_point( aes(color=as.factor(CHR)), size=1) +
  #设置颜色
  scale_color_manual(values = rep(c( "orange"), 21 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  #scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  #geom_hline(yintercept = c(6, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), size= 1.2, linetype = c("dotted", "twodash")) + 
  #设置主题
  theme_bw() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size = 10),
    axis.title.x=element_text(size = 10),
  )+
  ylab('XP-CLR ratio')+xlab('Chromosome Position')

#SCA Top5 nlr----
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Manhattan/gene/")
data1 <- read.table("SCA_gene_7chr.txt",sep="\t",header=F,stringsAsFactors = F) 
data1$SNP <- c(1:88)
gwasR <- data1[,c(6,1,2,4)]
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
# 1)计算chr长度
#chr_len <- gwasR %>% 
#  group_by(CHR) %>% 
#  summarise(chr_len=max(BP))
chr_len <- tibble(
  CHR=1:21,
  chr_len= c(594102056,689851870,495453186,780798557,801256715,651852609,750843639,830829764,615552423,744588157,673617499,509857067,709773743,713149757,566080677,618079260,720988478,473592718,736706236,750620385,638686055)
)
# 2)计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 
# 3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(gwasR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)

#查看转化后的数据
head(Snp_pos,2)
#1) 准备X轴标签位置--在每条chr的中间
#X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
start <- chr_pos$total
stop <- c(594102056,1283953926,1779407112,2560205669,3361462384,4013314993,4764158632,5594988396,6210540819,6955128976,7628746475,8138603542,8848377285,9561527042,10127607719,10745686979,11466675457,11940268175,12676974411,13427594796,14066280851)
X_axis <- tibble(
  CHR = 1:21,
  center = (start+stop)/2)

#2）绘制“改良版”Manhattan图
p3 <- ggplot(Snp_pos, aes(x=BPcum, y=P)) +
  #设置点的大小，透明度
  geom_point( aes(color=as.factor(CHR)), size=1) +
  #设置颜色
  scale_color_manual(values = rep(c( "orange"), 21 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  #scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  #geom_hline(yintercept = c(6, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), size= 1.2, linetype = c("dotted", "twodash")) + 
  #设置主题
  theme_bw() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size = 10),
    axis.title.x=element_text(size = 10),
  )+
  ylab('XP-CLR ratio')+xlab('Chromosome Position')

#EA Top5 nlr----
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Manhattan/gene/")
data1 <- read.table("EA_gene_7chr.txt",sep="\t",header=F,stringsAsFactors = F) 
data1$SNP <- c(1:70)
gwasR <- data1[,c(6,1,2,4)]
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
# 1)计算chr长度
#chr_len <- gwasR %>% 
#  group_by(CHR) %>% 
#  summarise(chr_len=max(BP))
chr_len <- tibble(
  CHR=1:21,
  chr_len= c(594102056,689851870,495453186,780798557,801256715,651852609,750843639,830829764,615552423,744588157,673617499,509857067,709773743,713149757,566080677,618079260,720988478,473592718,736706236,750620385,638686055)
)
# 2)计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 
# 3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(gwasR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)

#查看转化后的数据
head(Snp_pos,2)
#1) 准备X轴标签位置--在每条chr的中间
#X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
start <- chr_pos$total
stop <- c(594102056,1283953926,1779407112,2560205669,3361462384,4013314993,4764158632,5594988396,6210540819,6955128976,7628746475,8138603542,8848377285,9561527042,10127607719,10745686979,11466675457,11940268175,12676974411,13427594796,14066280851)
X_axis <- tibble(
  CHR = 1:21,
  center = (start+stop)/2)

#2）绘制“改良版”Manhattan图
p4 <- ggplot(Snp_pos, aes(x=BPcum, y=P)) +
  #设置点的大小，透明度
  geom_point( aes(color=as.factor(CHR)), size=1) +
  #设置颜色
  scale_color_manual(values = rep(c( "orange"), 21 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  #scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  #geom_hline(yintercept = c(6, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), size= 1.2, linetype = c("dotted", "twodash")) + 
  #设置主题
  theme_bw() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size = 10),
    axis.title.x=element_text(size = 10),
  )+
  ylab('XP-CLR ratio')+xlab('Chromosome Position')

#South Top5 nlr----
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Manhattan/gene/")
data1 <- read.table("South_gene_7chr.txt",sep="\t",header=F,stringsAsFactors = F) 
data1$SNP <- c(1:89)
gwasR <- data1[,c(6,1,2,4)]
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
# 1)计算chr长度
#chr_len <- gwasR %>% 
#  group_by(CHR) %>% 
#  summarise(chr_len=max(BP))
chr_len <- tibble(
  CHR=1:21,
  chr_len= c(594102056,689851870,495453186,780798557,801256715,651852609,750843639,830829764,615552423,744588157,673617499,509857067,709773743,713149757,566080677,618079260,720988478,473592718,736706236,750620385,638686055)
)
# 2)计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 
# 3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(gwasR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)

#查看转化后的数据
head(Snp_pos,2)
#1) 准备X轴标签位置--在每条chr的中间
#X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
start <- chr_pos$total
stop <- c(594102056,1283953926,1779407112,2560205669,3361462384,4013314993,4764158632,5594988396,6210540819,6955128976,7628746475,8138603542,8848377285,9561527042,10127607719,10745686979,11466675457,11940268175,12676974411,13427594796,14066280851)
X_axis <- tibble(
  CHR = 1:21,
  center = (start+stop)/2)

#2）绘制“改良版”Manhattan图
p5 <- ggplot(Snp_pos, aes(x=BPcum, y=P)) +
  #设置点的大小，透明度
  geom_point( aes(color=as.factor(CHR)), size=1) +
  #设置颜色
  scale_color_manual(values = rep(c( "orange"), 21 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  #scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  #geom_hline(yintercept = c(6, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), size= 1.2, linetype = c("dotted", "twodash")) + 
  #设置主题
  theme_bw() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size = 10),
    axis.title.x=element_text(size = 10),
  )+
  ylab('XP-CLR ratio')+xlab('Chromosome Position')
#North Top5 nlr----
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Manhattan/gene/")
data1 <- read.table("North_gene_7chr.txt",sep="\t",header=F,stringsAsFactors = F) 
data1$SNP <- c(1:92)
gwasR <- data1[,c(6,1,2,4)]
colnames(gwasR) <- c("SNP", "CHR", "BP","P")
# 1)计算chr长度
#chr_len <- gwasR %>% 
#  group_by(CHR) %>% 
#  summarise(chr_len=max(BP))
chr_len <- tibble(
  CHR=1:21,
  chr_len= c(594102056,689851870,495453186,780798557,801256715,651852609,750843639,830829764,615552423,744588157,673617499,509857067,709773743,713149757,566080677,618079260,720988478,473592718,736706236,750620385,638686055)
)
# 2)计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 
# 3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(gwasR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)

#查看转化后的数据
head(Snp_pos,2)
#1) 准备X轴标签位置--在每条chr的中间
#X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
start <- chr_pos$total
stop <- c(594102056,1283953926,1779407112,2560205669,3361462384,4013314993,4764158632,5594988396,6210540819,6955128976,7628746475,8138603542,8848377285,9561527042,10127607719,10745686979,11466675457,11940268175,12676974411,13427594796,14066280851)
X_axis <- tibble(
  CHR = 1:21,
  center = (start+stop)/2)

#2）绘制“改良版”Manhattan图
p6 <- ggplot(Snp_pos, aes(x=BPcum, y=P)) +
  #设置点的大小，透明度
  geom_point( aes(color=as.factor(CHR)), size=1) +
  #设置颜色
  scale_color_manual(values = rep(c( "orange"), 21 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  #scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  #geom_hline(yintercept = c(6, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), size= 1.2, linetype = c("dotted", "twodash")) + 
  #设置主题
  theme_bw() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size = 10),
    axis.title.x=element_text(size = 10),
  )+
  ylab('XP-CLR ratio')+xlab('Chromosome Position')
library("cowplot")
plot_grid(p1,p2, p3,p4,p5,p6,
          labels = c("WA", "EU", "SCA","EA","S","N"),
          ncol = 1, nrow = 6,label_size = 6)

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


#20210902 FuGWAS----
#准备文件来自：01_GWAS_gwas.sh
#Manhattan plot
setwd("/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/07_气孔导度数据/20211007/Manhattan/logP2")
path <- "/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/07_气孔导度数据/20211007/Manhattan/logP2" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F)})
#标注气孔相关基因在曼哈顿上的位置
#shell:yafei@66:/data1/home/yafei/009_GWAS/WEGA_out/stoma/Manhattan/logP2
grep -f Gene_id.txt /data1/home/yafei/009_GWAS/gene/gene_v1.1_Lulab.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$4-1000000"\t"$5+1000000"\t"$6"\t"$7"\t"$8"\t"$9}' | sed '1i ##gff-version 3' > Related_gene_1M.gff3
grep -f Fu_known.gene /data1/home/yafei/009_GWAS/gene/gene_v1.1_Lulab.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$4-1000000"\t"$5+1000000"\t"$6"\t"$7"\t"$8"\t"$9}' | sed '1i ##gff-version 3' | sort -k1,1n -k4,4n > Fu_gene_1M.gff3

for i in `ls *txt`
do
awk '{print $2"\t"$3-1"\t"$3}' $i | sed '1d ' > ${i::-3}bed
done

for i in `ls *bed`; do bedtools intersect -a ../Related_gene_5k.gff3 -b $i -wb; done | awk '{print $10"\t"$12}'|sed 's/\t/-/' > 5k.snp
for i in `ls *bed`; do bedtools intersect -a ../Related_gene.gff3 -b $i -wb; done | awk '{print $10"\t"$12}'|sed 's/\t/-/' > snp
for i in `ls *bed`; do bedtools intersect -a ../Related_gene_1M.gff3 -b $i -wb; done| awk '{print $10"\t"$12}'|sed 's/\t/-/' > 1M.snp

#基因上下游5M的位点
snp <- read.table("5M.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]
#基因上下游1M的位点
snp <- read.table("1M.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]
#Fu已知基因上下游1M的位点
snp <- read.table("Fu.1M.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]
#基因上下游3M的位点
snp <- read.table("3M.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]
#基因上下游5k的位点
highsnp <- c("4-28222453","5-346668758","14-196323275","14-196324770","14-196327488","14-196329328","14-196329549","14-196330321","14-196330993","14-196331076","14-196331701","18-39729024","21-30861571")
#基因区的位点
highsnp <- c("14-196327488","21-30861571")
#-lop的值大于5的位置
snp <- read.table("logP5_50k.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]

pdf("stoma_logP5_high_50k.pdf",height = 5,width = 15)
for (i in c(1:20)){
  all <- data[[i]]
  colnames(all) <- c("SNP", "CHR", "BP","P")
  #sub <- data[order(data$P),]
  #sub$logP <- -log10(sub$P)
  #sub2 <- sub[which(sub$logP>2),1:4]
  manhattan(all, col = c("#fdbf6f","#BEBADA"), highlight = highsnp,suggestiveline=FALSE,genomewideline=F,logp=T, ylim=c(2,8))
  #manhattan(all, col = c("#fdbf6f","#fdbf6f"), suggestiveline=FALSE, genomewideline=F, logp=T, ylim=c(2,8))
}
dev.off()

#标注-logP的值在5以上的信号位点在曼哈顿上的位置
#shell:yafei@66:/data1/home/yafei/009_GWAS/WEGA_out/stoma/logP5

#提取-lopP5 bed 上下游50k区域的富集基因
#提取-lopP5 bed 上下游50k的区域
for i in `cat txt.names`
do 
awk '{print $2"\t"$3"\t"$4"\t"$7"\t"(-log($7)/log(10))}' $i | awk '{if($5>5 && $4!="NaN") print $0}' |awk '{print $2"\t"$3-50000"\t"$3+50000}' |sed '1d' > logP5/${i::-15}.50k.bed
done
#把-lopP5 bed 上下游50k的区域在1M之内的合并在一起
for i in {1..7}
do
for j in {A,B,D}
do
bedtools merge -d 1000000 -i ${i}${j}.50k.bed > ${i}${j}.50k.merge.bed  
done
done
for i in `ls *50k.merge.bed`
do
awk '{print $0"\t"FILENAME"\t"NR}' $i >> All_50k_region.bed
done
#曼哈顿图注释表
bedtools intersect -a /data1/home/yafei/009_GWAS/gene/gene_v1.1_Lulab.gff3 -b All_50k_region.bed -wb | awk 'split($9, array, ";") {print $1"\t"$4"\t"$5"\t"array[1]"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' | sed '1i gene_chr\tgene_start\tgene_end\tgene_id\tsnpBlock_chr\tsnpBlock_start\tsnpBlock_end\tfileName\tsnpBlock_id'> snpBlock_annotation.txt

#mlm bed
for i in `cat txt.names`
do 
awk '{if($4!="NaN")print $3"\t"$4-1"\t"$4}' $i | sed '1d'  > logP5/${i::-15}.mlm.bed 
done

#intersect
for i in {1..7}
do
for j in {A,B,D}
do
bedtools intersect -a ${i}${j}.50k.bed -b ${i}${j}.mlm.bed -wb
done
done | awk '{print $4"\t"$6}'|sort -k1,1n -k2,2n | uniq | sed 's/\t/-/' > logP5_50k.snp



#QQ plot
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

