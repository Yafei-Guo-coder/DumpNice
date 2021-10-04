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



