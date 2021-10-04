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


