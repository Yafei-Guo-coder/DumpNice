#fst-piR
library(qqman)
library(tidyverse)
library(ggplot2)
library(gdata)
require(gridExtra)
#path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/09_GWAS/LD01/Manhuttan" ##文件目录
#fileNames <- dir(path)  ##获取该路径下的文件名
#filePath <- sapply(fileNames, function(x){ 
#  paste(path,x,sep='/')})   ##生成读取文件路径
#data <- lapply(filePath, function(x){
#  read.table(x, header=T,stringsAsFactors = F)})

setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/fst-piR/input")
file <- c("cultivar-landrace.A","cultivar-landrace.B","cultivar-landrace.D","domemmer-wildemmer.A","domemmer-wildemmer.B","freethresh-domemmer.A","freethresh-domemmer.B","landrace-freethresh.A","landrace-freethresh.B","landrace-strangulata.D")
input <- paste(file,".select",sep="")
output <- paste(file,".pdf",sep="")
for (i in c(1:length(input))) {
  p <- list()
#logpi-ratio----
  data <- read.table(input[i],header=T,stringsAsFactors = F)
  data <- data[which(data$piR != "Inf"),]
  a <- sample(c(1:dim(data)[1]),10000)
  b <- data[a[order(a)],]
  thresh1 <- as.numeric(quantile(log(data$piR),0.975))
  thresh2 <- as.numeric(quantile(log(data$piR),0.025))
  #colnames(b)[c(1,2,5)] <- c("CHR", "BP","P")
  colnames(b)[c(1,2,7)] <- c("CHR", "BP","P")
  gwasResults <- b
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
  p[[1]] <- ggplot(don, aes(x=BPcum, y=log(P))) +
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
    scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
    # Custom the theme:
    theme_bw() +
    ylab("log(Pi ratio)")+
    #ylab("weiight_fst")+
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=12),
      axis.text.y=element_text(size=12),
      axis.title.y=element_text(size = 12),
      axis.title.x=element_text(size = 12),
    )+
    scale_y_continuous(limits = c(-5,5))+
    geom_hline(yintercept=thresh1,color='#E69F00',linetype = "dashed")+
    geom_hline(yintercept=thresh2,color='#E69F00',linetype = "dashed")
#fst----
  #data <- read.table("input/landrace-strangulata.D.select",header=T,stringsAsFactors = F)
  #data <- data[which(data$piR != "Inf"),]
  #a <- sample(c(1:dim(data)[1]),10000)
  b <- data[a[order(a)],]
  thresh1 <- as.numeric(quantile(data$weighted.fst,0.95))
  colnames(b)[c(1,2,5)] <- c("CHR", "BP","P")
  #colnames(b)[c(1,2,7)] <- c("CHR", "BP","P")
  gwasResults <- b
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
  p[[2]] <- ggplot(don, aes(x=BPcum, y=P)) +
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
    scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
    # Custom the theme:
    theme_bw() +
    #ylab("log(piR)")+
    ylab("Weiight Fst")+
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=12),
      axis.text.y=element_text(size=12),
      axis.title.y=element_text(size = 12),
      axis.title.x=element_text(size = 12),
    )+
    scale_y_continuous(limits = c(0,1))+
    geom_hline(yintercept=thresh1,color='#E69F00',linetype = "dashed")
  #geom_hline(yintercept=thresh2,color='#E69F00',linetype = "dashed")
  pdf(output[i],width = 7.3,height = 4)
  grid.arrange(p[[1]],p[[2]],nrow=2)
  dev.off()
}
#scatter plot----
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/fst-piR/input")
file <- c("cultivar-landrace.A","cultivar-landrace.B","cultivar-landrace.D","domemmer-wildemmer.A","domemmer-wildemmer.B","freethresh-domemmer.A","freethresh-domemmer.B","landrace-freethresh.A","landrace-freethresh.B","landrace-strangulata.D")
input <- paste(file,".select",sep="")
output <- paste(file,".scatter.pdf",sep="")
for (i in c(1:length(input))) {
  p <- list()
data <- read.table(input[i],header=T,stringsAsFactors = F)
data <- data[which(data$piR != "Inf"),]
thresh1 <- as.numeric(quantile(log(data$piR),0.975))
thresh2 <- as.numeric(quantile(log(data$piR),0.025))
thresh3 <- as.numeric(quantile(data$weighted.fst,0.95))
a <- sample(c(1:dim(data)[1]),20000)
b <- data[a[order(a)],]
sub1 <- b[which((b$weighted.fst > thresh3) & (log(b$piR) < thresh2)),]
sub2 <- b[which((b$weighted.fst > thresh3) & (log(b$piR) > thresh1)),]
p <- ggplot(b, aes(log(piR), weighted.fst)) +
  geom_point(size=1,alpha=0.5,color="grey")+
  theme_bw() +
  geom_hline(yintercept=thresh3,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=thresh1,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=thresh2,color='#E69F00',linetype = "dashed")+
  geom_point(data=sub1,size=1,alpha=0.5,color="darkred")+
  geom_point(data=sub2,size=1,alpha=0.5,color="darkred")+
  xlab("log(Pi ratio)")+
  ylab("weighted Fst")
pdf(output[i],width = 3.5, height = 2.6)
print(p)
dev.off()
}

