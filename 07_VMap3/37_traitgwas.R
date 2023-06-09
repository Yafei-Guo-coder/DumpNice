#-----qq----
library(ggplot2)
library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Desktop/traitgwas/pdf")
path <- "/Users/guoyafei/Desktop/traitgwas/plot/tasselgwas/cultivar"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F)})
p <- list()
q <- list()
data <- read.table("cultivar.test.txt",header=F,stringsAsFactors = F)
colnames(data) <- c("SNP","CHR","BP","P")
for ( i in c(1:7)){
  #qq plot
  n1 <- strsplit( names(data)[i], ".shuf50k.")[[1]][1]
  filename1 <- paste(n1,"cultivar.tassel.qq",sep=".")
  filename2 <- paste(n1,"cultivar.tassel",sep=".")
  #colnames(data[[i]]) <- c("CHR","SNP","BP","P")
  colnames(data[[i]]) <- c("SNP","CHR","BP","P")
  qq_dat <- data.frame(obs=-log10(sort(data$P,decreasing =FALSE)),
                       exp=-log10( ppoints(length(data$P))))
  q[[i]]q <- ggplot(data=qq_dat,aes(exp,obs))+
    geom_point(alpha=0.7,color="#7F7F7FFF")+
    geom_abline(color="#D62728FF")+
    xlab("Expected -log10(P-value)")+
    ylab("Observed -log10(P-value)")+
    scale_x_continuous(limits = c(0,7))+
    scale_y_continuous(limits = c(0,7))+
    #ggtitle(name)+
    theme(
      plot.title = element_text(color="red", size=20, face="bold.italic"),
      axis.title = element_text(size=12,face="bold"),
      axis.text = element_text(face="bold",size=8,color = "black"),
      #axis.line = element_line(size=0.8,color="black"),
      axis.ticks= element_line(size=0.8,colour = "black"),
      panel.grid =element_blank(),
      panel.border = element_rect(fill=NA,size = 0.8),
      panel.background = element_blank())
  pdf(paste(filename1,".pdf",sep=""),height = 3,width = 3)
  print(q[[i]])
  dev.off()
  #manhuttan plot
  gwasResults <- data[[i]]
  don <- gwasResults %>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) 
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  p[[i]] p <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1) +
    scale_color_manual(values = rep(c("#999999","#999999", "#0072B2", "#0072B2","#E69F00", "#E69F00"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +   
    theme_bw() +
    ylim(0,7)+
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=8),
      axis.text.y=element_text(size=10),
      axis.title.y=element_text(size = 10),
      axis.title.x=element_text(size = 10),
    )+
    scale_y_continuous(limits = c(0,7))+
    geom_hline(yintercept = -log10(0.00001), colour="red",linetype=2, size=1)
  pdf(paste(filename2,".pdf",sep=""),height = 2,width = 9.6)
  
  print(p[[i]])
  dev.off()
}

pdf("cultivar.test.pdf",height = 2,width = 9.6)
print(p)
dev.off()
