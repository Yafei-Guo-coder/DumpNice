#FastCall2
setwd("/Users/guoyafei/Desktop/FC2/compare/2.0")

library(ggplot2)
input <- paste("chr",1:42,".all.txt",sep="")
output <- paste("chr",1:42,".all.pdf",sep="")
for(i in c(1:42)){
  data <- read.table(input[i], header=F, stringsAsFactors = F)
  pdf(output[i])
  p <- ggplot(data,aes(V2,fill=V1))+
    geom_histogram(bins=10,aes(x=V2, y=..density..),position="identity",alpha = 0.5)+
    theme_bw()
  plot(p)
  dev.off()
}


data <- read.table("/Users/guoyafei/Desktop/FC2/compare/2.0_raw/file2.0_raw.all.txt", header=F,stringsAsFactors = F)
ggplot(data,aes(V2,fill=V1))+
  geom_histogram(bins=10,aes(x=V2, y=..density..),position = 'dodge',alpha = 0.9)+
  scale_fill_manual(values = c("#FDC086","#BEAED4")) +
  theme_bw()

data <- read.table("/Users/guoyafei/Desktop/FC2/compare/2.1/file2.1.all.txt", header=F,stringsAsFactors = F)
ggplot(data,aes(V2,fill=V1))+
  geom_histogram(bins=10,aes(x=V2, y=..density..),position = 'dodge',alpha = 0.9)+
  scale_fill_manual(values = c("#FDC086","#BEAED4")) +
  theme_bw()

data <- read.table("/Users/guoyafei/Desktop/FC2/compare/filegatk.all.txt", header=F,stringsAsFactors = F)
ggplot(data,aes(V2,fill=V1))+
  geom_histogram(bins=15,aes(x=V2, y=..density..),position="identity",alpha = 0.5)+
  theme_bw()