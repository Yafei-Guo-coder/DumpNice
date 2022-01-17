#taxa QC
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
display.brewer.all()
brewer.pal(8, "Set1")
p <- list()
for (i in c(1:21)){
  name <- paste("chr",i,sep="")
  filename <- paste(i,"_taxa.txt",sep="")
  data <- read.table(filename,header=T,stringsAsFactors = F)
  data$HeterozygousProportion <- as.numeric(data$HeterozygousProportion)
  data$MissRate <- as.numeric(data$MissRate)
  p[[i]] <- ggplot(data=data,aes(x = HeterozygousProportion))+
    #geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
    geom_density(color = "black", fill = "#FB8072", alpha=0.8)+
    xlab("HeterozygousProportion")+
    ylab("")+
    scale_x_continuous(limits = c(0,0.5))+
    #scale_y_continuous(limits = c(0,7))+
    ggtitle(name)+
    theme(
      plot.title = element_text(color="red", size=20, face="bold.italic"),
      axis.title = element_text(size=12,face="bold"),
      axis.text = element_text(face="bold",size=8,color = "black"),
      axis.ticks= element_line(size=0.8,colour = "black"),
      panel.grid =element_blank(),
      panel.border = element_rect(fill=NA,size = 0.8),
      panel.background = element_blank())
}
pdf("A_taxa_het.pdf",width = 9,height = 9)
grid.arrange(p[[1]],p[[4]],p[[7]],p[[10]],p[[13]],p[[16]],p[[19]],nrow=3)
dev.off()
pdf("B_taxa_het.pdf",width = 9,height = 9)
grid.arrange(p[[2]],p[[5]],p[[8]],p[[11]],p[[14]],p[[17]],p[[20]],nrow=3)
dev.off()
pdf("D_taxa_het.pdf",width = 9,height = 9)
grid.arrange(p[[3]],p[[6]],p[[9]],p[[12]],p[[15]],p[[18]],p[[21]],nrow=3)
dev.off()

p <- list()
for (i in c(1:21)){
  name <- paste("chr",i,sep="")
  filename <- paste(i,"_taxa.txt",sep="")
  data <- read.table(filename,header=T,stringsAsFactors = F)
  data$HeterozygousProportion <- as.numeric(data$HeterozygousProportion)
  data$MissRate <- as.numeric(data$MissRate)
  p[[i]] <- ggplot(data=data,aes(x = MissRate))+
    geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
    xlab("MissRate")+
    ylab("")+
    scale_x_continuous(limits = c(0,0.5))+
    #scale_y_continuous(limits = c(0,7))+
    ggtitle(name)+
    theme(
      plot.title = element_text(color="red", size=20, face="bold.italic"),
      axis.title = element_text(size=12,face="bold"),
      axis.text = element_text(face="bold",size=8,color = "black"),
      axis.ticks= element_line(size=0.8,colour = "black"),
      panel.grid =element_blank(),
      panel.border = element_rect(fill=NA,size = 0.8),
      panel.background = element_blank())
}
pdf("A_taxa_MissRate.pdf",width = 9,height = 9)
grid.arrange(p[[1]],p[[4]],p[[7]],p[[10]],p[[13]],p[[16]],p[[19]],nrow=3)
dev.off()
pdf("B_taxa_MissRate.pdf",width = 9,height = 9)
grid.arrange(p[[2]],p[[5]],p[[8]],p[[11]],p[[14]],p[[17]],p[[20]],nrow=3)
dev.off()
pdf("D_taxa_MissRate.pdf",width = 9,height = 9)
grid.arrange(p[[3]],p[[6]],p[[9]],p[[12]],p[[15]],p[[18]],p[[21]],nrow=3)
dev.off()

library(ggplot2)
library(gridExtra)
pdf("all_taxa.pdf",width = 9,height = 9)
data <- read.table("all_taxa.txt",header=T,stringsAsFactors = F)
data$HeterozygousProportion <- as.numeric(data$HeterozygousProportion)
data$MissRate <- as.numeric(data$MissRate)
ggplot(data=data,aes(x = HeterozygousProportion))+
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  xlab("HeterozygousProportion")+
  ylab("")+
  scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,7))+
  theme(
    plot.title = element_text(color="red", size=20, face="bold.italic"),
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())
ggplot(data=data,aes(x = MissRate))+
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  xlab("MissRate")+
  ylab("")+
  
  scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,7))+
  theme(
    plot.title = element_text(color="red", size=20, face="bold.italic"),
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())
dev.off()