library(ggplot2)
library(dplyr)
library(ggridges)
library(RColorBrewer)
library(ggrepel)
#display.brewer.all()
col1 = c(brewer.pal(9,'Pastel1'),brewer.pal(10,'Set3'))
col1 <- c(rep("#FB8072", times=7),rep("#80B1D3", times=6),rep("#FDB462", times=6))

#读取FST文件
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density")
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/fst"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})
names <- read.table("nameMap.txt",header=F,stringsAsFactors = F)
all <- data[[1]]
all$Pop <- NA
colnames(all)[1:6] = c("Chr","Pos1","Pos2","Site","WEIGHTED_FST","MEAN_FST")
for(i in c(1:length(data))){
  colnames(data[[i]]) = c("Chr","Pos1","Pos2","Site","WEIGHTED_FST","MEAN_FST")
  data[[i]]$Pop <- names[i,1]
  all <- rbind(all,data[[i]])
}
all = all[which(all$MEAN_FST >= 0),]
all = all[which(all$Pop!="NA"),]
unique(all$Pop)
all$Pop = factor(all$Pop, levels=rev(unique(all$Pop)))

#提取已克隆的基因的位置
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/fst_clone"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})

#画图
pdf("fst_select_gene.pdf",height = 5,width = 10)
for(i in c(1:18)){
  gene[[i]]$Pop <- data[[i]][1,7]
  p <- ggplot(data[[i]], aes(MEAN_FST, colour = Pop)) +
    scale_fill_manual() +
    stat_density(alpha = 0.2) +
    theme_classic()+
    #theme(axis.title.y = element_blank()) +
    xlab("Fst") + ylab("Proportion") + xlim(0,0.5) +
    #geom_vline(xintercept = 1, color = 'gray', size = 0.5) + 
    geom_point(data = gene[[i]], aes(MEAN_FST, 0), color = 'red') +
    geom_text_repel(data = gene[[i]],aes(MEAN_FST, 0, label = gene[[i]]$Name))+
    ggtitle(names[i,2])+
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position="none",legend.text = element_blank(),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
  print(p)
}
dev.off()

#提取abioticgene的位置
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/abioticgene"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene_abio <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})

#提取bioticgene的位置
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/bioticgene"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene_bio <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})

#提取backgroundgene的位置
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/backgroudgene"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene_back <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})

#画图
pdf("all.pdf",height = 5,width = 10)
for(i in c(1:18)){
  gene_abio[[i]]$Pop <- "abioticgene"
  gene_bio[[i]]$Pop <- "bioticgene"
  gene_back[[i]]$Pop <- "backgroundgene"
  all <- rbind(gene_abio[[i]],gene_bio[[i]],gene_back[[i]])
  p <- ggplot(all, aes(MEAN_FST, fill= Pop)) +
    geom_density(alpha = 0.2) +
    theme_classic()+
    #theme(axis.title.y = element_blank()) +
    xlab("Fst") + ylab("Proportion") +
    ggtitle(names[i,2]) + xlim(0,0.5) +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position="none",legend.text = element_blank(),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
  print(p)
}
dev.off()


