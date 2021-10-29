library(ggplot2)
library(dplyr)
library(ggridges)
library(RColorBrewer)
library(ggrepel)
require(gridExtra)
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
#names <- read.table("nameMap.txt",header=F,stringsAsFactors = F)
#提取已克隆的基因的位置----
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/fst_clone"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})
#画图:密度图&标注关注基因----
pdf("fst_select_gene.pdf",height = 5,width = 10)
for(i in c(1:18)){
  gene[[i]]$Pop <- data[[i]][1,7]
  p <- ggplot(data[[i]], aes(MEAN_FST, colour = Pop)) +
    scale_fill_manual() +
    stat_density(alpha = 0.2) +
    theme_classic()+
    #theme(axis.title.y = element_blank()) +
    xlab("Fst") + ylab("Proportion") + xlim(0,0.6) +
    #geom_vline(xintercept = 1, color = 'gray', size = 0.5) + 
    geom_point(data = gene[[i]], aes(MEAN_FST, 0), color = 'red') +
    geom_text_repel(data = gene[[i]],aes(MEAN_FST, 0, label = gene[[i]]$Name)) +
    ggtitle(names[i,2]) +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position="none",legend.text = element_blank(),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
  print(p)
}
dev.off()
#提取abioticgene的位置----
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
#画图:叠加密度图----
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
#提取resistgene的位置----
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/resistgene"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene_resis <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})
#提取flowegene的位置
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/flowergene"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene_flow <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})
#提取developgene的位置
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/developgene"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene_deve <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})
#画图:叠加密度图----
pdf("test.pdf",height = 5,width = 10)
for(i in c(1:18)){
  data[[i]]$Pop <- "Overall"
  gene_resis[[i]]$Pop <- "resistgene"
  gene_flow[[i]]$Pop <- "flowegene"
  gene_deve[[i]]$Pop <- "developgene"
  all <- rbind(data[[i]],gene_resis[[i]],gene_flow[[i]],gene_deve[[i]])
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
#画图:密度图&标注关注基因位置-----
pdf("fst_develop_gene.pdf",height = 5,width = 10)
for(i in c(1:18)){
  gene_deve[[i]]$Pop <- data[[i]][1,7]
  p <- ggplot(data[[i]], aes(MEAN_FST, colour = Pop)) +
    scale_fill_manual() +
    stat_density(alpha = 0.2) +
    theme_classic()+
    #theme(axis.title.y = element_blank()) +
    xlab("Fst") + ylab("Proportion") + xlim(0,0.6) +
    #geom_vline(xintercept = 1, color = 'gray', size = 0.5) + 
    geom_point(data = gene_deve[[i]], aes(MEAN_FST, 0), color = 'red') +
    geom_text_repel(data = gene_deve[[i]],aes(MEAN_FST, 0, label = gene_deve[[i]]$Gene_id)) +
    ggtitle(names[i,2]) +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position="none",legend.text = element_blank(),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
  print(p)
}
dev.off()

#分页画图-----
names <- read.table("nameMap.txt",header=F,stringsAsFactors = F)
#all <- data.frame(CHROM="", BIN_START="", BIN_END="", MEAN_FST="", Gene_start="", Gene_end="",Gene_id="",Pop="", site="", lineage="", stringsAsFactors=FALSE)
p <- list()
for(i in c(1:18)){
  data[[i]]$Pop <- "Overall"
  sub <- data[[i]][,c(1,2,3,6,2,3,4,7)]
  colnames(sub) <- c("CHROM","BIN_START","BIN_END","MEAN_FST","Gene_start","Gene_end","Gene_id","Pop")
  gene_resis[[i]]$Pop <- "resistgene"
  gene_flow[[i]]$Pop <- "flowegene"
  gene_deve[[i]]$Pop <- "developgene"
  d <- median(gene_deve[[i]][which(gene_deve[[i]]$MEAN_FST >0),4])
  c <- median(gene_flow[[i]][which(gene_flow[[i]]$MEAN_FST >0),4])
  a <- median(sub[which(sub$MEAN_FST >0),4])
  b <- median(gene_resis[[i]][which(gene_resis[[i]]$MEAN_FST >0),4])
  all <- rbind(sub,gene_resis[[i]],gene_flow[[i]],gene_deve[[i]])
  p[[i]] <- ggplot(all, aes(MEAN_FST, fill=Pop)) +
    #facet_grid(lineage ~ site)+
    geom_density(alpha = 0.2) +
    theme_classic()+
    #theme(axis.title.y = element_blank()) +
    xlab("Fst") + ylab("Proportion") +
    ggtitle(names[i,2]) + xlim(0,0.6) +
    geom_vline(xintercept = c(d,c,a,b), color = c("#F8766D","#7CAE00","#008FC4","#C77CFF"), size= 0.9, linetype = "dotted") + 
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"), legend.position = "none",legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  #print(p)
}
pdf("Pop1_sub_fst.pdf",height = 10,width = 12)
grid.arrange(p[[1]],p[[7]],p[[13]],p[[2]],p[[8]],p[[14]],p[[3]],p[[9]],p[[15]],nrow=3)
dev.off()
pdf("Pop2_sub_fst.pdf",height = 10,width = 12)
grid.arrange(p[[4]],p[[10]],p[[16]],p[[5]],p[[11]],p[[17]],p[[6]],p[[12]],p[[18]],nrow=3)
dev.off()
