library(ggplot2)
library(dplyr)
library(ggridges)
library(RColorBrewer)
library(ggrepel)
require(gridExtra)
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
#提取266clonegene的位置
#path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/fst_266out"
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/fst_266out"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene_266 <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})
#提取已克隆的基因的位置----
names <- read.table("nameMap.txt",header=F,stringsAsFactors = F)
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density/fst_266out"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
gene <- lapply(filePath, function(x){
  read.table(x, header=T, stringsAsFactors = F, sep="\t")})
#分页画图----

p <- list()
for(i in c(2)){
  #rownames(gene[[i]]) <- gene[[i]]$NAME
  data[[i]]$Pop <- "Overall"
  sub <- data[[i]][, c(1,2,3,6,2,3,4,7)]
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
    geom_density(alpha = 0.2) +
    theme_classic()+
    theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
    xlab("Fst") + ylab("Proportion") +
    ggtitle(names[i,2]) + xlim(0,0.6) +
    geom_vline(xintercept = c(d,c,a,b), color = c("#F8766D","#7CAE00","#008FC4","#C77CFF"), size= 0.9, linetype = "dotted") + 
    #geom_point(data = gene[[i]], aes(MEAN_FST, 0), color = 'red') +
    #geom_text_repel(data = gene[[i]], aes(MEAN_FST, 0, label=rownames(gene[[i]]))) +
    #geom_label_repel(data = gene[[i]],aes(MEAN_FST, Y), label=rownames(gene[[i]]),segment.colour = NA,colour="white", segment.colour="black") +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"), legend.position = "none",legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y =element_blank() )
}

pdf("test1.pdf",height = 10,width = 12)
grid.arrange(p[[1]],p[[7]],p[[13]],p[[2]],p[[8]],p[[14]],p[[3]],p[[9]],p[[15]],nrow=3)
dev.off()
pdf("test2.pdf",height = 10,width = 12)
grid.arrange(p[[4]],p[[10]],p[[16]],p[[5]],p[[11]],p[[17]],p[[6]],p[[12]],p[[18]],nrow=3)
dev.off()

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

for(i in c(1:18)){
  print(gene_resis[[i]][which(gene_resis[[i]]$MEAN_FST > 0.2),])
  print(gene_flow[[i]][which(gene_flow[[i]]$MEAN_FST > 0.2),])
  print(gene_deve[[i]][which(gene_deve[[i]]$MEAN_FST > 0.2),])
  print("ok")
}

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
p <- list()
for(i in c(1,7,13,8,14,3,9,15)){
  a <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Abiotic_stimulus"),4])
  b <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Disease_resistance"),4])
  c <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Growing"),4])
  d <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Processing_quality"),4])
  e <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Rhythm"),4])
  point <- gene_266[[i]][which(gene_266[[i]]$Y != 0),]
  data[[i]]$Anno <- "Z"
  f <- median(data[[i]][which(data[[i]]$MEAN_FST >0),6])
  all <- rbind(data[[i]][,c(1,2,3,6,7)],gene_266[[i]][,c(1,2,3,4,8)])
  all <- all[which(all$MEAN_FST >0),]
  p[[i]] <- ggplot(all, aes(MEAN_FST, fill= Anno, color=Anno)) +
    geom_density(alpha = 0.2) +
    scale_color_manual(values = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4")) +
    scale_fill_manual(values = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4")) +
    #scale_colour_discrete(breaks = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6"), labels = c("Disease_resistance","Rhythm","Growing","Flowering","Gibberellin_related"))+
    #粉；樱桃红；青绿色；橙黄色；青色
    theme_classic() +
    xlab("Fst") + ylab("Proportion") +
    ggtitle(names[i,2]) + 
    xlim(0,0.6) +
    geom_vline(xintercept = c(a,b,c,d,e,f), color = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4"), size= 0.9, linetype = "dotted") + 
    geom_point(data = point, aes(MEAN_FST, 0), color = 'red') +
    geom_label_repel(data = point,aes(MEAN_FST, Y), label=point$Name,segment.colour = NA,colour="white", segment.colour="black") +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position = "none", legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y =element_blank())
}

for(i in c(11,4,10,5,17,6,12,18)){  
  a <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Abiotic_stimulus"),4])
  b <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Disease_resistance"),4])
  c <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Growing"),4])
  d <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Processing_quality"),4])
  e <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Rhythm"),4])
  point <- gene_266[[i]][which(gene_266[[i]]$Y != 0),]
  
  data[[i]]$Anno <- "Z"
  f <- median(data[[i]][which(data[[i]]$MEAN_FST >0),6])
  all <- rbind(data[[i]][,c(1,2,3,6,7)],gene_266[[i]][,c(1,2,3,4,8)])
  all <- all[which(all$MEAN_FST >0),]
  
  p[[i]] <- ggplot(all, aes(MEAN_FST, fill= Anno, color=Anno)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_color_manual(values = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4")) +
    scale_fill_manual(values = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4")) +
    xlab("Fst") + ylab("Proportion") +
    ggtitle(names[i,2]) + 
    xlim(0,0.8) +
    geom_vline(xintercept = c(a,b,c,d,e,f), color = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4"), size= 0.9, linetype = "dotted") + 
    geom_point(data = point, aes(MEAN_FST, 0), color = 'red') +
    geom_label_repel(data = point,aes(MEAN_FST, Y), label=point$Name,segment.colour = NA,colour="white", segment.colour="black") +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position = "none", legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y =element_blank())
}

for(i in c(2,8)){
  a <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST > 0 & gene_266[[i]]$Anno == "Abiotic_stimulus"),4])
  b <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST > 0 & gene_266[[i]]$Anno == "Disease_resistance"),4])
  c <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST > 0 & gene_266[[i]]$Anno == "Growing"),4])
  d <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST > 0 & gene_266[[i]]$Anno == "Processing_quality"),4])
  e <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST > 0 & gene_266[[i]]$Anno == "Rhythm"),4])
  point <- gene_266[[i]][which(gene_266[[i]]$Y != 0),]
  
  data[[i]]$Anno <- "Z"
  f <- median(data[[i]][which(data[[i]]$MEAN_FST >0),6])
  all <- rbind(data[[i]][,c(1,2,3,6,7)],gene_266[[i]][,c(1,2,3,4,8)])
  all <- all[which(all$MEAN_FST >0),]
  
  p[[i]] <- ggplot(all, aes(MEAN_FST, fill= Anno, color=Anno)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_color_manual(values = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4")) +
    scale_fill_manual(values = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4")) +
    xlab("Fst") + ylab("Proportion") +
    ggtitle(names[i,2]) + 
    xlim(0,0.6) +
    geom_vline(xintercept = c(a,b,c,d,e,f), color = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4"), size= 0.9, linetype = "dotted") + 
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position = "none", legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y =element_blank())
}

for(i in c(14,18)){  
  #a <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Abiotic_stimulus"),4])
  b <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Disease_resistance"),4])
  c <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Growing"),4])
  #d <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Processing_quality"),4])
  e <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Rhythm"),4])
  point <- gene_266[[i]][which(gene_266[[i]]$Y != 0),]
  
  data[[i]]$Anno <- "Z"
  f <- median(data[[i]][which(data[[i]]$MEAN_FST >0),6])
  all <- rbind(data[[i]][,c(1,2,3,6,7)],gene_266[[i]][,c(1,2,3,4,8)])
  all <- all[which(all$MEAN_FST >0 & all$Anno != "Abiotic_stimulus" & all$Anno != "Processing_quality"),]
  
  p[[i]] <- ggplot(all, aes(MEAN_FST, fill= Anno, color=Anno)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_color_manual(values = c("#F8766D","#879F00","#00B0F6","#780ec4")) +
    scale_fill_manual(values = c("#F8766D","#879F00","#00B0F6","#780ec4")) +
    xlab("Fst") + ylab("Proportion") +
    ggtitle(names[i,2]) + 
    xlim(0,1) +
    geom_vline(xintercept = c(b,c,e,f), color = c("#F8766D","#879F00","#00B0F6","#780ec4"), size= 0.9, linetype = "dotted") + 
    geom_point(data = point, aes(MEAN_FST, 0), color = 'red') +
    geom_label_repel(data = point,aes(MEAN_FST, Y), label=point$Name,segment.colour = NA,colour="white", segment.colour="black") +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position = "none", legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y =element_blank())
}
 
for(i in c(16)){ 
  #a <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Abiotic_stimulus"),4])
  b <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Disease_resistance"),4])
  c <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Growing"),4])
  #d <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Processing_quality"),4])
  e <- median(gene_266[[i]][which(gene_266[[i]]$MEAN_FST >0 & gene_266[[i]]$Anno == "Rhythm"),4])
  point <- gene_266[[i]][which(gene_266[[i]]$Y != 0),]
  
  data[[i]]$Anno <- "Z"
  f <- median(data[[i]][which(data[[i]]$MEAN_FST >0),6])
  all <- rbind(data[[i]][,c(1,2,3,6,7)],gene_266[[i]][,c(1,2,3,4,8)])
  all <- all[which(all$MEAN_FST >0 & all$Anno != "Abiotic_stimulus" & all$Anno != "Processing_quality"),]

  p[[i]] <- ggplot(all, aes(MEAN_FST, fill= Anno, color=Anno)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_color_manual(values = c("#F8766D","#879F00","#00B0F6","#780ec4")) +
    scale_fill_manual(values = c("#F8766D","#879F00","#00B0F6","#780ec4")) +
    xlab("Fst") + ylab("Proportion") +
    ggtitle(names[i,2]) + 
    xlim(0,0.8) +
    geom_vline(xintercept = c(b,c,e,f), color = c("#F8766D","#879F00","#00B0F6","#780ec4"), size= 0.9, linetype = "dotted") + 
    #geom_point(data = point, aes(MEAN_FST, 0), color = 'red') +
    #geom_label_repel(data = point,aes(MEAN_FST, Y), label=point$Name,segment.colour = NA,colour="white", segment.colour="black") +
    theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position = "none", legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y =element_blank())
}
pdf("Test3.pdf",height = 10,width = 12)
grid.arrange(p[[1]],p[[7]],p[[13]],p[[2]],p[[8]],p[[14]],p[[3]],p[[9]],p[[15]],nrow=3)
dev.off()
pdf("Test4.pdf",height = 10,width = 12)
grid.arrange(p[[4]],p[[10]],p[[16]],p[[5]],p[[11]],p[[17]],p[[6]],p[[12]],p[[18]],nrow=3)
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

ggplot(all, aes(MEAN_FST, fill=Pop)) +
  geom_density(alpha = 0.2) +
  theme_classic()+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  xlab("Fst") + ylab("Proportion") +
  ggtitle(names[i,2]) + xlim(0,0.6) +
  geom_vline(xintercept = c(d,c,a,b), color = c("#F8766D","#7CAE00","#008FC4","#C77CFF"), size= 0.9, linetype = "dotted") + 
  geom_point(data = gene[[1]], aes(MEAN_FST, 0), color = 'red') +
  geom_text_repel(data = gene[[1]],aes(MEAN_FST, 0, label = gene[[1]]$Gene_id)) +
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"), legend.position = "none",legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

for(i in c(1:18)){
  print(gene_resis[[i]][which(gene_resis[[i]]$MEAN_FST > 0.2),])
  print(gene_flow[[i]][which(gene_flow[[i]]$MEAN_FST > 0.2),])
  print(gene_deve[[i]][which(gene_deve[[i]]$MEAN_FST > 0.2),])
  print("ok")
}

a <- median(gene_266[[8]][which(gene_266[[8]]$MEAN_FST >0 & gene_266[[8]]$Anno == "Abiotic_stimulus"),4])
b <- median(gene_266[[8]][which(gene_266[[8]]$MEAN_FST >0 & gene_266[[8]]$Anno == "Disease_resistance"),4])
c <- median(gene_266[[8]][which(gene_266[[8]]$MEAN_FST >0 & gene_266[[8]]$Anno == "Growing"),4])
d <- median(gene_266[[8]][which(gene_266[[8]]$MEAN_FST >0 & gene_266[[8]]$Anno == "Processing_quality"),4])
e <- median(gene_266[[8]][which(gene_266[[8]]$MEAN_FST >0 & gene_266[[8]]$Anno == "Rhythm"),4])
point <- gene_266[[8]][which(gene_266[[8]]$Y != 0),]
data[[8]]$Anno <- "Z"
f <- median(data[[8]][which(data[[8]]$MEAN_FST >0),6])
all <- rbind(data[[8]][,c(1,2,3,6,7)],gene_266[[8]][,c(1,2,3,4,8)])
all <- all[which(all$MEAN_FST >0),]
p[[8]] <- ggplot(all, aes(MEAN_FST, fill= Anno, color=Anno)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_color_manual(values = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4")) +
  scale_fill_manual(values = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4")) +
  xlab("Fst") + ylab("Proportion") +
  ggtitle(names[8,2]) + 
  xlim(0,0.6) +
  geom_vline(xintercept = c(a,b,c,d,e,f), color = c("#E76BF3","#F8766D","#879F00","#f8c36d","#00B0F6","#780ec4"), size= 0.9, linetype = "dotted") + 
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.position = "none", legend.text = element_text(size = 10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y =element_blank())
