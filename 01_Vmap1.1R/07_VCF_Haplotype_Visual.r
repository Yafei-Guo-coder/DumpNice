#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5")
#annotation_col <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/E6_Landrace_locate_225.txt",header=T,stringsAsFactors = T)
#rownames(annotation_col) = c(1:225)
#A
annotation_col <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/AB_anno.txt",header=T,stringsAsFactors = T,sep="\t")
#rownames(annotation_col) = c(1:464)
rownames(annotation_col) = c(1:361)
#B
#annotation_col <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/VIP_gene/B_anno.txt",header=T,stringsAsFactors = T,sep="\t")
#rownames(annotation_col) = c(1:383)
#D
annotation_col <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/D_anno.txt",header=T,stringsAsFactors = T,sep="\t")
#rownames(annotation_col) = c(1:283)
rownames(annotation_col) = c(1:250)
seq <- annotation_col[,1]
#plot haplotype heatmap
anno <- annotation_col[,9,drop=FALSE]
brewer.pal(9, "YlGnBu")[c(7)]
"#225EA8"
brewer.pal(9, "Blues")[c(2)]
"#DEEBF7"
brewer.pal(9, "YlOrRd")[c(4)]
"#FEB24C"
brewer.pal(9, "YlOrRd")[c(8)]
"#BD0026"
cols <- c("#225EA8","#DEEBF7","#FEB24C","#BD0026")
#cols <- c("white", "#FFFFCC", "#FEB24C", "#E31A1C")
#cols <- c(brewer.pal(8,"BrBG")[c(1,3,4)], brewer.pal(8,"Set2"))
#AB
ann_color = list(
  #Growing_Habit = c(Facultative = "yellow", Spring="orange", Winter="blue"),
  Region_sub2 = c(Wild_emmer = "#8C510A", Domesticated_emmer = "#DFC27D", Freethreshing = "#F6E8C3", EU="#66C2A5", WA= "#FC8D62",North1="#8DA0CB", North2 ="#E78AC3",South="#A6D854", Tibet="#FFD92F",Other ="#B3B3B3"))
#D
ann_color = list(
  #Growing_Habit = c(Facultative = "yellow", Spring="orange", Winter="blue"),
  Region_sub2 = c(Strangulata = "#8C510A", EU="#66C2A5", WA= "#FC8D62",North1="#8DA0CB", North2 ="#E78AC3",South="#A6D854", Tibet="#FFD92F",Other ="#B3B3B3"))

path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/TXT" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})
#data <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/pos.txt",header=F,stringsAsFactors = F)
#seq <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/anno.txt",header=F,stringsAsFactors = F)
#seq1 <- seq[,1]
tit <- c("Ppd-1 (Tibet Vs South)","Ppd-A1 (North2 Vs South)","Rht-B1 (North2 Vs South & EU Vs South)","Sr33 (WA Vs South)","Sr45 (North2 Vs South & WA Vs EU & EU Vs South)","VRN2-2 (EU Vs South & WA Vs South)") 
pdf("A.pdf")
for (i in c(2)){
  all <- as.data.frame(data[[i]])
  
  all <- all[,seq]
  #colnames(all) <- c(1:464)
  colnames(all) <- c(1:361)
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = anno, annotation_colors = ann_color,annotation_names_col = F, main=tit[i])
}
dev.off()

pdf("B.pdf")
for (i in c(3)){
  all <- as.data.frame(data[[i]])
  all <- all[,seq]
  #colnames(all) <- c(1:383)
  colnames(all) <- c(1:361)
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation_col = anno,annotation_colors = ann_color, annotation_names_col = F, main=tit[i])
}
dev.off()
pdf("D.pdf")
for (i in c(1,4,5,6)){
  all <- as.data.frame(data[[i]])
  all <- all[,seq]
  #colnames(all) <- c(1:283)
  colnames(all) <- c(1:250)
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation_col = anno,annotation_colors = ann_color, annotation_names_col = F, main=tit[i])
}
dev.off()
#annotation_colors = ann_colors
#接下来可以使用06_Qmatrix_PieMap.r来画在地图上的单倍型分布。
