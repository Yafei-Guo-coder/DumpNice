#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/VIP_gene/")

#annotation_col <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/E6_Landrace_locate_225.txt",header=T,stringsAsFactors = T)
#rownames(annotation_col) = c(1:225)
#A
annotation_col <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/VIP_gene/A_anno.txt",header=T,stringsAsFactors = T,sep="\t")
rownames(annotation_col) = c(1:464)
#B
annotation_col <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/VIP_gene/B_anno.txt",header=T,stringsAsFactors = T,sep="\t")
rownames(annotation_col) = c(1:383)
#D
annotation_col <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/VIP_gene/D_anno.txt",header=T,stringsAsFactors = T,sep="\t")
rownames(annotation_col) = c(1:283)

seq <- annotation_col[,1]

#plot haplotype heatmap
anno <- annotation_col[,9,drop=FALSE]
cols <- brewer.pal(8, "PuOr")[c(6,5,8,1)]
cols <- c("white", "#FFFFCC", "#FEB24C", "#E31A1C")
#ann_colors = list(
#  Growing_Habit = c(Facultative = "yellow", Spring="orange", Winter="blue"),
#  Region = c(AF = "#8DD3C7", AM = "#FFFFB3", CA="#BEBADA",EU="#FB8072", NE_A= "#FFED6F",NW_A="#80B1D3", Other ="#FDB462",SA="#B3DE69", SE_A="#FCCDE5",SW_A ="#D9D9D9",Tibet ="#BC80BD",WA="#CCEBC5")
#)
path <- "/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/VIP_gene/TXT/Pos" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})
#data <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/pos.txt",header=F,stringsAsFactors = F)
#seq <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/anno.txt",header=F,stringsAsFactors = F)
#seq1 <- seq[,1]
pdf("A.pdf")
for (i in c(1,2,5,6,8)){
  all <- as.data.frame(data[[i]])
  all <- all[,seq]
  colnames(all) <- c(1:464)
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation_col = anno, cutree_cols = 8, main=i)
}
dev.off()
pdf("B.pdf")
for (i in c(7,9)){
  all <- as.data.frame(data[[i]])
  all <- all[,seq]
  colnames(all) <- c(1:383)
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation_col = anno, cutree_cols = 8, main=i)
}
dev.off()
pdf("D.pdf")
for (i in c(3,4)){
  all <- as.data.frame(data[[i]])
  all <- all[,seq]
  colnames(all) <- c(1:283)
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation_col = anno, cutree_cols = 8, main=i)
}
dev.off()
#annotation_colors = ann_colors
#接下来可以使用06_Qmatrix_PieMap.r来画在地图上的单倍型分布。
