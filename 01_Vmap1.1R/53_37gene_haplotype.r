#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
cols <- c("#225EA8","#DEEBF7","#FEB24C","#BD0026")
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/VCF_withAn/225_TXT")
#annotation_col <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/E6_Landrace_locate_225.txt",header=T,stringsAsFactors = T)
#rownames(annotation_col) = c(1:225)
annotation <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/select_taxa4.txt", header=T, row.names=1,  stringsAsFactors = F, sep="\t")
anno <- annotation[,7,drop=FALSE]

#--------------------------------------------------------------------------------------------------------------------------------------
#AB
AB_name <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/VCF_withAn/225_TXT/AB_name.txt",header=F,stringsAsFactors = F)
AB_anno <- anno[which(rownames(anno) %in% AB_name[,1]),1,drop=F]
ann_color = list(
  heatmap = c(Wild_emmer = "#8C510A", Domesticated_emmer = "#DFC27D", Freethreshing = "#F6E8C3", EU="#66C2A5", WA= "#FC8D62",CA="#8DA0CB", EA ="#E78AC3",SA="#A6D854",Tibet="#FFD92F"))
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/VCF_withAn/225_TXT/AB" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})
AB_file <-read.table("AB_file_name.txt",header=F,stringsAsFactors = F)
AB_tit <- AB_file[,1]
pdf("AB.pdf")
for (i in c(1:length(data))){
  all <- as.data.frame(data[[i]])
  colnames(all) <- AB_name[,1]
  all <- all[, rownames(AB_anno)]
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE,color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = AB_anno, annotation_colors = ann_color,annotation_names_col = F,main=AB_tit[i])
}
dev.off()

#--------------------------------------------------------------------------------------------------------------------------------------
#D
ann_color = list(
  #Growing_Habit = c(Facultative = "yellow", Spring="orange", Winter="blue"),
  heatmap = c(Strangulata = "#8C510A", EU="#66C2A5", WA= "#FC8D62",CA="#8DA0CB", EA ="#E78AC3",SA="#A6D854",Tibet="#FFD92F"))

path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/VCF_withAn/225_TXT/D" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})

D_name <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/VCF_withAn/225_TXT/D_name.txt",header=F,stringsAsFactors = F)
D_anno <- anno[which(rownames(anno) %in% D_name[,1]),1,drop=F]
D_file <-read.table("D_file_name.txt",header=F,stringsAsFactors = F)
D_tit <- D_file[,1]
pdf("D.pdf")
for (i in c(1:length(data))){
  all <- as.data.frame(data[[i]])
  colnames(all) <- D_name[,1]
  all <- all[, rownames(D_anno)]
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE,color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = D_anno, annotation_colors = ann_color,annotation_names_col = F,main=D_tit[i])
}
dev.off()

#annotation_colors = ann_colors
#接下来可以使用06_Qmatrix_PieMap.r来画在地图上的单倍型分布。

#plot haplotype heatmap
brewer.pal(9, "YlGnBu")[c(7)]
"#225EA8"
brewer.pal(9, "Blues")[c(2)]
"#DEEBF7"
brewer.pal(9, "YlOrRd")[c(4)]
"#FEB24C"
brewer.pal(9, "YlOrRd")[c(8)]
"#BD0026"
cols <- c("#225EA8","#DEEBF7","#FEB24C","#BD0026")
