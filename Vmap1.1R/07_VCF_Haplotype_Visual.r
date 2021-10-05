#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/gene5k")
annotation_col <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/E6_Landrace_locate_225.txt",header=T,stringsAsFactors = T)
rownames(annotation_col) = c(1:225)
#labels_col = rep(c(""), 225)
anno <- annotation_col[,8,drop=FALSE]

colors <- brewer.pal(n = 10, name = "Set3")[c(1:10)]

seq <- annotation_col[,1]

ann_colors = list(
  Growing_Habit = c(Facultative = "yellow", Spring="orange", Winter="blue"),
  Region = c(AF = "#66C2A5", AM = "#FC8D62", CA="#8DA0CB",EAE="#FFD92F", EAM="#E78AC3", EAN ="#A6D854",EAS="#FFD92F", EU="#E78AC3",SA ="#A6D854",WA ="#A6D854")
)

path <- "/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/haplotype5k" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})

data <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/pos.txt",header=F,stringsAsFactors = F)
seq <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Gene/anno.txt",header=F,stringsAsFactors = F)
seq1 <- seq[,1]
pdf("3.pdf")
for (i in c(3)){
  all <- as.data.frame(data[[i]])
  all <- all[,seq]
  colnames(all) <- c(1:225)
  pheatmap(all, show_rownames=FALSE, show_colnames=FALSE, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = FALSE, cluster_row = FALSE, annotation_col = anno, annotation_colors = ann_colors, main=i)
}
dev.off()

all <-data[,seq1]
anno <- seq[,3,drop=FALSE]
colnames(all) <- c(1:58)
pheatmap(all, show_rownames=FALSE, show_colnames=FALSE, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = FALSE, cluster_row = FALSE, annotation_col = anno, main="Sr33")




