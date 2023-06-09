# 相关系数
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
setwd("/Users/guoyafei/Desktop/btr1/ibs")
input <- c("btr1-A.10k.all.ibs","btr1-A.1M.all.ibs","btr1-A.200k.all.ibs","btr1-A.2k.all.ibs","btr1-A.500k.all.ibs","btr1-B.10k.all.ibs","btr1-B.1M.all.ibs","btr1-B.200k.all.ibs","btr1-B.2k.all.ibs","btr1-B.500k.all.ibs")
output <- c("btr1-A.10k.all.ibs.pdf","btr1-A.1M.all.ibs.pdf","btr1-A.200k.all.ibs.pdf","btr1-A.2k.all.ibs.pdf","btr1-A.500k.all.ibs.pdf","btr1-B.10k.all.ibs.pdf","btr1-B.1M.all.ibs.pdf","btr1-B.200k.all.ibs.pdf","btr1-B.2k.all.ibs.pdf","btr1-B.500k.all.ibs.pdf")

for(i in c(1:length(input))){
  test <- read.table(input[i],row.names = 1, header = T)
  pdf(output[i],width=10,height=10)
  pheatmap(test, show_rownames=T, show_colnames=T, cluster_col = T, clustering_method="ward.D2",cluster_row = T,main="ward.D2")
  pheatmap(test, show_rownames=T, show_colnames=T, cluster_col = T, clustering_method="ward.D2",cluster_row = T,main="complete")
  pheatmap(test, show_rownames=T, show_colnames=T, cluster_col = T, clustering_method="ward.D2",cluster_row = T,main="average")
  pheatmap(test, show_rownames=T, show_colnames=T, cluster_col = T, clustering_method="ward.D2",cluster_row = T,main="mcquitty")
  pheatmap(test, show_rownames=T, show_colnames=T, cluster_col = T, clustering_method="ward.D2",cluster_row = T,main="median")
  pheatmap(test, show_rownames=T, show_colnames=T, cluster_col = T, clustering_method="ward.D2",cluster_row = T,main="centroid")
  pheatmap(test, show_rownames=T, show_colnames=T, cluster_col = T, clustering_method="ward.D2",cluster_row = T,main="single")
  dev.off()
}
