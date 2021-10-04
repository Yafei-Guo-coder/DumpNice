#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
setwd("/Users/guoyafei/Documents/个人项目/傅老师/20210311/GWAS_sign_gene/")
#VRN-A1
annotation_col <- read.table("Annotation_col.txt",header=T,stringsAsFactors = T)
rownames(annotation_col) = c(1:349)
labels_col = rep(c(""), 349)
ann_colors = list(
  Growing_Habit = c(Facultative = "#228B22", Spring="#D95F02", Winter="#1E90FF"),
  Region = c(AM = "#228B22",EA = "#D95F02",SCA="#1E90FF",WA="#E7298A", IND="#000000", JAP="#FFD700")
)
file17 <- read.table("HaploType/17_haplo.txt",header=F,stringsAsFactors = F)
file19 <- read.table("HaploType/19_haplo.txt",header=F,stringsAsFactors = F)
file47 <- read.table("HaploType/47_haplo.txt",header=F,stringsAsFactors = F)
file48 <- read.table("HaploType/48_haplo.txt",header=F,stringsAsFactors = F)
file49 <- read.table("HaploType/49_haplo.txt",header=F,stringsAsFactors = F)
file50 <- read.table("HaploType/50_haplo.txt",header=F,stringsAsFactors = F)
file51 <- read.table("HaploType/51_haplo.txt",header=F,stringsAsFactors = F)
file52 <- read.table("HaploType/52_haplo.txt",header=F,stringsAsFactors = F)
colnames(file17) <- c(1:349)
colnames(file19) <- c(1:349)
colnames(file47) <- c(1:349)
colnames(file48) <- c(1:349)
colnames(file49) <- c(1:349)
colnames(file50) <- c(1:349)
colnames(file51) <- c(1:349)
colnames(file52) <- c(1:349)
setwd("/Users/guoyafei/Documents/个人项目/傅老师/20210311/GWAS_sign_gene/")
data <- read.table("Gene_VariantConponent.txt",header=T,stringsAsFactors = F)

data$group_A1 <- as.factor(data$group_A1)

data$group_D1 <- as.factor(data$group_D1)

data$group_3 <- as.factor(data$group_3)

data$group_4 <- as.factor(data$group_4)

data$group_5 <- as.factor(data$group_5)

data$group_6 <- as.factor(data$group_6)

data$group_7 <- as.factor(data$group_7)

data$group_8 <- as.factor(data$group_8)


Cluster_Method<-c( "ward.D", "single", "complete", "average", "mcquitty", "median", "centroid")
pdf("plot.pdf")
pheatmap(file17, clustering_method = Cluster_Method[1], show_colnames=F,show_rownames=T,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-A1")
ggplot(data, aes(x = group_A1, y = weight_A1))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="GA2ox3-A1") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

pheatmap(file19, clustering_method = Cluster_Method[1],show_colnames=F,show_rownames=T,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=5, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-D1")
ggplot(data, aes(x = group_D1, y = weight_A1))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="GA2ox3-D1") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

pheatmap(file47, clustering_method = Cluster_Method[1],show_colnames=F,show_rownames=T,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=4, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="NGR5-1")
ggplot(data, aes(x = group_3, y = weight_A1))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="NGR5-1") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

pheatmap(file48, clustering_method = Cluster_Method[1],show_colnames=F,show_rownames=T,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=4, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="NGR5-2")
ggplot(data, aes(x = group_4, y = weight_A1))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="NGR5-2") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

pheatmap(file49, clustering_method = Cluster_Method[1],show_colnames=F,show_rownames=T,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=5, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="NGR5-3")
ggplot(data, aes(x = group_5, y = weight_A1))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="NGR5-3") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

pheatmap(file50, clustering_method = Cluster_Method[1],show_colnames=F,show_rownames=T,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=5, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="PIF-1")

ggplot(data, aes(x = group_6, y = weight_A1))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="PIF_1") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))
pheatmap(file51, clustering_method = Cluster_Method[1],show_colnames=F,show_rownames=T,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=7, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="PIF-2")
ggplot(data, aes(x = group_7, y = weight_A1))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="PIF_2") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

pheatmap(file52, clustering_method = Cluster_Method[1],show_colnames=F,show_rownames=T,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="PIF-3")
ggplot(data, aes(x = group_8, y = weight_A1))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="PIF_3") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

dev.off()


pdf("SNP_geno_trait.pdf")
dev.off()
