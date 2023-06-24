#plot
setwd("/Users/guoyafei/Desktop/coXiao/单倍型")
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)

data <- read.table("gene.pos-geno.txt",header=T,check.names=F,stringsAsFactors = F)

annotation_col <- read.table("haploinfo.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- annotation_col$ID
anno <- annotation_col[which(rownames(annotation_col) %in% colnames(data)),c(3,4,5),drop=FALSE]

anno2 <- anno[,3,drop=F]
anno2$type <- factor(anno$`Common name`,levels = c("Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Polish wheat","Durum","Rivet wheat","Khorasan wheat","Persian wheat","Spelt","Indian dwarf wheat","Tibetan semi wild","Macha","Club wheat","Vavilovii","Yunan wheat","Xinjiang wheat","Landrace","Cultivar"))
anno2$Ctnt <- factor(anno$Region,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM"))

#anno4 <- anno3[order(anno3$type,anno3$Ctnt),]

### extract data & plot ------>
data2 <- data[,which(colnames(data) %in% rownames(anno))]
data3 <- data2[,rownames(anno)]
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
pdf(output[i],width=12,height = 8 )
pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, cluster_col = F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"),cluster_row = FALSE, annotation_col = anno2, annotation_names_col = F)
dev.off()

#haploinfo2
annotation_col <- read.table("haploinfo2.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- annotation_col$ID
anno <- annotation_col[which(rownames(annotation_col) %in% colnames(data)),c(4,5,6),drop=FALSE]

anno2 <- anno[,3,drop=F]
anno2$type <- factor(anno$type,levels = c("Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","free-threshing tetraploids","Spelt","Indian dwarf wheat","Tibetan semi wild","Macha","Club wheat","Vavilovii","Yunan wheat","Xinjiang wheat","Landrace","Cultivar"))
anno2$Ctnt <- factor(anno$Region,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM"))

#anno4 <- anno3[order(anno3$type,anno3$Ctnt),]

### extract data & plot ------>
data2 <- data[,which(colnames(data) %in% rownames(anno))]
data3 <- data2[,rownames(anno)]
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
pdf(output[i],width=12,height = 8 )
pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, cluster_col = F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"),cluster_row = FALSE, annotation_col = anno2, annotation_names_col = F)
dev.off()

#渗入分析
library(ggplot2)
setwd("/Users/guoyafei/Desktop/coXiao/渗入")
data <- read.table("chr2A.withAnc-1.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-All.csv",header=T,stringsAsFactors = F, sep=",")

data <- read.table("chr2A.withAnc-D-1.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-D-All.csv",header=T,stringsAsFactors = F, sep=",")

data <- read.table("chr2A.withAnc-U-1.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-U-All.csv",header=T,stringsAsFactors = F, sep=",")

data <- data[which(data$fd >= 0 & data$fd < 1 & data$D > 0),]
#data[which(data$fd < 0 | data$fd > 1 | data$D <= 0),10] <- 0

ggplot(data, aes(x=start, y=fd)) +
  #geom_point(size=0.5)+
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("fd") +
  #geom_hline(yintercept=0.5,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759452125,color='red')
