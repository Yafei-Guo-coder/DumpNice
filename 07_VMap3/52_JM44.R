#JM44########
setwd("~/Desktop/JM44/")
data <- read.table("chr1D.21chr.bed",header=F,stringsAsFactors = F)

sub <- data[which(data$start > 740000000),]

ggplot(data, aes(x=V2, y=V6)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("XP-CLR score") +
  geom_hline(yintercept=2.52 ,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=412217776 ,color='red')

############################################################ 单倍型分析 ################################################
#plot
setwd("/Users/guoyafei/Desktop/JM44")
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)

#data <- read.table("TraesCS1D02G317211-geno.txt", header=T, check.names=F, stringsAsFactors = F)
#A
data <- read.table("TraesCS1A02G317311-geno.txt", header=T, check.names=F, stringsAsFactors = F)
data_A <- data[c(3,14),]
#data <- read.table("TraesCS1A02G466500LC-geno.txt", header=T, check.names=F, stringsAsFactors = F)
#B
data <- read.table("TraesCS1B02G329711-geno.txt", header=T, check.names=F, stringsAsFactors = F)
data_B1 <- data[c(1,3,16),]
data <- read.table("TraesCS1B02G329992-geno.txt", header=T, check.names=F, stringsAsFactors = F)
data_B2 <- data[45,,drop=F]
#D
data <- read.table("TraesCS1D02G317211-geno.txt", header=T, check.names=F, stringsAsFactors = F)
data_D1 <- data[c(6,10),]
data <- read.table("TraesCS1D02G317301-geno.txt", header=T, check.names=F, stringsAsFactors = F)
data_D2 <- data[c(2,5,10,13,15),]


#zuoA
data <- read.table("TraesCS1B02G329992-geno.txt", header=T, check.names=F, stringsAsFactors = F)
#zuoB
data <- read.table("TraesCS1B02G320000-geno.txt", header=T, check.names=F, stringsAsFactors = F)
#zuoD
data <- read.table("TraesCS1D02G308600-geno.txt", header=T, check.names=F, stringsAsFactors = F)

data <- rbind(data_A,data_B1,data_B2)
data <- rbind(data_D1,data_D2)
annotation_col <- read.table("/Users/guoyafei/Desktop/JM44/abd_186anno.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")
annotation_col <- read.table("/Users/guoyafei/Desktop/JM44/hap_class.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")

rownames(annotation_col) <- annotation_col$VcfID
anno <- annotation_col[which(rownames(annotation_col) %in% colnames(data)),c(2,9),drop=FALSE]
anno$Subspecies <- factor(anno$Subspecies,levels = c("Landrace","Cultivar"))
data2 <- data[,which(colnames(data) %in% rownames(anno))]
data3 <- data2[,rownames(anno)]
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
#pdf("/Users/guoyafei/Desktop/coTong/test_tong_region2_all.pdf",width=12,height = 8 )
pheatmap(data3, show_rownames=FALSE, show_colnames=T, cluster_col = F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"),cluster_row = FALSE, annotation_col = anno, annotation_names_col = F)
#dev.off()

#zuoA
data <- read.table("zuoA-geno.txt", header=T, check.names=F, stringsAsFactors = F)
#zuoB
data <- read.table("zuoB-geno.txt", header=T, check.names=F, stringsAsFactors = F)
#zuoD
data <- read.table("zuoD-geno.txt", header=T, check.names=F, stringsAsFactors = F)

#annotation_col <- read.table("/Users/guoyafei/Desktop/coTong/haploinfo_V3.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")
annotation_col <- read.table("ABD_annotation.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")
#annotation_col <- read.table("ABD_only_annotation.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")

rownames(annotation_col) <- annotation_col$ID

subA <- data_A[,which(colnames(data_A) %in% rownames(annotation_col))]
subB1 <- data_B1[,which(colnames(data_B1) %in% rownames(annotation_col))]
subB1[subB1 == 2] <- -1
subB1[subB1 == 0] <- 2
subB1[subB1 == -1] <- 0
subB2 <- data_B2[,which(colnames(data_B2) %in% rownames(annotation_col))]
subB2[subB2 == 2] <- -1
subB2[subB2 == 0] <- 2
subB2[subB2 == -1] <- 0
subD1 <- data_D1[,which(colnames(data_D1) %in% rownames(annotation_col))]
subD2 <- data_D2[,which(colnames(data_D2) %in% rownames(annotation_col))]
data <- rbind(subA,subB1,subB2,subD1,subD2)

#anno <- annotation_col[which(rownames(annotation_col) %in% colnames(data)),c(3,4,5),drop=FALSE]
anno <- annotation_col[which(rownames(annotation_col) %in% colnames(data)),c(2,3,4,5,6,9,13),drop=FALSE]
#anno2 <- anno[,3,drop=F]
anno2 <- anno[order(anno$seq_year),c(5,6,7)]
#anno2$type_final <- factor(anno2$type_final,levels = c("Wild_emmer","Domesticated_emmer","Free_threshing_tetraploids","Landrace","Cultivar"))
anno2$type_final <- factor(anno2$type_final,levels = c("Landrace","Cultivar","Strangulata"))

#anno2$type <- factor(anno$type_final,levels = c("Wild_emmer","Domesticated_emmer","Free_threshing_tetraploids","Landrace","Cultivar"))
#anno2$Ctnt <- factor(anno$`Ctnt(9)`,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM","-"))
#anno2$Ctnt <- factor(anno$Region,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM"))
#anno4 <- anno3[order(anno3$type,anno3$Ctnt),]
### extract data & plot ------>
data2 <- data[,which(colnames(data) %in% rownames(anno2))]
data3 <- data2[,rownames(anno2)]
#cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
#pdf("/Users/guoyafei/Desktop/coTong/test_tong_region2_all.pdf",width=12,height = 8 )
pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, cluster_col = T, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"),cluster_row = FALSE, annotation_col = anno2, annotation_names_col = F)
#dev.off()

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
pdf("test3.pdf",width=60,height = 40)
pheatmap(data3, show_rownames=FALSE, show_colnames=T, cluster_col = F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"),cluster_row = FALSE, annotation_col = anno2, annotation_names_col = F)
dev.off()
############################################################ xpclr选择分析 #############################################
#折线图
library(ggplot2)
setwd("/Users/guoyafei/Desktop/JM44/选择")
data <- read.table("landrace_cultivar_final.txt",header=F,stringsAsFactors = F, sep="\t")
quantile(data$V4,0.95)

setwd("/Users/guoyafei/Desktop/JM44/xpclr/v1")
data <- read.table("landrace_cultivar_final.txt",header=F,stringsAsFactors = F, sep="\t")
quantile(data$V4,0.95)
#3.23

setwd("/Users/guoyafei/Desktop/JM44/xpclr/v3")
data <- read.table("cultivar_landrace.21chr.Alineage.smooth.txt",header=F,stringsAsFactors = F, sep="\t")
quantile(data$V6,0.95)
#3.52
data <- read.table("cultivar_landrace.21chr.Blineage.smooth.txt",header=F,stringsAsFactors = F, sep="\t")
quantile(data$V6,0.95)
#3.52
data <- read.table("cultivar_landrace.21chr.Dlineage.smooth.txt",header=F,stringsAsFactors = F, sep="\t")
quantile(data$V6,0.95)
#2.65

#alpha_cs:1A
sub <- data[which(data$V1 == "chr1A"),]
sub <- data[which(data$V1 == "1"),]
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.52, color='#E69F00',linetype = "dashed") 
  #geom_vline(xintercept=3910525,color='red')+
  #geom_vline(xintercept=4095727,color='red')+
  #geom_vline(xintercept=213000000,color='red')
  
#alpha_cs:6A
sub <- data[which(data$V1 == "chr6A"),]
sub <- data[which(data$V1 == "16"),]
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.52, color='#E69F00',linetype = "dashed")
  #geom_vline(xintercept=24921897,color='red')+
  #geom_vline(xintercept=283000000,color='red')


#alpha_cs:6B
sub <- data[which(data$V1 == "chr6B"),]
sub <- data[which(data$V1 == "17"),]
ggplot(sub, aes(x=V2, y=V6)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.23, color='#E69F00',linetype = "dashed")
#geom_vline(xintercept=24921897,color='red')+
#geom_vline(xintercept=283000000,color='red')

#alpha_cs:1B
sub <- data[which(data$V1 == "chr1B"),]
sub <- data[which(data$V1 == "2"),]
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.52, color='#E69F00',linetype = "dashed")
  #geom_vline(xintercept=5003864,color='red')+
  #geom_vline(xintercept=5056936,color='red')+
  #geom_vline(xintercept=5686611,color='blue')+
  #geom_vline(xintercept=238000000,color='blue')

sub <- data[which(data$V1 == "chr1B" & data$V2 < 10000000),]
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.23, color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=5003864,color='red')+
  geom_vline(xintercept=5056936,color='red')+
  geom_vline(xintercept=5686611,color='blue')

#gamma_cs:1D
sub <- data[which(data$V1 == "chr1D"),]
sub <- data[which(data$V1 == "3"),]
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=2.65, color='#E69F00',linetype = "dashed")
  #geom_vline(xintercept=381137,color='red')+
  #geom_vline(xintercept=170000000,color='red')


sub <- data[which(data$V1 == "chr1D"  & data$V2 < 10000000),]
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.23, color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=381137,color='red')