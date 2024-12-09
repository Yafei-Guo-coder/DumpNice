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

setwd("/Users/guoyafei/Desktop/JM44/xpclr/g6")
data <- read.table("cultivar_landrace.Alineage.smooth.21chr.txt", header=F,stringsAsFactors = F)
data <- read.table("cultivar_landrace.Blineage.smooth.21chr.txt", header=F,stringsAsFactors = F)
data <- read.table("cultivar_landrace.Dlineage.smooth.21chr.txt", header=F,stringsAsFactors = F)

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

pos_1A_L <- c(4095727,4202215,6517557,6561934)
pos_1A_H <- c(508723999,508924985)
pos_1A_G <- c(3094761,3776501,3795341,3845783,3846832,3910395,3910525,4033473,4041458,15924605)
pos_1B_L <- c(5519697,5584166,5686611,6436429,22185903)
pos_1B_H <- c(555765127,555933489)
pos_1B_G <- c(4965759,4991949,5003864,5056936,5099315,5137194,5142157,21843089)
pos_1D_L <- c(84815,99355,3709624,3825716,4341746,5139253,5356898,7155533)
pos_1D_H <- c(412160786,412217776)
pos_1D_G <- c(250463,259657,327118,339837,381137)
pos_6A_G <- c(24921897,25070522,25107550,25203705,25314279,25408473,25472841,25526000,25581178)
pos_6B_G <- c(43384370,43436446,43560471,43571670,43820081,43858073,43868418,43922841,43932231,44061847,62653760,62657777)

setwd("/Users/guoyafei/Desktop/JM44/xpclr/g6")
data_A <- read.table("cultivar_landrace.Alineage.smooth.21chr.txt", header=F,stringsAsFactors = F)
quantile(data$V4,0.95)
#3.594556
data_B <- read.table("cultivar_landrace.Blineage.smooth.21chr.txt", header=F,stringsAsFactors = F)
quantile(data$V4,0.95)
#2.908107
data_D <- read.table("cultivar_landrace.Dlineage.smooth.21chr.txt", header=F,stringsAsFactors = F)
quantile(data$V4,0.95)

inputA <- paste("cultivar_",rep(c("1950s_1960s.","1960s_1970s.","1970s_1980s.","1980s_1990s.","1990s_2000s.","2000s_2010s."), each=1),"Alineage.","smooth.21chr.txt",sep="")
threshA <- c(3.53, 3.66, 3.11, 2.37, 3.36, 2.22)
outputA <- paste("cultivar_",rep(c("1950s_1960s.","1960s_1970s.","1970s_1980s.","1980s_1990s.","1990s_2000s.","2000s_2010s."), each=1),"Alineage.","smooth.21chr",sep="")
inputB <- paste("cultivar_",rep(c("1950s_1960s.","1960s_1970s.","1970s_1980s.","1980s_1990s.","1990s_2000s.","2000s_2010s."), each=1),"Blineage.","smooth.21chr.txt",sep="")
threshB <- c(3.96, 3.92, 3.02, 2.65, 2.91, 2.05)
outputB <- paste("cultivar_",rep(c("1950s_1960s.","1960s_1970s.","1970s_1980s.","1980s_1990s.","1990s_2000s.","2000s_2010s."), each=1),"Blineage.","smooth.21chr",sep="")
inputD <- paste("cultivar_",rep(c("1950s_1960s.","1960s_1970s.","1970s_1980s.","1980s_1990s.","1990s_2000s.","2000s_2010s."), each=1),"Dlineage.","smooth.21chr.txt",sep="")
threshD <- c(3.44, 3.47, 3.36, 2.56, 2.70, 2.58)
outputD <- paste("cultivar_",rep(c("1950s_1960s.","1960s_1970s.","1970s_1980s.","1980s_1990s.","1990s_2000s.","2000s_2010s."), each=1),"Dlineage.","smooth.21chr",sep="")

#period 1A----
output <- paste(outputA,".1A.pdf",sep="")
for(i in c(1:length(inputA))){
  data <- read.table(inputA[i], header=F,stringsAsFactors = F)
  sub <- data[which(data$V1 == "1"),]
  pdf(output[i], width=11.5,height=2.5)
  p <- ggplot(sub) +
    geom_vline(xintercept=pos_1A_L,color='red',size=0.1,alpha=0.5)+ 
    geom_vline(xintercept=pos_1A_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
    geom_vline(xintercept=pos_1A_G,color='#A07D35',size=0.1,alpha=0.5)+ #棕
    geom_line(aes(x=V2, y=V6),size=0.3) +
    geom_hline(yintercept=threshA[i], color='#E69F00',linetype = "dashed") +
    theme_classic() +
    #theme(axis.title.y = element_blank()) 
    xlab("Position") + ylab("xpclr score") 
  print(p)
  dev.off()
}

#period 1B----
output <- paste(outputB,".1B.pdf",sep="")
for(i in c(1:length(inputB))){
  data <- read.table(inputB[i], header=F,stringsAsFactors = F)
  sub <- data[which(data$V1 == "2"),]
  pdf(output[i], width=11.5,height=2.5)
  p <- ggplot(sub, aes(x=V2, y=V6)) +
    geom_line(size=0.5) +
    theme_classic() +
    #theme(axis.title.y = element_blank()) 
    xlab("Position") + ylab("xpclr score") +
    geom_hline(yintercept=threshB[i], color='#E69F00',linetype = "dashed") +
    geom_vline(xintercept=pos_1B_L,color='red',size=0.1,alpha=0.5)+ 
    geom_vline(xintercept=pos_1B_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
    geom_vline(xintercept=pos_1B_G,color='#A07D35',size=0.1,alpha=0.5) 
  print(p)
  dev.off()
}

#period 1D----
output <- paste(outputD,".1D.pdf",sep="")
for(i in c(1:length(inputD))){
  data <- read.table(inputD[i], header=F,stringsAsFactors = F)
  sub <- data[which(data$V1 == "3"),]
  
  pdf(output[i], width=11.5,height=2.5)
  p <- ggplot(sub, aes(x=V2, y=V6)) +
    geom_line(size=0.5) +
    theme_classic() +
    #theme(axis.title.y = element_blank()) 
    xlab("Position") + ylab("xpclr score") +
    geom_hline(yintercept=threshD[i], color='#E69F00',linetype = "dashed") +
    geom_vline(xintercept=pos_1D_L,color='red',size=0.1,alpha=0.5)+ 
    geom_vline(xintercept=pos_1D_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
    geom_vline(xintercept=pos_1D_G,color='#A07D35',size=0.1,alpha=0.5)
  print(p) 
  dev.off()
}

#period 6A----
output <- paste(outputA,".6A.pdf",sep="")
for(i in c(1:length(inputA))){
  data <- read.table(inputA[i], header=F,stringsAsFactors = F)
  sub <- data[which(data$V1 == "16"),]
  
  pdf(output[i], width=11.5,height=2.5)
  
  p <- ggplot(sub, aes(x=V2, y=V6)) +
    geom_line(size=0.5) +
    theme_classic() +
    #theme(axis.title.y = element_blank()) 
    xlab("Position") + ylab("xpclr score") +
    geom_hline(yintercept=threshA[i], color='#E69F00',linetype = "dashed") +
    geom_vline(xintercept=pos_6A_G,color='#A07D35',size=0.1,alpha=0.5)
  print(p)
  dev.off()
}

#period 6B----
output <- paste(outputB,".6B.pdf",sep="")
for(i in c(1:length(inputB))){
  data <- read.table(inputB[i], header=F,stringsAsFactors = F)
  sub <- data[which(data$V1 == "17"),]
  
  pdf(output[i], width=11.5,height=2.5)
  
  p <- ggplot(sub, aes(x=V2, y=V6)) +
    geom_line(size=0.5) +
    theme_classic() +
    #theme(axis.title.y = element_blank()) 
    xlab("Position") + ylab("xpclr score") +
    geom_hline(yintercept=threshB[i], color='#E69F00',linetype = "dashed") +
    geom_vline(xintercept=pos_6B_G,color='#A07D35',size=0.1,alpha=0.5) 
  print(p)
  dev.off()
}










#3.364177
#landrace_cultivar 1A----
sub <- data_A[which(data_A$V1 == "1"),]
pdf("landrace-cultivar-1A.pdf", width=11.5,height=2.5)
ggplot(sub) +
  geom_vline(xintercept=pos_1A_L,color='red',size=0.1,alpha=0.5)+ 
  geom_vline(xintercept=pos_1A_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
  geom_vline(xintercept=pos_1A_G,color='#A07D35',size=0.1,alpha=0.5) +#棕
  geom_line(aes(x=V2, y=V4),size=0.3) +
  geom_hline(yintercept=3.594556, color='#E69F00',linetype = "dashed") +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") 
dev.off()
#landrace_cultivar 1B----
sub <- data_B[which(data_B$V1 == "2"),]
pdf("landrace-cultivar-1B.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=2.908107, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_1B_L,color='red',size=0.1,alpha=0.5)+ 
  geom_vline(xintercept=pos_1B_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
  geom_vline(xintercept=pos_1B_G,color='#A07D35',size=0.1,alpha=0.5) 
dev.off()
#landrace_cultivar 1D----
sub <- data_D[which(data_D$V1 == "3"),]
pdf("landrace-cultivar-1D.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.364177, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_1D_L,color='red',size=0.1,alpha=0.5)+ 
  geom_vline(xintercept=pos_1D_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
  geom_vline(xintercept=pos_1D_G,color='#A07D35',size=0.1,alpha=0.5) 
dev.off()
#landrace_cultivar 6A----
sub <- data_A[which(data_A$V1 == "16"),]
pdf("landrace-cultivar-6A.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.594556, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_6A_G,color='#A07D35',size=0.1,alpha=0.5)
dev.off()
#landrace_cultivar 6B----
sub <- data_B[which(data_B$V1 == "17"),]
pdf("landrace-cultivar-6B.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=2.908107, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_6B_G,color='#A07D35',size=0.1,alpha=0.5) 
dev.off()

data <- read.table("/Users/guoyafei/Desktop/JM44/xpclr/g6/wildeinkorn_domeinkorn_final.txt", header=F,stringsAsFactors = F)
quantile(data$V4,0.95)
#3.38848 



data <- read.table("/Users/guoyafei/Desktop/JM44/xpclr/g6/wildeinkorn_domeinkorn_final.txt", header=F,stringsAsFactors = F)
quantile(data$V4,0.95)
#3.38848 
#WEI-DEI 1A----
sub <- data[which(data$V1 == "chr1A"),]
pdf("WEI-DEI-1A.pdf", width=11.5,height=2.5)
ggplot(sub) +
  geom_vline(xintercept=pos_1A_L,color='red',size=0.1,alpha=0.5)+ 
  geom_vline(xintercept=pos_1A_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
  geom_vline(xintercept=pos_1A_G,color='#A07D35',size=0.1,alpha=0.5) +#棕
  geom_line(aes(x=V2, y=V4),size=0.3) +
  geom_hline(yintercept=3.38848, color='#E69F00',linetype = "dashed") +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") 
dev.off()

#WEI-DEI 6A----
sub <- data[which(data$V1 == "chr6A"),]
pdf("WEI-DEI-6A.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.38848, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_6A_G,color='#A07D35',size=0.1,alpha=0.5)
dev.off()




data <- read.table("/Users/guoyafei/Desktop/JM44/xpclr/g6/wildemmer_domemmer_final.txt", header=F,stringsAsFactors = F)
data_A <- data[which(data$V1 %in% c("chr1A","chr2A","chr3A","chr4A","chr5A","chr6A","chr7A")),]
data_B <- data[which(data$V1 %in% c("chr1B","chr2B","chr3B","chr4B","chr5B","chr6B","chr7B")),]
quantile(data_A$V4,0.95)
#2.005881 
quantile(data_B$V4,0.95)
#2.464232
#WE-DE 1A----
sub <- data_A[which(data_A$V1 == "chr1A"),]
pdf("WE_DE-1A.pdf", width=11.5,height=2.5)
ggplot(sub) +
  geom_vline(xintercept=pos_1A_L,color='red',size=0.1,alpha=0.5)+ 
  geom_vline(xintercept=pos_1A_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
  geom_vline(xintercept=pos_1A_G,color='#A07D35',size=0.1,alpha=0.5) +#棕
  geom_line(aes(x=V2, y=V4),size=0.3) +
  geom_hline(yintercept=2.005881, color='#E69F00',linetype = "dashed") +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") 
dev.off()
#WE-DE 1B----
sub <- data_B[which(data_B$V1 == "chr1B"),]
pdf("WE_DE-1B.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=2.464232, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_1B_L,color='red',size=0.1,alpha=0.5)+ 
  geom_vline(xintercept=pos_1B_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
  geom_vline(xintercept=pos_1B_G,color='#A07D35',size=0.1,alpha=0.5) 
dev.off()

#WE-DE 6A----
sub <- data_A[which(data_A$V1 == "chr6A"),]
pdf("WE_DE-6A.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=2.005881, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_6A_G,color='#A07D35',size=0.1,alpha=0.5)
dev.off()
#WE-DE 6B----
sub <- data_B[which(data_B$V1 == "chr6B"),]
pdf("WE_DE-6B.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=2.464232, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_6B_G,color='#A07D35',size=0.1,alpha=0.5) 
dev.off()

data <- read.table("/Users/guoyafei/Desktop/JM44/xpclr/g6/domemmer_durum_final.txt", header=F,stringsAsFactors = F)
data_A <- data[which(data$V1 %in% c("chr1A","chr2A","chr3A","chr4A","chr5A","chr6A","chr7A")),]
data_B <- data[which(data$V1 %in% c("chr1B","chr2B","chr3B","chr4B","chr5B","chr6B","chr7B")),]
quantile(data_A$V4,0.95)
#3.093504
quantile(data_B$V4,0.95)
#4.044948
#DE-FT 1A----
sub <- data_A[which(data_A$V1 == "chr1A"),]
pdf("DE_FT-1A.pdf", width=11.5,height=2.5)
ggplot(sub) +
  geom_vline(xintercept=pos_1A_L,color='red',size=0.1,alpha=0.5)+ 
  geom_vline(xintercept=pos_1A_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
  geom_vline(xintercept=pos_1A_G,color='#A07D35',size=0.1,alpha=0.5) +#棕
  geom_line(aes(x=V2, y=V4),size=0.3) +
  geom_hline(yintercept=3.093504, color='#E69F00',linetype = "dashed") +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") 
dev.off()
#DE-FT 1B----
sub <- data_B[which(data_B$V1 == "chr1B"),]
pdf("DE_FT-1B.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=4.044948, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_1B_L,color='red',size=0.1,alpha=0.5)+ 
  geom_vline(xintercept=pos_1B_H,color='#56B4E9',size=0.1,alpha=0.5)+ #蓝
  geom_vline(xintercept=pos_1B_G,color='#A07D35',size=0.1,alpha=0.5) 
dev.off()

#DE-FT 6A----
sub <- data_A[which(data_A$V1 == "chr6A"),]
pdf("DE_FT-6A.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.093504, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_6A_G,color='#A07D35',size=0.1,alpha=0.5)
dev.off()
#DE-FT 6B----
sub <- data_B[which(data_B$V1 == "chr6B"),]
pdf("DE_FT-6B.pdf", width=11.5,height=2.5)
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5) +
  theme_classic() +
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=4.044948, color='#E69F00',linetype = "dashed") +
  geom_vline(xintercept=pos_6B_G,color='#A07D35',size=0.1,alpha=0.5) 
dev.off()


#alpha_cs:1A----
sub <- data[which(data$V1 == "chr1A"),]
sub <- data[which(data$V1 == "17"),]
ggplot(sub, aes(x=V2, y=V4)) +
  geom_line(size=0.5)+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("xpclr score") +
  geom_hline(yintercept=3.36, color='#E69F00',linetype = "dashed") 
  #geom_vline(xintercept=3094761,color='red')+
  #geom_vline(xintercept=4965759,color='red')+
  #geom_vline(xintercept=5137194,color='red')+
  #geom_vline(xintercept=5142157,color='red')

  #geom_vline(xintercept=3910525,color='red')+
  #geom_vline(xintercept=4095727,color='red')+
  #geom_vline(xintercept=213000000,color='red')

#alpha_cs:6A----
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

#alpha_cs:6B----
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

#alpha_cs:1B----
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

#gamma_cs:1D----
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
####分类画图#####
setwd("/Users/guoyafei/Desktop/JM44/plot")
library(RColorBrewer)
data <-read.table("plot.txt", header=F,stringsAsFactors = F)
data <-read.table("plot2.txt", header=F,stringsAsFactors = F)
data <- read.table("plot_PicMin.txt", header=T, stringsAsFactors = F)
data <- read.table("plot_periodPicMin.txt", header=T, stringsAsFactors = F)
data <- read.table("~/Desktop/plot.txt",header=F,stringsAsFactors = F)
pdf("test1.pdf", height=5, width=3)
ggplot(data, aes(y = p, group=as.factor(data$type), fill=as.factor(data$type)))+
  geom_boxplot(notch = F) +
  theme_classic() + 
  #scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
  scale_fill_manual(values = c(brewer.pal(11, "Set3")[c(4,5)],"#DEEBF7","#FEB24C"))+
  ylab("p")+
  theme(legend.title="")+
  #ggtitle(colnames(out)[m])+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
dev.off()
pdf("test.pdf", height=14, width=8)
ggplot(data, aes(y = V4, group=as.factor(data$V5), fill=as.factor(data$V5)))+
  geom_boxplot(notch = F) + 
  facet_grid(V1~.)+
  theme_classic() + 
  #scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
  scale_fill_manual(values = c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C"))+
  ylab("XPCLR")+
  theme(legend.title="")+
  #ggtitle(colnames(out)[m])+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
dev.off()


#### phyperTest #####
setwd("/Users/guoyafei/Desktop/network/2018science/phyperTest")
shuf_integra <- read.table("shuf/integra_abiotic.txt", header=T, stringsAsFactors = F)
shuf_phyper <- read.table("shuf/phyper_abiotic.txt", header=T, stringsAsFactors = F)




