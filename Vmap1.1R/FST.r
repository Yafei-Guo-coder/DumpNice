#画热图 
#文件目录：yafei@203:/data2/yafei/Project3/FST_group/VmapData/heatmap
#A_mean.fst	B_mean.fst	D_mean.fst	AB_mean.fst
#A_weight.fst	B_weight.fst	D_weight.fst	AB_weight.fst

library("corrplot")

#weighted FST

#A lineage
data <- read.table("A_weight.fst",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$Alineage
data<- data[,-1]
data<- as.matrix(data)
pdf("weight_A_fst.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()

#B lineage
data <- read.table("B_weight.fst",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$Blineage
data<- data[,-1]
data<- as.matrix(data)
pdf("weight_B_fst.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()

#D lineage
data <- read.table("D_weight.fst",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$Dlineage
data<- data[,-1]
data<- as.matrix(data)
pdf("weight_D_fst.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()

#AB lineage
data <- read.table("AB_weight.fst",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$ABlineage
data<- data[,-1]
data<- as.matrix(data)
pdf("weight_AB_fst.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()

#mean FST

#A lineage
data <- read.table("A_mean.fst",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$Alineage
data<- data[,-1]
data<- as.matrix(data)
pdf("mean_A_fst.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()

#B lineage
data <- read.table("B_mean.fst",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$Blineage
data<- data[,-1]
data<- as.matrix(data)
pdf("mean_B_fst.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()

#D lineage
data <- read.table("D_mean.fst",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$Dlineage
data<- data[,-1]
data<- as.matrix(data)
pdf("mean_D_fst.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()

#AB lineage
data <- read.table("AB_mean.fst",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$ABlineage
data<- data[,-1]
data<- as.matrix(data)
pdf("mean_AB_fst.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()
