**统计每个文件里的未分离位点
setwd("/data1/home/yafei/Project3/Download")
EC <- read.table("new_EC",header=T,stringsAsFactors=F)
TR <- read.table("new_TR",header=T,stringsAsFactors=F)
WW <- read.table("new_WW",header=T,stringsAsFactors=F)

bing <- read.table("/data1/home/yafei/Project3/Vmap1.1/out.txt",header=F,stringsAsFactors=F)
EC_id <- EC[,3]
WW_id <- WW[,3]
TR_id <- TR[,3]

noSep_EC <- read.table("EC_noSep.txt",header=F,stringsAsFactors=F)
noSep_TR <- read.table("TR_noSep.txt",header=F,stringsAsFactors=F)
noSep_WW <- read.table("WW_noSep.txt",header=F,stringsAsFactors=F)

#1. 算各文件的分离位点
EC_sep <- EC_id[-which(EC_id %in% noSep_EC[,1])]
WW_sep <- WW_id[-which(WW_id %in% noSep_WW[,1])]
TR_sep <- TR_id[-which(TR_id %in% noSep_TR[,1])]
#2. 求其他文件的分离位点
TR_WW_V <- c(TR_sep,WW_sep,bing[,1])
EC_WW_V <- c(EC_sep,WW_sep,bing[,1])
EC_TR_V <- c(EC_sep,TR_sep,bing[,1])
#3. 求实际各文件的非分离位点
noSep_real_EC <- noSep_EC[-which(noSep_EC[,1] %in% TR_WW_V),1]
noSep_real_TR <- noSep_TR[-which(noSep_TR[,1] %in% EC_WW_V),1]
noSep_real_WW <- noSep_WW[-which(noSep_WW[,1] %in% EC_TR_V),1]
#4. 求实际各文件的分离位点
Sep_real_EC <- EC_id[-which(EC_id %in% noSep_real_EC)]
Sep_real_TR <- TR_id[-which(TR_id %in% noSep_real_TR)]
Sep_real_WW <- WW_id[-which(WW_id %in% noSep_real_WW)]
#画韦恩图
#三个文件
Length_EC<-length(Sep_real_EC)
Length_TR<-length(Sep_real_TR)
Length_WW<-length(Sep_real_WW)
Length_EC_TR<-length(intersect(Sep_real_EC,Sep_real_TR))
Length_WW_TR<-length(intersect(Sep_real_WW,Sep_real_TR))
Length_EC_WW<-length(intersect(Sep_real_EC,Sep_real_WW))
library(VennDiagram)
T<-venn.diagram(list(EC_86=Sep_real_EC,TR_28=Sep_real_TR,WW_741=Sep_real_WW),filename=NULL,lwd=1,lty=2,col=c('red','green','blue'),fill=c('red','green','blue'),cat.col=c('red','green','blue'),reverse=TRUE)
pdf("Venn.pdf")
grid.draw(T)
dev.off()
#两个文件
pdf("Venn2.pdf")
draw.pairwise.venn(area1=Length_A,area2=Length_B,cross.area=Length_AB,category=c('A','B'),lwd=rep(1,1),lty=rep(2,2),col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),rotation.degree=90)
dev.off()

write.table(Sep_real_EC,"Sep_real_EC", quote=F,row.names=F)
write.table(Sep_real_TR,"Sep_real_TR", quote=F,row.names=F)
write.table(Sep_real_WW,"Sep_real_WW", quote=F,row.names=F)

#输出hapscanner的输入文件：pos.txt & posAllele.txt
Real_WW <- read.table("Sep_real_WW", header=F,stringsAsFactors=F)
selec <- WW[which(WW$ID %in% Real_WW[,1]),] 
posAllele <- selec[,c(1,2,4,5)]
colnames(posAllele) <- c("Chr","Pos","Ref","Alt")
pos <- selec[,c(1,2)]
colnames(pos) <- c("Chr","Pos")

for(i in 1:42){
a <- posAllele[which(posAllele$Chr == i),]
fileA <- paste("/data1/home/yafei/Project3/Download/Hapscanner/chr",i,"_posAllele.txt",sep="")
write.table(a,fileA,row.names=F,quote=F,sep="\t")
b <- pos[which(pos$Chr == i),]
fileB <- paste("/data1/home/yafei/Project3/Download/Hapscanner/chr",i,"_pos.txt",sep="")
write.table(b,fileB,row.names=F,quote=F,sep="\t")
}



