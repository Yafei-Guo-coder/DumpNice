#画图文件：yafei@203
#/data2/yafei/Project3/CS_Vmap/A_ibs.txt
#/data2/yafei/Project3/CS_Vmap/B_ibs.txt
#/data2/yafei/Project3/CS_Vmap/D_ibs.txt

#工作目录：yafei@203:/data2/yafei/Project3/CS_Vmap

#热图  
library("corrplot")
#A lineage
ibs <- read.table("A_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

#A:0.5 B:0.5 D:0.7
pdf("IBS_A_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.5))
dev.off()

#B lineage
ibs <- read.table("B_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

pdf("IBS_B_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.5))
dev.off()

#D lineage
ibs <- read.table("D_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

pdf("IBS_A_CS_heat.pdf",width=10,height=10)
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=0.6,number.digits=0.3,tl.cex=1,cl.cex=1.2,cl.lim = c(0, 0.7))
dev.off()
