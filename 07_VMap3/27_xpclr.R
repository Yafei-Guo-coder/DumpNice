#xpclr
setwd("/Users/guoyafei/Documents/02_VmapIII/22_xpclr")

################################################ points intersect #############################################
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)

setwd("/Users/guoyafei/Documents/02_VmapIII/22_xpclr")
shuf5k <- read.table("shuf5k/AB.year.7.shuf5k", header=T, stringsAsFactors = F)
rownames(shuf5k) <- shuf5k$gene

path <- "/Users/guoyafei/Documents/02_VmapIII/22_xpclr/module" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data1 <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep="\t")})

all1 <- read.table("allgene/AB.year.7", header=T, stringsAsFactors = F)
rownames(all1) <- all1$gene

#AB-year
thresh <- c(3.32304533697249, 2.85365870779203, 2.86502700188772, 2.35072437470687, 2.82586032822358, 3.06320636816039, 2.67427291961055)
#AB-region
#thresh <- c(2.89502846023601, 2.16224947799078, 2.47694881306687, 2.01117096911165, 2.28196828145439)
#D-year
#thresh <- c(2.91572659731368, 2.8903418869464, 2.89755603608904, 2.51492298054788, 3.17230587582881, 3.04380068663059, 2.63912099621883)
#D-region
#thresh <- c(2.81978999736161, 2.42461298752664, 2.2807268299181, 2.31906477839178, 2.3057143939575)

for (i in c(1:80)) {
  name <- as.data.frame(data1[[i]])
  rownames(name) <- name$V1
  sub <- all1[which(all1$gene %in% name$V1),]
  out <- paste(names(data1)[i],"pdf",sep=".")
  p <- list()
  p[[1]] <- ggplot(shuf5k, aes(X1950, X1960)) +
    geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
    geom_point(data = sub, aes(X1950, X1960, size=1, alpha = 0.5,colour = "red")) +
    geom_vline(xintercept=thresh[1],color='red',linetype = "dotted")+
    geom_hline(yintercept=thresh[2],color='red',linetype = "dotted")+
    #geom_text_repel(data = sub,aes(V2, V3, label = name),max.overlaps = 100) +
    xlab("1950-1960")+
    ylab("1960-1970")+
    theme_bw()+
    theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  p[[2]] <- ggplot(shuf5k, aes(X1960, X1970)) +
    geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
    geom_point(data = sub, aes(X1960, X1970, size=2, alpha = 0.5,colour = "red")) +
    geom_vline(xintercept=thresh[2],color='red',linetype = "dotted")+
    geom_hline(yintercept=thresh[3],color='red',linetype = "dotted")+
    #geom_text_repel(data = sub,aes(V2, V3, label = name),max.overlaps = 100) +
    xlab("1960-1970")+
    ylab("1970-1980")+
    theme_bw()+
    theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  p[[3]] <- ggplot(shuf5k, aes(X1970, X1980)) +
    geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
    geom_point(data = sub, aes(X1970, X1980, size=2, alpha = 0.5,colour = "red")) +
    geom_vline(xintercept=thresh[3],color='red',linetype = "dotted")+
    geom_hline(yintercept=thresh[4],color='red',linetype = "dotted")+
    #geom_text_repel(data = sub,aes(V2, V3, label = name),max.overlaps = 100) +
    xlab("1970-1980")+
    ylab("1980-1990")+
    theme_bw()+
    theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  p[[4]] <- ggplot(shuf5k, aes(X1980, X1990)) +
    geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
    geom_point(data = sub, aes(X1980, X1990, size=2, alpha = 0.5,colour = "red")) +
    geom_vline(xintercept=thresh[4],color='red',linetype = "dotted")+
    geom_hline(yintercept=thresh[5],color='red',linetype = "dotted")+
    #geom_text_repel(data = sub,aes(V2, V3, label = name),max.overlaps = 100) +
    xlab("1980-1990")+
    ylab("1990-2000")+
    theme_bw()+
    theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  p[[5]] <- ggplot(shuf5k, aes(X1990, X2000)) +
    geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
    geom_point(data = sub, aes(X1990, X2000, size=2, alpha = 0.5,colour = "red")) +
    geom_vline(xintercept=thresh[5],color='red',linetype = "dotted")+
    geom_hline(yintercept=thresh[6],color='red',linetype = "dotted")+
    #geom_text_repel(data = sub,aes(V2, V3, label = name),max.overlaps = 100) +
    xlab("1990-2000")+
    ylab("2000-2010")+
    theme_bw()+
    theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  p[[6]] <- ggplot(shuf5k, aes(X2000, X2010)) +
    geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
    geom_point(data = sub, aes(X2000, X2010, size=2, alpha = 0.5,colour = "red")) +
    geom_vline(xintercept=thresh[6],color='red',linetype = "dotted")+
    geom_hline(yintercept=thresh[7],color='red',linetype = "dotted")+
    #geom_text_repel(data = sub,aes(V2, V3, label = name),max.overlaps = 100) +
    xlab("2000-2010")+
    ylab("2010-2020")+
    theme_bw()+
    theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  pdf(out,width = 12,height = 8)
  grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],nrow=2)
  dev.off()
}




