#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require(rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)

setwd("/Users/guoyafei/Desktop/群体结构")
#load Qmatrix files
sfiles <- list.files(path="/Users/guoyafei/Desktop/群体结构/Q矩阵", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))

pdf("Q.pdf",width = 24,height = 8)
p1 <- plotQ(slist[2:6],returnplot=T,exportplot=F,basesize=11,
            sortind="all",showindlab=F,showyaxis=T,showticks=T,sharedindlab=T,clustercol=brewer.pal(7, "Set2"))
grid.arrange(p1$plot[[1]],p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],nrow=5)
dev.off()

#画PC散点图（随机位点和适应性位点）
library(scatterpie)
library(RColorBrewer)
color <- brewer.pal(9, "Set1")


setwd("/Users/guoyafei/Desktop/群体结构")

for (i in c(3:7)) {
  infile <- paste("random.Q",i,".txt",sep="")
  outfile <- paste("random.Q",i,".pdf",sep="")
  data <- read.table(infile, header=T, stringsAsFactors = F)
  data <- as.data.frame(data)
  cols <- paste("Q",c(1:i),sep="")
  pdf(outfile,height=4,width=5)
  p <- ggplot() + 
    geom_scatterpie(aes(x=PC1, y=PC2), data=data,cols = cols, color=NA) + 
    coord_equal()+
    scale_fill_brewer(palette = "Set2")+
    theme_classic()
  print(p)
  dev.off()
  
}

for (i in c(3:7)) {
  infile <- paste("adap.Q",i,".txt",sep="")
  outfile <- paste("adap.Q",i,".pdf",sep="")
  data <- read.table(infile, header=T, stringsAsFactors = F)
  data <- as.data.frame(data)
  cols <- paste("Q",c(1:i),sep="")
  pdf(outfile,height=4,width=5)
  p <- ggplot() + 
    geom_scatterpie(aes(x=PC1, y=PC2), data=data,cols = cols, color=NA) + 
    coord_equal()+
    scale_fill_brewer(palette = "Set2")+
    theme_classic()
  print(p)
  dev.off()
  
}
  