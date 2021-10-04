library(pheatmap)
library(RColorBrewer)
##Top1Up----
display.brewer.all()
brewer.pal(9,'BuPu')
cols = c("#e9d5f0","#88419D")
data2 <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/HeatMap/Top1Up.txt",header=T,
                   row.names= 1, stringsAsFactors=F,sep="\t")
data <- as.matrix(t(data2[,1:6]))
data = as.matrix(data)
pheatmap(data,cluster_rows = F,cluster_cols = F,border_color=NA,color = cols,fontsize=7)
#data11 = t(data)
#pheatmap(data11,cluster_rows = F,cluster_cols = F,color = cols,fontsize=6)
##Top5Up----
display.brewer.all()
brewer.pal(9,'YlOrRd')
cols = c("#f7e6d0","#FEB24C")
data <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/HeatMap/Top5Up.txt",header=T,
                   row.names= 1, stringsAsFactors=F,sep="\t")
data <- as.matrix(t(data2[,1:6]))
data = as.matrix(data)
pheatmap(data,cluster_rows = F,cluster_cols = F,border_color=NA,color = cols,fontsize=8)

##Top5----
display.brewer.all()
data2 <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/HeatMap/Top5.txt",header=T,
                   row.names= 1, stringsAsFactors=F,sep="\t")
data <- as.matrix(t(data2[,1:6]))
#A
A <- data[,1:16]
brewer.pal(9,'BuPu')
colsA = c("#e9d5f0","#88419D")
pheatmap(A,cluster_rows = F,cluster_cols = T,border_color=NA,color = colsA,fontsize=8)
#B
B <- data[,17:25]
brewer.pal(9,'Greens')
colsB = c("#a3d9b8","#006D2C")
pheatmap(B,cluster_rows = F,cluster_cols = T,border_color=NA,color = colsB,fontsize=8)
#D
D <- data[,26:36]
brewer.pal(9,'YlOrRd')
colsD = c("#f7e6d0","#FEB24C")
pheatmap(D,cluster_rows = F,cluster_cols = T,border_color=NA,color = colsD,fontsize=8)
par(mfcol=c(1,3)) 

##Top5_R----
display.brewer.all()
data2 <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/HeatMap/Top5_R_new_noEA.txt",header=T,
                    row.names= 1, stringsAsFactors=F,sep="\t")
data <- as.matrix(data2)
#非nlr基因的颜色
colsA = c("#F7FCFD","#8C96C6","#F7FCF5","#74C476","#FFFFCC","#FD8D3C")
pheatmap(data,cluster_rows = F,cluster_cols = F,border_color=NA,color = colsA,fontsize=8)
#nlr基因的颜色
data <- as.matrix(t(data2))
colsA = c("#006D2C","#a3d9b8","#88419D","#e9d5f0","#FEB24C","#f7e6d0")
pheatmap(data,cluster_rows = F,cluster_cols = F,border_color=NA,color = colsA,fontsize=8)
#P值
xlab <-factor(xlab,levels=c("Top5%","Top1%","Top0.5%","Top0.1%"))
P <- c(106.1767,69.862,65.09623,9.649433)
data <- as.data.frame(cbind(P,xlab))
ggplot(data = data, mapping = aes(x = xlab, y = P, group = 1,colour = "darkred")) + geom_line() + ylab('-log(P)')+xlab('Signal region ratio')+scale_x_discrete(limits= c("Top5%","Top1%","Top0.5%","Top0.1%"))+theme_classic()