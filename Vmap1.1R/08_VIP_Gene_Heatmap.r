library(pheatmap)
library(RColorBrewer)
display.brewer.all()
##Top5_R----
data2 <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/HeatMap/heatmap_format2.txt",header=T,
                    row.names= 1, stringsAsFactors=F,sep="\t")
data <- as.matrix(data2)
x=vector()
for ( i in c(1:87)){
  x <- c(x,strsplit(colnames(data), "type_")[[i]][2])
}
colnames(data)<-x
region <- c("Strang_WA", "WA_EU", "WA_CA", "WA_NW_A", "WA_NE_A", "WA_SA", "WA_SW_A", "WA_Tibet", "WA_SE_A", "CA_SA", "SA_SW_A", "SA_Tibet", "NW_A_NE_A", "SW_A_NE_A", "SE_A_NE_A")
data <- data[region,]
#非nlr基因的颜色
colsA = c("#F7FCFD","#8C96C6","#F7FCF5","#74C476","#FFFFCC","#FD8D3C")
pheatmap(data,cluster_rows = F,cluster_cols = T,border_color=NA,color = colsA,fontsize=8)
#nlr基因的颜色
data <- as.matrix(t(data2))
colsA = c("#006D2C","#a3d9b8","#88419D","#e9d5f0","#FEB24C","#f7e6d0")
pheatmap(data,cluster_rows = F,cluster_cols = F,border_color=NA,color = colsA,fontsize=8)
#P值
xlab <-factor(xlab,levels=c("Top5%","Top1%","Top0.5%","Top0.1%"))
P <- c(106.1767,69.862,65.09623,9.649433)
data <- as.data.frame(cbind(P,xlab))
ggplot(data = data, mapping = aes(x = xlab, y = P, group = 1,colour = "darkred")) + geom_line() + ylab('-log(P)')+xlab('Signal region ratio')+scale_x_discrete(limits= c("Top5%","Top1%","Top0.5%","Top0.1%"))+theme_classic()