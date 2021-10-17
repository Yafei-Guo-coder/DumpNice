library(pheatmap)
library(RColorBrewer)
display.brewer.all()
##Top5_R----
data2 <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/HeatMap/V3/Top5_heatmap_format3.txt",
                    header=T,
                    row.names= 1, stringsAsFactors=F,sep="\t")
data <- as.matrix(data2)
x=vector()
for ( i in c(1:47)){
  x <- c(x,strsplit(colnames(data), "type_")[[i]][2])
}
colnames(data)<-x
#region <- c("Strang_WA", "WA_EU", "WA_CA", "WA_NW_A", "WA_NE_A", "WA_SA", "WA_SW_A", "WA_Tibet", "WA_SE_A", "CA_SA", "SA_SW_A", "SA_Tibet", "NW_A_NE_A", "SW_A_NE_A", "SE_A_NE_A")
#region <- c("neg_North1_North2", "neg_WA_North1", "neg_WA_North2", "WA_EU", "Strang_WA","WA_South","Tibet_South", "North2_South")
region <- c("neg_WA_North2", "Strang_WA","WA_EU","WA_South","Tibet_South", "North2_South")
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

#Go analysis
#Working directory
#xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage/Top5%/Go

for i in A B D
do
for j in `cat names`
do
grep -v -f neg_WA_North2_smooth_${i}.top5.txt ${j}_smooth_${i}.top5.txt > ${j}.go.gene.txt
done
done



