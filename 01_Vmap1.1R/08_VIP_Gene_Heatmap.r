library(pheatmap)
library(RColorBrewer)
display.brewer.all()
##Top5 VIP gene & Top1 NLR gene
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene")
data2 <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/heatmap_format2.txt",
                    header=T,
                    row.names= 1, stringsAsFactors=F,sep="\t")
data2 <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Xpclr/V3/Top1nlr/heatmap_format2_changeName.txt",
                    header=T,
                    row.names= 1, stringsAsFactors=F,sep="\t")
data <- as.matrix(data2)
x=vector()
for ( i in c(1:dim(data)[2])){
  x <- c(x,strsplit(colnames(data), "type_")[[i]][2])
}
colnames(data)<-x
#region <- c("Strang_WA", "WA_EU", "WA_CA", "WA_NW_A", "WA_NE_A", "WA_SA", "WA_SW_A", "WA_Tibet", "WA_SE_A", "CA_SA", "SA_SW_A", "SA_Tibet", "NW_A_NE_A", "SW_A_NE_A", "SE_A_NE_A")
#region <- c("neg_North1_North2", "neg_WA_North1", "neg_WA_North2", "WA_EU", "Strang_WA","WA_South","Tibet_South", "North2_South")
region <- c( "Strang_WA","WA_EU","WA_South","Tibet_South", "North2_South","EU_South")
region <- c("Strang_WA","WA_EU","WA_South","neg_WA_North1","neg_WA_North2")
data <- data[region,]
colsA = c("#F7FCFD","#8C96C6","#F7FCF5","#74C476","#FFFFCC","#FD8D3C")
pheatmap(data,cluster_rows = F,cluster_cols = T,border_color=NA,color = colsA,fontsize=8)

#画top5% nlr基因的曼哈顿图----
#204@xuebo:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%/nlr
#cat neg_WA_North2_smooth_* > manhattan/neg_WA_North2.txt
#cd manhattan/
#ls neg_WA_North2_smooth_*bed
#cat neg_WA_North2_smooth_*bed > nlr/manhattan/neg_WA_North2_top5.bed
#grep -f neg_WA_North2.txt ../../../nlr_wheat.bed | awk '{print $1"\t"$2"\t"$3}' > neg_WA_North2_gene.bed
#for i in `cat names`
#do
#bedtools intersect -a ${i}_top5_21chr.bed -b ${i}_gene.bed -wo | awk '{print $7"\t"$8"\t"$9"\t"$5}' |sed 's/chr//' |sort > ${i}.manhattan.txt
#done
library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Xpclr/V3/Top1nlr/manhattan")
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Xpclr/V3/Top1nlr/manhattan" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})
tit <- c("EU_South","WA_North2","North2_South","Strang_WA","Tibet_South","WA_EU","WA_South")
pdf("nlr.xpclr.pdf",width = 7.5,height = 4)
for(i in c(1:length(data))){
  data1 <- data[[i]]
  data1$SNP <- c(1:dim(data[[i]])[1])
  gwasR <- data1[,c(5,1,2,4)]
  colnames(gwasR) <- c("SNP", "CHR", "BP","P")
  # 1)计算chr长度
  chr_len <- tibble(
    CHR=1:21,
    chr_len= c(594102056,689851870,495453186,780798557,801256715,651852609,750843639,830829764,615552423,744588157,673617499,509857067,709773743,713149757,566080677,618079260,720988478,473592718,736706236,750620385,638686055)
  )
  # 2)计算每条chr的初始位置
  chr_pos <- chr_len  %>% 
    mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) 
  # 3)计算累计SNP的位置
  Snp_pos <- chr_pos %>%
    left_join(gwasR, ., by="CHR") %>%
    arrange(CHR, BP) %>%
    mutate( BPcum = BP + total)
  
  #1) 准备X轴标签位置--在每条chr的中间
  #X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )
  start <- chr_pos$total
  stop <- c(594102056,1283953926,1779407112,2560205669,3361462384,4013314993,4764158632,5594988396,6210540819,6955128976,7628746475,8138603542,8848377285,9561527042,10127607719,10745686979,11466675457,11940268175,12676974411,13427594796,14066280851)
  X_axis <- tibble(
    CHR = 1:21,
    center = (start+stop)/2)
  #2）绘制“改良版”Manhattan图
  p1 <- ggplot(Snp_pos, aes(x=BPcum, y=P)) +
    geom_point( aes(color=as.factor(CHR)), size=1) +
    scale_color_manual(values = rep(c( "orange"), 21 )) +
    scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
    theme_bw() +
    theme(
      legend.position="none",
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=10),
      axis.text.y=element_text(size=10),
      axis.title.y=element_text(size = 10),
      axis.title.x=element_text(size = 10),
    )+
    ylab('XP-CLR ratio')+xlab('Chromosome Position')+
    ggtitle(tit[i])+
    theme(plot.title = element_text(color="red", size=15, face="bold"))
  print(p1)
}
dev.off()
#P值----
xlab <-factor(xlab,levels=c("Top5%","Top1%","Top0.5%","Top0.1%"))
P <- c(106.1767,69.862,65.09623,9.649433)
data <- as.data.frame(cbind(P,xlab))
ggplot(data = data, mapping = aes(x = xlab, y = P, group = 1,colour = "darkred")) + geom_line() + ylab('-log(P)')+xlab('Signal region ratio')+scale_x_discrete(limits= c("Top5%","Top1%","Top0.5%","Top0.1%"))+theme_classic()
