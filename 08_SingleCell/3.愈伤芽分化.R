#单倍型热图
#画单倍型图
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
cols <- c("#225EA8","#DEEBF7","#FEB24C","#BD0026")
setwd("/Users/yafeiguo/Documents/项目/拟南芥愈伤/haplotype/snpeff")
annotation <- read.table("../accession.anno.txt", header=T, stringsAsFactors = F, sep="\t")
rownames(annotation) <- annotation$AccessionID

anno <- annotation[,c(2:5,9),drop=FALSE]
anno <- anno[order(annotation$LatSeq),]

data <- read.table("/Users/yafeiguo/Documents/项目/拟南芥愈伤/haplotype/snpeff/AT1G74650-vcf-geno.txt", header=T,check.names = F, stringsAsFactors = F)
#-1,0,2,6
data <- read.table("/Users/yafeiguo/Documents/项目/拟南芥愈伤/haplotype/snpeff/AT3G11260-vcf-geno.txt", header=T,check.names = F, stringsAsFactors = F)
#-1,0,2
data <- read.table("/Users/yafeiguo/Documents/项目/拟南芥愈伤/haplotype/snpeff/AT4G21750-vcf-geno.txt", header=T,check.names = F, stringsAsFactors = F)
#-1,0,2,6
data <- read.table("/Users/yafeiguo/Documents/项目/拟南芥愈伤/haplotype/snpeff/AT5G66700-vcf-geno.txt", header=T,check.names = F, stringsAsFactors = F)
#-1,0,2,6
ann_color = list(
  heatmap = c(Wild_emmer = "#8C510A", Domesticated_emmer = "#DFC27D", Freethreshing = "#F6E8C3", EU="#66C2A5", WA= "#FC8D62",CA="#8DA0CB", EA ="#E78AC3",SA="#A6D854",Tibet="#FFD92F"))

#for(i in names(table(data$File))){
  #sub <- data[which(data$File==i),4:253]
  #file <- paste(i,".pdf",sep="")
  
  pdf(file)
  AB_anno <- anno[which(rownames(anno) %in% colnames(data)),,drop=F]
  data <- data[!apply(data, 1, function(x) any(grepl("2/2|3/3|4/4", x))), ]
  out <- data[,4:1138]
  colnames(out) <- paste("col",colnames(out),sep="")
  rownames(AB_anno) <- paste("col",rownames(AB_anno),sep="")
  sub <- out[, rownames(AB_anno)]
  sub <- as.data.frame(lapply(sub, as.integer))
  pheatmap(sub, show_rownames=FALSE, show_colnames=FALSE,color = cols,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = AB_anno,annotation_names_col = F)
  #dev.off()
#}

  
  
  
#单倍型地图
library(reshape)
library(ggplot2)
library(RColorBrewer)
require (rworldmap)
setwd("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/06_Structure/")
data <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/225.taxa.txt", header=T,row.names = 1, sep="\t", stringsAsFactors = F)
taxa <- read.table("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/06_Structure/20210829/Cluster_Location/cluster_70.txt",header=T,row.names = 1, stringsAsFactors = F)
for(i in names(table(data$File))){
  sub <- data[which(data$File==i),]
  file <- paste(i,".pdf",sep="")
  pdf(file)
  for( j in c(1:dim(sub)[1])){
    tit <- sub[j,3]
    sample <- sub[j,4:228]
    loc <- taxa[names(sample),4:5]
    loc$type <- as.numeric(sample[1,])
    loc$value <- 1
    wheat_reshape <- cast(loc,Latitude_70+Longitude_70~type) 
    wheat_reshape2 <- as.data.frame(wheat_reshape)
    name <- names(table(loc$type))
    if((j-1)%%4 !=0){
      mapPies(wheat_reshape2,xlim=c(-10,130),ylim=c(10,70),nameX="Longitude_70",nameY="Latitude_70",nameZs=name,symbolSize=1.5,
              zColours=brewer.pal(8, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9],main="tit")
      legend(25,95,box.lty=0,bg="transparent","108 landrace GrowHabbit", col="black")
    }else{
      par(mfrow=c(2,2),oma=c(1,1,1,1), mar=c(0,1,0,1), cex=1)
      mapPies(wheat_reshape2,xlim=c(-10,130),ylim=c(10,70),nameX="Longitude_70",nameY="Latitude_70",nameZs=name,symbolSize=1.5,
              zColours=brewer.pal(12, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9], main="tit")
      legend(25,95,box.lty=0,bg="transparent",tit, col="black")
    }
  }
  dev.off()
}

#1001拟南芥样本的地理位置分布

setwd("/Users/yafeiguo/Documents/项目/拟南芥愈伤/haplotype")
data <- read.table("accession.anno.txt",header=T,stringsAsFactors=F,sep="\t")
#R version 4.2.3
library(ggplot2)
library(sp)
library(maptools)
library(maps)
library(RColorBrewer)
#register_google(key = "AIzaSyAFnTwa7kEJAqIeHU-Kw6qPZiPUU-zXCLI")
visit.x <- as.numeric(data$Long)
visit.y <- as.numeric(data$Lat)
mp<-NULL #定义一个空的地图
mapworld<-borders("world", colour = "gray50", fill="white") #绘制基本地图
mp<-ggplot()+mapworld+ylim(-60,90)

#利用ggplot呈现，同时地图纵坐标范围从-60到90
mp2<-mp+geom_point(aes(x=visit.x,y=visit.y, size=data$count), color="darkorange") + scale_size(range=c(1,1))

#绘制带点的地图，geom_point是在地图上绘制点，x轴为经度信息，y轴为纬度信息，size是将点的大小按照收集的个数确定，color为暗桔色，scale_size是将点变大一些
mp3<-mp2+labs(x="Longitude", y="Latitude")+
  theme(legend.position = "none",panelx.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 25),axis.text.y = element_text(size = 25),panel.border = element_blank()) #将图例去掉

color <- brewer.pal(8, "Dark2")
color <- brewer.pal(8, "Dark2")[c(6)]

mp+geom_point(aes(x=data$Long, y=data$Lat,color =data$AdmixtureGroup), size=1.0)+
  scale_size(range=c(1,1))+
  theme_classic()+
  #scale_color_manual(values = color) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
       # legend.position = "null",
        panel.border = element_blank())

pdf("castAll.pdf",width=20,height=8)
mp3
dev.off()

#1001拟南芥样本的kmeans聚类分布






"#a6cee3"
"#1f78b4"
"#b2df8a"
"#33a02c"
"#fb9a99"
"#e31a1c"
"#fdbf6f"
"#ff7f00"
"#cab2d6"
"#fbb4ae"
"#b3cde3"
"#ccebc5"
"#decbe4"
"#fed9a6"
"#ffffcc"
"#e5d8bd"
"#fddaec"
"#f2f2f2"
"#8dd3c7"
"#ffffb3"
"#bebada"
"#fb8072"
"#80b1d3"
"#fdb462"
"#b3de69"
"#fccde5"
"#d9d9d9"

















