##R: library prepare
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
#display.brewer.all()
setwd("/Users/guoyafei/Desktop/coXiao/TraesCS2A02G554200")
annotation_col <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3.info",header=T,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- annotation_col[,1]
#AB
fileAB <- c("8.297074952-297081886.vcf-geno","8.297075952-297080886.vcf-geno")
input <- paste(fileAB,".txt",sep="")
output <- paste(fileAB,".pdf",sep="")

for (i in c(1:length(input))) {
  all <- read.table(input[i],header=T,stringsAsFactors = F)
  #AB
  #colnames(all) <- c("ID","REF","ALT",paste("ABD_",c(0001:1143),sep=""),paste("AB_",c(001:212),sep=""),paste("ABD_",c(1144:1196),sep=""))
  data <- all[,4:1411]
  #D
  #colnames(all) <- c("ID","REF","ALT",paste("ABD_",c(0001:1143),sep=""),paste("D_",c(001:220),sep=""),paste("ABD_",c(1144:1196),sep=""))
  #data <- all[,4:1419]
  anno <- annotation_col[colnames(data),c(2,3,6),drop=FALSE]
  anno2 <- anno[which(anno$Common_name != "OtherHexaploids"),]
  
  anno2$type <- anno2$Common_name
  anno2[which(anno2$Common_name == "Polish_wheat" |anno2$Common_name == "Rivet_wheat" | anno2$Common_name == "Persian_wheat" |anno2$Common_name == "Khorasan_wheat"  |anno2$Common_name ==  "Durum"),4] <- "Freethreshing-Tetraploids"
  anno2[which(anno2$Common_name == "Wild_emmer_N"  |anno2$Common_name ==  "Wild_emmer_S" | anno2$Common_name == "Wild_emmer_SD"),4] <- "Wild_emmer"
  anno3 <- anno2[,c(1,3,4)]
  anno3$type <- factor(anno3$type,levels = c("Wild_emmer","Domesticated_emmer","Ispahanicum","Georgian_wheat","Freethreshing-Tetraploids","Spelt","Macha","Club_wheat","Tibetan_semi_wild", "Xinjiang_wheat", "Vavilovii","Indian_dwarf_wheat", "Yunan_wheat","Landrace","Cultivar"))
  anno3$Ctnt <- factor(anno3$Ctnt,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM"))
  
  anno4 <- anno3[order(anno3$type),]
  data2 <- data[,which(colnames(data) %in% rownames(anno4))]
  data3 <- data[,rownames(anno4)]
  cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
  pdf("test.pdf",width=150,height = 80 )
  pheatmap(data3, show_rownames=FALSE, show_colnames=T, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = T, cluster_row = FALSE, annotation_col = anno4, annotation_names_col = F)
  dev.off()
}


#画单倍型地理分布图
anno4$Latitude <- annotation_col[rownames(anno4),4]
anno4$Longitude <- annotation_col[rownames(anno4),5]


#准备文件（target.pie.txt）
library(cluster)
library(factoextra)
hap <- read.table("单倍型分类.txt",header=T,row.names=1,stringsAsFactors = F)
anno4$hap <- hap[rownames(anno4),1]
target <- anno4[which(anno4$type == "Landrace" & anno4$Latitude != "NA"),]

#根据经纬度给样本聚类
dataA <- target[,c(4,5)]
#df = scale(dataA,center = T,scale = T)
#colnames(df) <- colname
#按列进行标准化
#先求样本之间两两相似性
result <- dist(dataA, method = "euclidean")
#使用指定距离来计算数据矩阵行之间的距离
#euclidean：欧几里得距离
result_hc <- hclust(d = result, method = "ward.D2")
dataA$type <- cutree(result_hc, k=20)
lat_mean <- tapply(dataA[,1],dataA$type,mean,na.rm = TRUE)
lon_mean <- tapply(dataA[,2],dataA$type,mean,na.rm = TRUE)
dataA$cluster30_lat <- NA
dataA$cluster30_lon <- NA
for(i in 1:445) {
  for(j in 1:20){
    if(dataA[i,3] == j ){
      dataA[i,4] <- as.numeric(lat_mean[j])
      dataA[i,5] <- as.numeric(lon_mean[j])
    }
  }
}
write.table(dataA,"20_cluster.txt",sep="\t",row.names = F,quote=F)

target$lat_cluster <- dataA[rownames(target),4]
target$lon_cluster <- dataA[rownames(target),5]
write.table(target,"target.pie_20.txt", quote=F, sep="\t",row.names = T)

#画图
target <- read.table("target.pie.txt",header=T,row.names=1,sep="\t",stringsAsFactors = F)
library(reshape)
library(ggplot2)
library(RColorBrewer)
require (rworldmap)

pdf("pie2.pdf")
sample <- target[,c(7,8,6)]
sample$value <- 1
wheat_reshape <- cast(sample,lat_cluster+lon_cluster~hap) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-10,150),ylim=c(20,60),nameX="lon_cluster",nameY="lat_cluster",symbolSize=1.5,
              zColours=brewer.pal(8, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9])
legend(25,95,box.lty=0,bg="transparent","108 landrace GrowHabbit", col="black")
dev.off()

#等位基因频率
setwd("/Users/guoyafei/Desktop/coXiao/TraesCS2A02G554200")

data <- read.table("/Users/guoyafei/Desktop/coXiao/TraesCS2A02G554200/freq.snp2.txt", header=T,stringsAsFactors = F)
p <- vector()
cor <- vector()
for(i in c(2:39)){
  a <- cor.test(data[,i],data[,40])
  p <- append(p,a$p.value)
  cor <- append(cor,a$estimate)
}

sub <- data[c(-14,-17),]
for(i in c(2:39)){
  a <- cor.test(sub[,i],sub[,40])
  p <- append(p,a$p.value)
  cor <- append(cor,a$estimate)
}

all <- cbind(p,cor)
rownames(all) <- colnames(data)[2:39]
rownames(all) <- colnames(data)[2:39]

pdf("prec4-23pop.pdf")
ggplot(sub, aes(y = X8.297075044, x = bio_prec4)) +
  geom_point(size=3,alpha=0.8,color="#999999")+
  theme_classic()+
  geom_smooth(method = "lm", se=TRUE, 
              color="#0072B2", formula = y ~ x) 
dev.off()

ggplot(sub,aes(y = X8.297075044, x = bio_prec4)) + geom_point() + stat_smooth(method=lm)

pdf("temp2-23pop.pdf")
ggplot(data, aes(y = X8.297075044, x = bio_temp2)) +
  geom_point(size=3,alpha=0.8,color="#999999")+
  theme_classic()+
  geom_smooth(method = "lm", se=TRUE, 
              color="#0072B2", formula = y ~ x) 
dev.off()



#画25个群体的单倍型分布
setwd("/Users/guoyafei/Desktop/coXiao/TraesCS2A02G554200")
data <- read.table("25pop.pie.txt",header=T,stringsAsFactors = F)
pdf("pie3.pdf")
sample <- data[,c(8,9,6)]
sample$value <- 1
wheat_reshape <- cast(sample,mean_lat+mean_lon~hap) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-10,150),ylim=c(20,60),nameX="mean_lon",nameY="mean_lat",symbolSize=1.5,
        zColours=brewer.pal(8, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9])
legend(25,95,box.lty=0,bg="transparent","108 landrace GrowHabbit", col="black")
dev.off()


#画相关性最高的两个环境变量在群体里面的分布
library(foreign)
library(ggplot2)
library(plyr)
setwd("/Users/guoyafei/Desktop/coXiao/环境")

EA <- c(2,6,7,13,20)
CA <- c(5,8,11,12,25)
WA <- c(1,3,16,17,18)
SA <- c(4,10,19)
EU <- c(9,15,22,23)
AM <- c(21,24)
AF <- c(14)

#pop25
PC <- read.table("env.txt", header=T, stringsAsFactors = F)

#pop23
PC2 <- PC[which(PC$type != 21 & PC$type != 24),]
PC <- PC2

M <- aggregate(PC, by=list(PC$type),FUN=mean)
PC$type <- factor(PC$type, levels=M[order(M$prec4),][,3],ordered = TRUE)

#添加单倍型频率的数据
data <- read.table("/Users/guoyafei/Desktop/coXiao/TraesCS2A02G554200/freq.snp2.txt", header=T,stringsAsFactors = F)
sub <- data[,c(1,40)]
sub$type <- c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,25,3,4,5,6,7,8,9)

#pop23
sub2 <- sub[c(-14,-17),]
sub <- sub2
sub$type <- factor(sub$type, levels=M[order(M$prec4),][,3],ordered = TRUE)

PC$region <- NA
PC[which(PC$type %in% EA),dim(PC)[2]] <- "EA"
PC[which(PC$type %in% WA),dim(PC)[2]] <- "WA"
PC[which(PC$type %in% SA),dim(PC)[2]] <- "SA"
PC[which(PC$type %in% EU),dim(PC)[2]] <- "EU"
PC[which(PC$type %in% CA),dim(PC)[2]] <- "CA"
PC[which(PC$type %in% AM),dim(PC)[2]] <- "AM"
PC[which(PC$type %in% AF),dim(PC)[2]] <- "AF"
color <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#999999","#000000")
#乌拉尔图: "#999999"（灰色）, 野生一粒:"#E69F00"(橘黄色), 栽培一粒:"#56B4E9"(明兰色), AABBDD:"#009E73"(橄榄绿), AABB:"#F0E442"（明黄）,"#0072B2"(天蓝), EA："#D55E00"(橘红色), SCA："#CC79A7"(皮粉色))
pdf("prec4-23pop.pdf",width=6,height=3)
ggplot() +
  geom_boxplot(data=PC, aes(x=type,y = prec4,fill=region,color=region), alpha=0.8) +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  #scale_colour_discrete(breaks = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#999999","#000000"), labels = c('EU','WA','CA','SA','EA','AF','AM'))+
  ylab("bio_prec4") +
  xlab("Group") +
  geom_point(data=sub,aes(x = type, y = rescale(X8.297075044,c(0,160)),),size=2,alpha=0.8,color="#B22222")+
  scale_y_continuous(breaks=pretty_breaks(4),
                     sec.axis = sec_axis( ~rescale(.,c(0,0.643)),name = "Haplotype Frequency"))+
  
  #theme_classic()+
  #geom_smooth(method = "lm", se=TRUE, 
              #color="#0072B2", formula = y ~ x) 
  theme_classic()
dev.off()

