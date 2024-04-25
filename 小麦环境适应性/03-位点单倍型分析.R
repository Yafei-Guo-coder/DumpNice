################################################# 画单倍型热图 ###################################
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

###################################################### 画单倍型地理分布图 ##########################################
setwd("/Users/guoyafei/Desktop/coXiao/TraesCS2A02G554200")
data <- read.table("单倍型分类_V1.0.txt", header=T, sep = "\t", stringsAsFactors = F)

#准备文件（target.pie.txt）
library(cluster)
library(factoextra)

target <- data[which(data$Type != "Cultivar" & data$Latitude != "NA"),]

#根据经纬度给样本聚类
dataA <- target[,c(2,3)]
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
for(i in 1:136) {
  for(j in 1:20){
    if(dataA[i,3] == j ){
      dataA[i,4] <- as.numeric(lat_mean[j])
      dataA[i,5] <- as.numeric(lon_mean[j])
    }
  }
}
target$lat_cluster <- dataA[rownames(target),4]
target$lon_cluster <- dataA[rownames(target),5]
write.table(target,"V1.0.pie_20.txt", quote=F, sep="\t",row.names = F)

#画图
target <- read.table("/Users/guoyafei/Desktop/coXiao/TraesCS2A02G554200/V1.0.pie_20.txt",header=T,row.names=1,sep="\t",stringsAsFactors = F)
library(reshape)
library(ggplot2)
library(RColorBrewer)
require (rworldmap)

pdf("pie2.pdf")
sample <- target[,c(6,7,5)]
sample$value <- 1
wheat_reshape <- cast(sample,lat_cluster+lon_cluster~Haplotype) 
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

###################################################### haplotype network #################################################
setwd("/Users/guoyafei/Desktop/coZuo")
suppressMessages(library(dplyr))
library(gdata)
library (geosphere)
library(reshape)
library(ggplot2)
suppressMessages(library(stringr))
suppressMessages(library(pegas))
alignment <- read.dna("V1.37exonPos.21chr.beagle.min4.fasta",format = "fasta")
pal_269 <- c( "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
              "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
              "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
              "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
              "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
              "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
              "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
              "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
              
              "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
              "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
              "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
              "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
              "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
              "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
              "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
              "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
              
              "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
              "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
              "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
              "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
              "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
              "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
              "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
              "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",
              
              "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
              "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
              "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
              "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
              "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
              "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
              "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
              "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
              
              "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
              "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B")

samples <- samples %>% filter(ID %in% labels(alignment))
samples <- droplevels(samples);

# otherwise, the haplotype table will be screwy
rownames(samples) <- 1:nrow(samples)

# make the data table be in the same order as the aligment (or everything breaks)
samples <- samples[match(labels(alignment),samples$ID),]
hap <- haplotype(alignment,strict=F)

# Calculate haplotype network
rownames(hap) <-  paste0("hap_", seq(1, length(rownames(hap)), by=1) )
hap.net <- haploNet(hap)

hap.pies <- with(
  stack(setNames(attr(hap,'index'),1:length(attr(hap,'index')))),
  table(hap=as.numeric(as.character(ind)),pop=samples[values,"Common.name"])
)



rownames(hap.pies) <-  paste0("hap_",seq(1, length(rownames(hap)), by=1) )
hap.pies <- hap.pies[,c("Wild einkorn","Domesticated einkorn","Wild emmer","Domesticated emmer","Landrace","Cultivar")]


pal <- pal_269

if (length(pal) > ncol(hap.pies)) {
  pal <- pal[1:ncol(hap.pies)]
}
pdf(file = "1A.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7)

plot(hap.net, size=attr(hap.net, "freq")*0.01, bg=pal,
     scale.ratio =10, cex = 0.7, labels=T,
     pie=hap.pies, font=2, fast=F, legend =T, show.mutation=T, threshold=0)









