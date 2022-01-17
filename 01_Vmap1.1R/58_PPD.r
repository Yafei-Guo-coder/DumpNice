#PPD-D1
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/10_Gene")
data<- read.table("all.bin.depth.txt",header=T,row.names = 1,stringsAsFactors = F)
data <- read.table("ppd_331taxa.txt",header=T,stringsAsFactors = F)
sub <- data
rownames(data) <- round(data[,1])
sub <- data[,2:332]
names <- colnames(sub)
for(i in c(1:331)){
    sub[which(sub[,i] >0 & sub[,i] <10),i] <- 10
    sub[which(sub[,i] >10 & sub[,i] <30),i] <- 25
    sub[which(sub[,i] >30),i] <- 50
  }
annotation <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/select_taxa4.txt", header=T,  stringsAsFactors = F, sep="\t")
anno <- annotation[,c(2,9),drop=FALSE]

AB_anno <- anno[which(anno$VMap3 %in% colnames(sub)),]

rownames(AB_anno) <- AB_anno$VMap3

AB_anno <- AB_anno[,2,drop=F]
AB_anno <- AB_anno[!is.na(AB_anno$heatmap),1,drop=F]
labels_row =  rownames(sub)
labels_row[1:150] <- ""
labels_row[361:511] <- ""
labels_row[152:359] <- ""
ann_color = list(
  #Growing_Habit = c(Facultative = "yellow", Spring="orange", Winter="blue"),
  heatmap = c(Strangulata = "#8C510A", EU="#66C2A5", WA= "#FC8D62",CA="#8DA0CB", EA ="#E78AC3",SA="#A6D854",Tibet="#FFD92F", Cultivar="#999999"))
sub <- sub[, rownames(AB_anno)]
labels_row <- rep("",times=32)
labels_row[28] <- "16bp"
labels_row[29] <- "5bp"
labels_row[32] <- "2k"
pheatmap(sub, show_rownames=T, labels_row=labels_row,show_colnames=F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), annotation = AB_anno,cluster_col = F,annotation_names_col = F,annotation_colors = ann_color, cluster_row = FALSE) 

#Haplotype Map: PPD-D1
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/10_Gene")
data <- read.table("ppd_331taxa.txt", header = T, stringsAsFactors = F)
taxa <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/10_Gene/cluster_70.txt",header=T,row.names = 1, stringsAsFactors = F)
taxa <- taxa[!is.na(taxa$VMap3),]
rownames(taxa) <- taxa[,1]
pdf("ppd.pdf")
for( j in c(1:dim(sub)[1])){
    loc <- taxa[,5:6]
    loc$type <- as.numeric(data[j,rownames(taxa)])
    loc$value <- 1
    wheat_reshape <- cast(loc,Latitude_70+Longitude_70~type) 
    wheat_reshape2 <- as.data.frame(wheat_reshape)
    name <- names(table(loc$type))
    if((j-1)%%4 !=0){
      mapPies(wheat_reshape2,xlim=c(-10,130),ylim=c(10,70),nameX="Longitude_70",nameY="Latitude_70",nameZs=name,symbolSize=1.5,
              zColours=brewer.pal(8, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9],main="tit")
    }else{
      par(mfrow=c(2,2),oma=c(1,1,1,1), mar=c(0,1,0,1), cex=1)
      mapPies(wheat_reshape2,xlim=c(-10,130),ylim=c(10,70),nameX="Longitude_70",nameY="Latitude_70",nameZs=name,symbolSize=1.5,
              zColours=brewer.pal(12, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9], main="tit")
    }
    }
dev.off()

#fst
library(ggplot2)
data <- read.table("topP.txt1",header=T,stringsAsFactors=F)
col3 <- colorRampPalette(c("white","blue"))
data$pop1 <- factor(stbp$V5, levels = c("EU","WA","CA","EA","SA"))
data$pop2 <- factor(stbp$V6, levels = c("EU","WA","CA","EA","SA"))

pdf("gene.pdf",height = 6,width = 7)
ggplot(data = data) + 
  geom_tile(aes(pop1, pop2, fill = gene)) +
  scale_fill_gradientn(colours = col3(100), limits = c(-0.1, 0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("")+
  ylab("")
dev.off()

pdf("16bp.pdf",height = 6,width = 7)
ggplot(data = data) + 
  geom_tile(aes(pop1, pop2, fill = sixbp)) +
  scale_fill_gradientn(colours = col3(100), limits = c(-0.1, 0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("")+
  ylab("")
dev.off()

pdf("5bp.pdf",height = 6,width = 7)
ggplot(data = data) + 
  geom_tile(aes(pop1, pop2, fill = fivebp)) +
  scale_fill_gradientn(colours = col3(100), limits = c(-0.1, 0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("")+
  ylab("")
dev.off()

pdf("2k.pdf",height = 6,width = 7)
ggplot(data = data) + 
  geom_tile(aes(pop1, pop2, fill = twobp)) +
  scale_fill_gradientn(colours = col3(100), limits = c(-0.1, 0.5))+
  theme_classic()+
  
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("")+
  ylab("")
dev.off()
