#--------------manhuttun-------------
library(qqman)
library(tidyverse)
library(ggrepel)
#setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/09_GWAS/LD01")
#path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/09_GWAS/LD01/Manhuttan" ##文件目录

setwd("/Users/guoyafei/Desktop/毕业论文/growthHabit/V2")
path <- "/Users/guoyafei/Desktop/毕业论文/growthHabit/V2/fastGWA" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})

path2 <- "/Users/guoyafei/Desktop/毕业论文/growthHabit/V2/gene" ##文件目录
fileNames2 <- dir(path2)  ##获取该路径下的文件名
filePath2 <- sapply(fileNames2, function(x){ 
  paste(path2,x,sep='/')})   ##生成读取文件路径
data2 <- lapply(filePath2, function(x){
  read.table(x, header=F,stringsAsFactors = F)})

path3 <- "/Users/guoyafei/Desktop/毕业论文/growthHabit/V2/anno" ##文件目录
fileNames3 <- dir(path3)  ##获取该路径下的文件名
filePath3 <- sapply(fileNames3, function(x){ 
  paste(path3,x,sep='/')})   ##生成读取文件路径
data3 <- lapply(filePath3, function(x){
  read.table(x, header=F,stringsAsFactors = F)})
#threshold
library(gdata)
thresh <- read.xls("thresh.xlsx",sheet=1,row.name=1,na.strings=c("NA","#DIV/0!"))

for (i in c(1:3)){
  #for (i in seq(1,40,by=2)){
  filename <- strsplit( names(data)[i], "_")[[1]][1]
  p <- list()
  #for (j in c(0,1)){
  #a <- i+j
  gwasResults <- data[[i]]
  point <- data2[[i]]
  anno <- data3[[i]]
  colnames(gwasResults) <- c("SNP", "CHR", "BP","P")
  #gwasResults3 <- gwasResults2[which(gwasResults2$P < 0.25),]
  #other <- gwasResults2[which(gwasResults2$P > 0.01 & gwasResults2$P <0.2),]
  #other2 <- other[sample(nrow(other), 20000), ]
  #gwasResults <- rbind(gwasResults3,other)
  don <- gwasResults %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot)
  don$highlight <- "no"
  don[which(don$SNP %in% point[,4]), 7] <- "yes"
  
  don$anno <- "no"
  rownames(don) <- don$SNP
  don[anno$V4, 8] <- "yes"
  
  don$name <- "no"
  don[anno$V4, 9] <- anno[,6]
  
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  p[[i]] <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Add highlighted points
    geom_point(data=subset(don, highlight=="yes"), color="orange", size=1.3) +
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, anno=="yes"), aes(label=name), size=1, max.overlaps=200) +
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=15),
      axis.text.y=element_text(size=15),
      axis.title.y=element_text(size = 15),
      axis.title.x=element_text(size = 15),
    )+
    scale_y_continuous(limits = c(0,10))
  #geom_point(data=point,aes(x=BPcum,y=-log10(P)),color="red")
  #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
  #geom_hline(yintercept = -log10(thresh[a,1]), colour="red",linetype=2, size=1)+
  #geom_hline(yintercept = -log10(thresh[a,2]), colour="blue",linetype=2, size=1)
  #}
  pdf(paste(filename,".pdf",sep=""),height = 4.5,width = 9)
  #grid.arrange(p[[1]],p[[2]],nrow=2)
  print(p[[i]])
  dev.off()
}

#--------------Accession Map-------------
library(ggmap)
library(RColorBrewer)
setwd("/Users/guoyafei/Desktop/毕业论文/growthHabit")
taxa <- read.table("taxa_withGeo.txt", header = T, stringsAsFactors = F)

taxa$Latitude <- as.numeric(taxa$Latitude)
taxa$Longitude <- as.numeric(taxa$Longitude)

#绘图:全部的
mp <- NULL
mapworld <- borders("world",colour = "gray90",fill="gray90") 
mp <- ggplot() + 
  mapworld + 
  #xlim(0,90) +
  ylim(-60,90) + 
  theme_classic()

color <- brewer.pal(8, "Dark2")[c(4,6,3)]
color <- brewer.pal(8, "Dark2")[c(6)]

mp+geom_point(aes(x=taxa$Longitude, y=taxa$Latitude,color = taxa$GrowthHabit), size=1.5)+
  scale_size(range=c(1,1))+
  scale_color_manual(values = color) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank())  

#--------------显著性检验------
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
brewer.pal(9, "Set3")
setwd("/Users/guoyafei/Desktop/毕业论文/growthHabit/显著性检验")
data <- read.table("permu400.snp.txt", header=F, stringsAsFactors = F)

ggplot(data, aes(x = V2)) +
  geom_histogram(binwidth = 1, fill = "lightblue", colour = "black")+
  theme_classic()

ggplot(data, aes(y = V1, group = V2)) +
  geom_histogram(binwidth = 1, fill = "lightblue", colour = "black")+
  theme_classic()+
  ylim(1,10)

ggplot(data,aes(x=V1,group=V2,fill=V2)) +
  geom_density()+
  scale_fill_brewer(palette = "Set3")+
  theme_classic()+
  geom_vline(xintercept = 11, color = '#8DD3C7', size = 0.5) + 
  geom_vline(xintercept = 28, color = '#FFFFB3', size = 0.5) + 
  geom_vline(xintercept = 25, color = '#BEBADA', size = 0.5) + 
  geom_vline(xintercept = 14, color = '#FB8072', size = 0.5) 

ggplot(data,aes(x=V1,group=V2,fill=V2)) +
  geom_density()+
  #scale_fill_brewer(palette = "Set3")+
  theme_classic()+
  geom_vline(xintercept = 5, color = '#8DD3C7', size = 0.5) + 
  geom_vline(xintercept = 139, color = '#FFFFB3', size = 0.5) + 
  geom_vline(xintercept = 44, color = '#BEBADA', size = 0.5) + 
  geom_vline(xintercept = 2, color = '#FB8072', size = 0.5) 

################# VIP4 gene haplotype distribution #########################
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
cols <- c("#225EA8","#DEEBF7","#FEB24C","#BD0026")
setwd("/Users/guoyafei/Desktop/毕业论文/growthHabit/VIP4")

data <- read.table("VIP4.geno.txt", header= T,stringsAsFactors = F,check.names = F)
geo <- read.table("../taxa_withGeo.txt", header=T, stringsAsFactors = F)
split_list <- strsplit(geo$ID,"_")
second_elements <- sapply(split_list, `[`, 2)
rownames(geo) <- second_elements
anno <- geo[,7,drop=FALSE]
ann_color = list(heatmap = c(Wild_emmer = "#8C510A", Domesticated_emmer = "#DFC27D", Freethreshing = "#F6E8C3", EU="#66C2A5", WA= "#FC8D62",CA="#8DA0CB", EA ="#E78AC3",SA="#A6D854",Tibet="#FFD92F"))
all <- data[, rownames(anno)]
pheatmap(all, show_rownames=FALSE, show_colnames=FALSE,color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = anno, annotation_names_col = F)


data <- read.table("VIP4.geno.all.txt", header= T,stringsAsFactors = F,check.names = F)
geo <- read.table("../taxa_withGeo.all.txt", header=T, stringsAsFactors = F)
rownames(geo) <- geo$ID
anno <- geo[,c(3,7),drop=FALSE]
all <- data[, rownames(anno)]
pheatmap(all, show_rownames=FALSE, show_colnames=FALSE,color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = anno, annotation_names_col = F)









