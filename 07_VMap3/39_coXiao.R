#plot
setwd("/Users/guoyafei/Desktop/coXiao/单倍型")
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)


########################################################## 单倍型分析 ####################################################

data <- read.table("/Users/guoyafei/Desktop/coTong/7.177767585-460000000.shuf1k.txt", header=T, check.names=F, stringsAsFactors = F)

#annotation_col <- read.table("/Users/guoyafei/Desktop/coXiao/单倍型/haploinfo.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")
annotation_col <- read.table("/Users/guoyafei/Desktop/coTong/haploinfo_V3.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- annotation_col$ID
#anno <- annotation_col[which(rownames(annotation_col) %in% colnames(data)),c(3,4,5),drop=FALSE]
anno <- annotation_col[which(rownames(annotation_col) %in% colnames(data)),c(2,4,5),drop=FALSE]
#anno2 <- anno[,3,drop=F]
anno2 <- anno[,2,drop=F]
#anno2$type <- factor(anno$`Common name`,levels = c("Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Polish wheat","Durum","Rivet wheat","Khorasan wheat","Persian wheat","Spelt","Indian dwarf wheat","Tibetan semi wild","Macha","Club wheat","Vavilovii","Yunan wheat","Xinjiang wheat","Landrace","Cultivar"))
anno2$type <- factor(anno$type_final,levels = c("Wild_emmer","Domesticated_emmer","Free_threshing_tetraploids","Landrace","Cultivar"))
anno2$Ctnt <- factor(anno$`Ctnt(9)`,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM","-"))

#anno2$Ctnt <- factor(anno$Region,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM"))

#anno4 <- anno3[order(anno3$type,anno3$Ctnt),]

### extract data & plot ------>
data2 <- data[,which(colnames(data) %in% rownames(anno))]
data3 <- data2[,rownames(anno)]
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
pdf("/Users/guoyafei/Desktop/coTong/test_tong_region2_all.pdf",width=12,height = 8 )
pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, cluster_col = F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"),cluster_row = FALSE, annotation_col = anno2[,c(2),drop=F], annotation_names_col = F)
dev.off()

#haploinfo2
annotation_col <- read.table("haploinfo2.txt",header=T,check.names=F,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- annotation_col$ID
anno <- annotation_col[which(rownames(annotation_col) %in% colnames(data)),c(4,5,6),drop=FALSE]

anno2 <- anno[,3,drop=F]
anno2$type <- factor(anno$type,levels = c("Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","free-threshing tetraploids","Spelt","Indian dwarf wheat","Tibetan semi wild","Macha","Club wheat","Vavilovii","Yunan wheat","Xinjiang wheat","Landrace","Cultivar"))
anno2$Ctnt <- factor(anno$Region,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM"))

#anno4 <- anno3[order(anno3$type,anno3$Ctnt),]

### extract data & plot ------>
data2 <- data[,which(colnames(data) %in% rownames(anno))]
data3 <- data2[,rownames(anno)]
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
pdf("test3.pdf",width=60,height = 40)
pheatmap(data3, show_rownames=FALSE, show_colnames=T, cluster_col = F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"),cluster_row = FALSE, annotation_col = anno2, annotation_names_col = F)
dev.off()

########################################################## 渗入分析 #################################################

library(ggplot2)
setwd("/Users/guoyafei/Desktop/coXiao/渗入/六倍体/")
#data <- read.table("chr2A.withAnc-1.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-All.csv",header=T,stringsAsFactors = F, sep=",")

#data <- read.table("chr2A.withAnc-D-1.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-D-All.csv",header=T,stringsAsFactors = F, sep=",")

#data <- read.table("chr2A.withAnc-U-1.csv",header=T,stringsAsFactors = F, sep=",")
#data <- read.table("chr2A.withAnc-U-All.csv",header=T,stringsAsFactors = F, sep=",")

setwd("/Users/guoyafei/Desktop/coXiao/渗入/二倍体/")
data <- read.table("chr2A.withAnc-All_D.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-All_F.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("/Users/guoyafei/Desktop/coXiao/渗入/二倍体/chr2A.withAnc-All_W.csv",header=T,stringsAsFactors = F, sep=",")

data <- read.table("chr2A.withAnc-D-All_D.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-D-All_F.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-D-All_W.csv",header=T,stringsAsFactors = F, sep=",")

data <- read.table("chr2A.withAnc-All_W_P1D.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-All_D_P1D.csv",header=T,stringsAsFactors = F, sep=",")
data <- read.table("chr2A.withAnc-All_F_P1D.csv",header=T,stringsAsFactors = F, sep=",")

data <- data[which(data$fd >= 0 & data$fd < 1 & data$D > 0),]
#data[which(data$fd < 0 | data$fd > 1 | data$D <= 0),10] <- 0
quantile(data$fd,0.95)
#0.317498 # 2.41M
sub <- data[which(data$start > 740000000),]
ggplot(sub, aes(x=start, y=fd)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("fd") +
  geom_hline(yintercept=0.0721 ,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759452125,color='red')

#0.3167 #2.36M
sub <- data[which(data$start > 740000000),]
ggplot(sub, aes(x=start, y=fd)) +
  geom_line()+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("fd") +
  geom_hline(yintercept=0.3167,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759452125,color='red')

######################################################### 选择和适应性分析 ##################################################
setwd("/Users/guoyafei/Desktop/coXiao/选择")
data <- read.table("prec4.chr2A.txt",header=F,stringsAsFactors = F, sep="\t")

pdf("baypass.pdf",width=14,height = 6)
ggplot(data, aes(x=V2, y=V5)) +
  geom_point(alpha=0.5)+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("fd") +
  #geom_hline(yintercept=0.3167,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759452125,color='red')
dev.off()
  
data <- read.table("SH_IA.chr2A.21chr.txt",header=F,stringsAsFactors = F, sep="\t")
pdf("xpclr.pdf",width=16,height = 6)
ggplot(data, aes(x=V2, y=V4)) +
  geom_line(linewidth=0.5)+
  theme_classic()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("fd") +
  geom_hline(yintercept=2.282,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759452125,color='red')
dev.off()

############################################################### 地图分析 #####################################################
library(leaflet)
library(rgdal)
library(RColorBrewer)
setwd("/Users/guoyafei/Desktop/coXiao/选择/")
tibet <- read.table("/Users/guoyafei/Desktop/coXiao/选择/all.ibs.3.txt",sep="\t",header=T,stringsAsFactors = F)
data <- na.omit(tibet)
type_var <- c("Domesticated einkorn","Wild einkorn","free-threshing tetraploids")
data <- data[which(data$type %in% type_var),]
point <- tibet[which(tibet$type == "Landrace" & tibet$TE == "1"  ),]
library(ggplot2)
library(RColorBrewer)
library(ggmap)
library(sp)
library(maptools)
library(maps)
library(psych)
library(reshape)
library("corrplot")

WE <- read.table("WE.txt", header=F,stringsAsFactors = F)
#类型1
mp <- NULL
mapworld <- borders("world",colour = "grey90",fill="white") 
mp <- ggplot()+mapworld+ylim(0,75)+xlim(-25,70)
mp2 <- mp+
  geom_point(data=data,aes(x=data$Longitude,y=data$Latitude,color=data$ibs,shape= data$type),size=3)+
  scale_size(range=c(1,1))+ 
  scale_colour_gradientn(colours = brewer.pal(11, "RdYlBu"))+
  scale_shape_manual(values = c(15, 16, 17))+
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  theme(panel.grid = element_blank()) + 
  theme(legend.title=element_blank())+
  xlab("Lon")+
  ylab("Lat")+
  
  theme(legend.text = element_text(size = 15))+
  
  #theme(panel.grid =element_blank()) +   ## 删去网格线
  #theme(axis.text = element_blank()) +   ## 删去刻度标签
  #theme(axis.ticks = element_blank()) +   ## 删去刻度线
  geom_point(data=point,aes(x=point$Longitude,y=point$Latitude),color="#000000",shape=8,size=3 )+
  geom_point(data=WE,aes(x=WE$V3,y=WE$V2),color="#000000",shape=10,size=3 )
pdf("ibs2.pdf",height = 8,width=16)
print(mp2)
dev.off()
#类型2
ggplot(data[[1]], aes(Logititude, Latitude))+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(color=IBS,shape=Type),size=1.5)+
  #geom_point(aes(color=mean,shape=Value),size=3)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))


############################################################ 重新画地图 ######################################3#############
library(rgdal)
library(RColorBrewer)
library(ggplot2)
library(ggmap)
library(sp)
library(maptools)
library(maps)
library(psych)
library(reshape)
library(leaflet)
library(rgdal)

setwd("/Users/guoyafei/Desktop/coXiao/选择/")
tibet <- read.table("/Users/guoyafei/Desktop/coXiao/选择/all.ibs.3.txt",sep="\t",header=T,stringsAsFactors = F)
data <- na.omit(tibet)
type_var <- c("Domesticated einkorn","Wild einkorn","free-threshing tetraploids")
data <- data[which(data$type %in% type_var),]
point <- tibet[which(tibet$type == "Landrace" & tibet$TE == "1"  ),]

WE <- read.table("WE.txt", header=F,stringsAsFactors = F)
colnames(WE) <- c("type","Latitude","Longitude")

mypalette <- colorBin( palette="YlOrBr", domain=data$ibs, na.color="transparent")
pchIcons <- function(pch = 0:14, width = 30, height = 30, ...) {
  n <- length(pch)
  files <- character(n)
  # create a sequence of png images
  for (i in seq_len(n)) {
    f <- tempfile(fileext = ".png")
    png(f, width = width, height = height, bg = "transparent")
    par(mar = c(0, 0, 0, 0))
    plot.new()
    points(.5, .5, pch = pch[i], cex = min(width, height) / 8, ...)
    dev.off()
    files[i] <- f
  }
  files
}

shapes <- sample(0:14, 10)
iconFiles <- pchIcons(shapes, 40, 40, col = "steelblue", lwd = 2)

data$group <- NA
data[which(data$type == "Wild einkorn"),9] <- 2
data[which(data$type == "Domesticated einkorn"),9] <- 9
data[which(data$type == "free-threshing tetraploids"),9] <- 6
sub <- data[,c(3,4,9,2)]
colnames(sub) <- c("lat","lng","group","ibs")
# 添加不同形状和颜色的点
my_map <- leaflet() %>%
  addTiles()  %>% 
  addProviderTiles("Esri.WorldStreetMap")%>%
  addCircleMarkers(data=sub,~lng, ~lat, fillColor = ~mypalette(sub$ibs), fillOpacity = 1,
                   color="white", radius=15, stroke=F
  )%>%
  addCircleMarkers(data=point, ~Longitude, ~Latitude,
                   color="black", radius=12
  ) %>%
  addCircleMarkers(data=WE, ~Longitude, ~Latitude,
                   color="green", radius=12
  ) %>%
  addMarkers(
    data = sub,
    lng = ~lng,
    lat = ~lat,
    icon = ~ icons(
      iconUrl = iconFiles[group],
      iconWidth = 35,
      iconHeight = 35,
      iconAnchorX = 0,
      iconAnchorY = 0)
  )%>%
  addLegend( data=sub,pal=mypalette, values=~sub$ibs, opacity=1, title = "ibs distance", position = "bottomright" )

ggplot(data, aes(x=V2, y=V4)) +
  geom_point(size=0.5)+
  geom_line()+
  theme_bw()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("XP-CLR") 


library(tidyverse)
library(hrbrthemes)
sub <- data[which(data$type == "Landrace" | data$type == "Cultivar"),]
ggplot(sub, aes(x=type, y = ABD_0733, fill=type,group=type))+
  geom_boxplot(fill = '#f8766d')+
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic()+  
  theme(
    legend.position="none",
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
  ) 





