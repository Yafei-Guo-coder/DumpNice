#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require(rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)

setwd("/Users/guoyafei/Desktop/群体结构")
#load Qmatrix files
sfiles <- list.files(path="/Users/guoyafei/Desktop/群体结构/Q矩阵", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))

pdf("Q.pdf",width = 24,height = 8)
p1 <- plotQ(slist[2:6],returnplot=T,exportplot=F,basesize=11,
            sortind="all",showindlab=F,showyaxis=T,showticks=T,sharedindlab=T,clustercol=brewer.pal(7, "Set2"))
grid.arrange(p1$plot[[1]],p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],nrow=5)
dev.off()

#画PC散点图（随机位点和适应性位点）----
library(scatterpie)
library(RColorBrewer)
color <- brewer.pal(9, "Set1")
setwd("/Users/guoyafei/Desktop/群体结构")
for (i in c(3:7)) {
  infile <- paste("random.Q",i,".txt",sep="")
  outfile <- paste("random.Q",i,".pdf",sep="")
  data <- read.table(infile, header=T, stringsAsFactors = F)
  data <- as.data.frame(data)
  cols <- paste("Q",c(1:i),sep="")
  pdf(outfile,height=4,width=5)
  p <- ggplot() + 
    geom_scatterpie(aes(x=PC1, y=PC2), data=data,cols = cols, color=NA) + 
    coord_equal()+
    scale_fill_brewer(palette = "Set2")+
    theme_classic()
  print(p)
  dev.off()
}
#PC1:0.39 PC2:0.23
for (i in c(3:7)) {
  infile <- paste("adap.Q",i,".txt",sep="")
  outfile <- paste("adap.Q",i,".pdf",sep="")
  data <- read.table(infile, header=T, stringsAsFactors = F)
  data <- as.data.frame(data)
  cols <- paste("Q",c(1:i),sep="")
  pdf(outfile,height=4,width=5)
  p <- ggplot() + 
    geom_scatterpie(aes(x=PC1, y=PC2), data=data,cols = cols, color=NA) + 
    coord_equal()+
    scale_fill_brewer(palette = "Set2")+
    theme_classic()
  print(p)
  dev.off()
}
#PC1:0.527 PC2:0.16

#adapPC_fst_V3----
library(scatterpie)
library(RColorBrewer)
color <- brewer.pal(9, "Set1")
setwd("/Users/guoyafei/Desktop/群体结构/adapPC_fst_V3")
for (i in c("all_struct.Q5","prec_struct.Q5","soil_struct.Q5","solar_struct.Q5","temp_struct.Q5")) {
  infile <- paste(i,".txt",sep="")
  outfile <- paste(i,".pdf",sep="")
  data <- read.table(infile, header=T, stringsAsFactors = F)
  data <- as.data.frame(data)
  cols <- paste("Q",c(1:5),sep="")
  pdf(outfile,height=4,width=5)
  p <- ggplot() + 
    geom_scatterpie(aes(x=PC1, y=PC2), data=data,cols = cols, color=NA) + 
    coord_equal()+
    scale_fill_brewer(palette = "Set2")+
    theme_classic()
  print(p)
  dev.off()
}
#PC1: PC2:
i=5
sub <- data[which(data$PC1 > -0.025 & data$PC1 < 0.01),]
sub2 <- sub[which(sub$PC2 > -0.01 & sub$PC2 < 0.03),]
pdf("放大.pdf",height=4,width=5)
p <- ggplot() + 
  geom_scatterpie(aes(x=PC1, y=PC2), data=sub2, cols = cols, color=NA) + 
  coord_equal()+
  xlim(-0.025,0.01)+
  ylim(-0.01,0.03)+
  scale_fill_brewer(palette = "Set2")+
  theme_classic()
print(p)
dev.off()

#adap_soil7----
setwd("/Users/guoyafei/Desktop/群体结构/adap_soil7")
for (i in c(3:7)) {
  infile <- paste("adap.Q",i,".txt",sep="")
  outfile <- paste("adap.Q",i,".pdf",sep="")
  data <- read.table(infile, header=T, stringsAsFactors = F)
  data <- as.data.frame(data)
  cols <- paste("Q",c(1:i),sep="")
  pdf(outfile,height=4,width=5)
  p <- ggplot() + 
    geom_scatterpie(aes(x=PC1, y=PC2), data=sub2, cols = cols, color=NA) + 
    coord_equal()+
    xlim(-0.025,0.01)+
    ylim(-0.01,0.03)+
    scale_fill_brewer(palette = "Set2")+
    theme_classic()
  print(p)
  dev.off()
}
#PC1:0.26 PC2:0.25

#画所有样本地图上的25哥群体的群体结构组成----
#准备输入文件
setwd("/Users/guoyafei/Desktop/群体结构")
library(tidyr)
library(reshape2)
library ("geosphere")
geo <- read.table("/Users/guoyafei/Desktop/baypass/471_baypass_taxa.Info", header=T,stringsAsFactors = F)
data <- read.table("random.Q5.txt", header=T, stringsAsFactors = F)
data <- as.data.frame(data)
m <- merge(data,geo, by.x="IID", by.y="ID")
m$r <- 1
sub <- m[,c(7,8,9,10,11,15,52)]
sub_geo <- m[,c(12,13,15)]
M1 <- aggregate(sub, by=list(sub$type),FUN=sum)
M2 <- aggregate(sub_geo, by=list(sub_geo$type),FUN=mean)
all <- cbind(M1[,c(2:6,8)],M2[,c(2:4)])
write.table(all,"random.map.txt",sep="\t",row.names = F,quote=F)

data <- read.table("adap_soil7/adap.Q5.txt", header=T, stringsAsFactors = F)
data <- as.data.frame(data)
m <- merge(data,geo, by.x="IID", by.y="ID")
m$r <- 1
sub <- m[,c(7,8,9,10,11,15,52)]
sub_geo <- m[,c(12,13,15)]
M1 <- aggregate(sub, by=list(sub$type),FUN=sum)
M2 <- aggregate(sub_geo, by=list(sub_geo$type),FUN=mean)
all <- cbind(M1[,c(2:6,8)],M2[,c(2:4)])
write.table(all,"adapt_soil7.map.txt",sep="\t",row.names = F,quote=F)

#画图
setwd("/Users/guoyafei/Desktop/群体结构/")
library(scatterpie)
#R4.2
library(RColorBrewer)
color <- brewer.pal(9, "Set1")
data <- read.table("random.map.txt", header=T, stringsAsFactors = F)
cols <- paste("Q",c(1:5),sep="")
world <- map_data('world')
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="grey") +
  coord_quickmap()
pdf("random.map.pdf",height=5,width=10)
p + geom_scatterpie(aes(x=Longitude, y=Latitude, r=log(r)),
                    data=data, cols=cols, color=NA,alpha=.9)
#  geom_scatterpie_legend(log(data$r), x=-160, y=-55)
print(p)
dev.off()


