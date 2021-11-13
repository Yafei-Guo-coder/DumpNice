#Ppd-1和VRN基因随纬度变化的structure分布
#Working diirectory: 203:yafei:/data1/home/yafei/008_Software/snpEff/Xp-clr_6VIP
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)
require(gridExtra)
library(RColorBrewer)
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure")
annotation_col <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Taxa_Region_225.txt",header=T,stringsAsFactors = F)
rownames(annotation_col) <- annotation_col[,1]
#load Qmatrix files
#Ppd-1
sfiles <- list.files(path="/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/Ppd-1", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))

name <- read.table("Ppd-1.txt",header=F,stringsAsFactors = F)
name <- annotation_col[name[,1],]
k <- read.table("Qmatrix.txt",header=T,stringsAsFactors = F,fill = NA)
for(i in c(1:length(slist))){
  rownames(slist[[1]]) <- name[,1]
}

pdf("cut_8_piemap.pdf")
for (i in c(1:82)){
  all <- as.matrix(data[[i]])
  colnames(all) <- c(1:225)
  all <- all[,seq]
  data.e <- dist(t(all))
  model <- hclust(data.e,method = "complete")
  #plot(model)
  result <- as.numeric(cutree(model, k=8))
  out2 <- as.data.frame(cbind(mode,result))
  out2$value <- 1
  colnames(out2)[1:5] <- c("type","Latitude", "Logititude", "Sub.population.2", "value")
  wheat2 = out2[!is.na(out2$Latitude),]
  wheat = wheat2[!is.na(wheat2$Sub.population.2),]
  wheat_reshape <- cast(wheat,Latitude+Logititude~Sub.population.2) 
  wheat_reshape2 <- as.data.frame(wheat_reshape)
  mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(0,40),nameX="Logititude",nameY="Latitude",nameZs=c("1","2","3","4","5","6","7","8"),symbolSize=1,
          zColours=brewer.pal(8, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main=i)
}
dev.off()









