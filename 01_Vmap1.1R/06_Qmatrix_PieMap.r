#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5")
annotation_col <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Taxa_Region_225.txt",header=T,stringsAsFactors = F)
rownames(annotation_col) = c(1:325)
seq <- annotation_col[,1]
#cluster samples: 按照算Xp-clr的区域划分
out <- annotation_col[,c(4,5)]
lat_mean <- tapply(out[,1],out$Region_sub2,mean,na.rm = TRUE)
lon_mean <- tapply(out[,2],out$Region_sub2,mean,na.rm = TRUE)
out$cluster1 <- NA
out$cluster2 <- NA
for(i in 1:225) {
        for(j in names(lat_mean)){
                if(out[i,3] == j ){clas
                        out[i,4] <- as.numeric(lat_mean[j])
                        out[i,5] <- as.numeric(lon_mean[j])
                } 
        }
}
mode <- as.data.frame(cbind(annotation_col[,2],as.numeric(out[,4]),as.numeric(out[,5])),stringsAsFactors = F)
mode$V2 <- as.numeric(mode$V2)
mode$V3 <- as.numeric(mode$V3)
#map
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/TXT" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
        paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
        read.table(x, header=F,stringsAsFactors = F)})
pdf("cut_8_piemap.pdf")
for (i in c(1:82)){
        all <- as.matrix(data[[i]])
        colnames(all) <- c(1:361)
        all <- all[,seq]
        data.e <- dist(t(all))
        model <- hclust(data.e,method = "complete")
        #plot(model)
        result <- as.numeric(cutree(model, k=8))
        out2 <- as.data.frame(cbind(mode,result))
        out2$value <- 1
        colnames(out2)[1:4] <- c("Latitude", "Logititude", "Sub.population.2", "value")
        wheat2 = out2[!is.na(out2$Latitude),]
        wheat = wheat2[!is.na(wheat2$Sub.population.2),]
        wheat_reshape <- cast(test,Latitude+Logititude~Sub.population.2) 
        wheat_reshape2 <- as.data.frame(wheat_reshape)
        mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(0,40),nameX="Logititude",nameY="Latitude",nameZs=c("1","2","3"),symbolSize=1,
                zColours=brewer.pal(8, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main=i)
}
dev.off()