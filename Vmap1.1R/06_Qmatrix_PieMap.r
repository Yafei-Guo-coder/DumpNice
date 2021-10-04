library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
#map
require("RColorBrewer")
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/02_Structure_map/Select")
data <- read.table("select_taxa.txt",header=T,sep="\t",stringsAsFactors = F)
dataA <- data[,c(1,7,8,14,15)]
colnames(dataA)[1:5] <- c("type","Latitude", "Longitude", "Sub.population.6", "value")
wheat2 = dataA[!is.na(dataA$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.6),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.6) 
wheat_reshape2 <- as.data.frame(wheat_reshape)

#8
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4","5","6","7","8"),symbolSize=1,
        zColours=brewer.pal(8, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#6
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4","5","6"),symbolSize=1,
        zColours=brewer.pal(6, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#5
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4","5"),symbolSize=1,
        zColours=brewer.pal(5, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#4
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,36),nameX="Longitude",nameY="Latitude",nameZs=c('1','2',"3","4"),symbolSize=1,
        zColours=brewer.pal(4, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#3
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=brewer.pal(3, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")
#2
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2'),symbolSize=1,
        zColours=brewer.pal(2, "Set2"),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="A subgenome")

par(mfcol=c(5,1))