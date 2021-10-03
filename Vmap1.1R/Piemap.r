require(reshape)
require (rworldmap)
require(rworldxtra)
#Rht1----
#wheat1 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/wheatLonLat_Rht1.txt", head=T,sep="\t")
wheat1 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/haplo_Rht1.txt", head=T,sep="\t")
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) #对数据进行预处理
wheat_reshape2 <- as.data.frame(wheat_reshape)
#mapPies(wheat_reshape2,xlim=c(-10,120),ylim=c(10,50),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5'),mapRegion='world',symbolSize=1,
#       zColours=c('red','blue','yellow','green','pink'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5','6','7'),symbolSize=1,
        zColours=c('red','blue','yellow','green','pink','orange','black'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")

mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5','6','7','8','9','10'),symbolSize=1,
        zColours=c('red','blue','green','yellow','pink','orange','black','brown','purple','lightblue'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")

mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5'),symbolSize=1,
        zColours=rainbow(5),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")

#Rht-1A----
#wheat1 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/wheatLonLat_Rht1.txt", head=T,sep="\t")
wheat1 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/haplo_Rht-1A.txt", head=T,sep="\t")
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) #对数据进行预处理
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5'),symbolSize=1,
        zColours=rainbow(5),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")

#Rht-1A
#wheat1 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/wheatLonLat_Rht1.txt", head=T,sep="\t")
wheat1 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/haplo_Rht-1B.txt", head=T,sep="\t")
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.5),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.5) #对数据进行预处理
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5'),symbolSize=1,
        zColours=rainbow(5),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")


#VRN-1B
#wheat1 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/wheatLonLat_Rht1.txt", head=T,sep="\t")
wheat1 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/haplo_VRN-1B.txt.txt", head=T,sep="\t")
wheat1 <- wheat1[1:290,]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) #对数据进行预处理
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=rainbow(3),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")

#Migration----
wheat1 <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/13_Plots/08_HaploType/data.txt", header=T,sep="\t")
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$k.5),]
wheat_reshape <- wheat_reshape <- cast(wheat,cluster1_70+cluster2_70~k.5) #对数据进行预处理
wheat_reshape2 <- as.data.frame(wheat_reshape)
#mapPies(wheat_reshape2,xlim=c(-10,120),ylim=c(10,50),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5'),mapRegion='world',symbolSize=1,
#       zColours=c('red','blue','yellow','green','pink'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="cluster2_70",nameY="cluster1_70",nameZs=c('1','2','3','4','5'),symbolSize=1,
        zColours=c('red','blue','yellow','green','pink'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")

mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5','6','7','8','9','10'),symbolSize=1,
        zColours=c('red','blue','green','yellow','pink','orange','black','brown','purple','lightblue'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")

mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5'),symbolSize=1,
        zColours=rainbow(5),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9")
