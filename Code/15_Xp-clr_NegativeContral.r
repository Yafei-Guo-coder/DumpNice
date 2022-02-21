#XP-CLR negative contral
library(RColorBrewer)
library(ggplot2)
setwd("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment/02_XP-CLR/NegativeContral")
dist <- read.table("dis_matrix.txt",header=T,stringsAsFactors = F)
out <- strsplit(dist[,1], ":")
dist$ID1 <- NA
dist$ID2 <- NA
a<-1
for (i in c(1:length(out))){
  dist[a,6] <- out[[i]][1]
  dist[a,7] <- out[[i]][2]
  a <- a+1
}
#dim(dist[which(dist$geopolitical_dis>100 & dist$TEMP_dis<10),])
WA <- read.table("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment/02_XP-CLR/NegativeContral/group_V2/WA.txt",header =F,stringsAsFactors = F)
North1 <- read.table("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment/02_XP-CLR/NegativeContral/group_V2/North1.txt",header=F,stringsAsFactors = F)
North2 <- read.table("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment/02_XP-CLR/NegativeContral/group_V2/North2.txt",header=F,stringsAsFactors = F)
#Temperature
pdf("Temperature_dist.pdf",height = 10,width = 10)
#WA VS North1
dist$Type1 <- NA
dist[which((dist$ID1 %in% WA[,1] & dist$ID2 %in% North1[,1])) ,8] <- "Yes"
dist[which((dist$ID1 %in% North1[,1] & dist$ID2 %in% WA[,1])) ,8] <- "Yes"
p <- ggplot(dist, aes(geopolitical_dis,TEMP_dis, color=Type1)) +
  geom_point( size=3,alpha=0.2)+
  theme_classic()+
  ggtitle("WA VS North1")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
print(p)
#WA VS North2
dist$Type1 <- NA
dist[which((dist$ID1 %in% WA[,1] & dist$ID2 %in% North2[,1])) ,8] <- "Yes"
dist[which((dist$ID1 %in% North2[,1] & dist$ID2 %in% WA[,1])) ,8] <- "Yes"

p <- ggplot(dist, aes(geopolitical_dis,TEMP_dis, color=Type1)) +
  geom_point( size=3,alpha=0.2)+
  theme_classic()+
  ggtitle("WA VS North2")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

print(p)

#North1 VS North2
dist$Type1 <- NA
dist[which((dist$ID1 %in% North2[,1] & dist$ID2 %in% North1[,1])),8] <- "Yes"
dist[which((dist$ID1 %in% North1[,1] & dist$ID2 %in% North2[,1])),8] <- "Yes"

p <- ggplot(dist, aes(geopolitical_dis,TEMP_dis, color=Type1)) +
  geom_point( size=3,alpha=0.2)+
  theme_classic()+
  ggtitle("North1 VS North2")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
print(p)
dev.off()
