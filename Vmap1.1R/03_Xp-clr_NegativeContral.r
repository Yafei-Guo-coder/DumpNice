#XP-CLR negative contral
library(RColorBrewer)
library(ggplot2)
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/NegativeContral")
dist <- read.table("dis_matrix.txt",header=T,stringsAsFactors = F)
out <- strsplit(sub('_',':', dist[,1]) , ":")

dist$ID1 <- NA
dist$ID2 <- NA

a<-1
for (i in c(1:length(out))){
  dist[a,6] <- out[[i]][1]
  dist[a,7] <- out[[i]][2]
  a <- a+1
}

WA <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/NegativeContral/group/WA.txt",header =F,stringsAsFactors = F)
North1 <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/NegativeContral/group/North1.txt",header=F,stringsAsFactors = F)
North2 <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/NegativeContral/group/North2.txt",header=F,stringsAsFactors = F)
#precipitation
pdf("Precipitation_dist.pdf",height = 10,width = 10)
#WA VS North1
dist$Type1 <- NA
dist[which((dist$ID1 %in% WA[,1] & dist$ID2 %in% North1[,1])) ,8] <- "Yes"
dist[which((dist$ID1 %in% North1[,1] & dist$ID2 %in% WA[,1])) ,8] <- "Yes"

p <- ggplot(dist, aes(geopolitical_dis,PREC_dis, color=Type1)) +
  geom_point( size=3,alpha=0.2)+
  theme_classic()+
  ggtitle("WA VS North1")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

print(p)
#WA VS North2
dist$Type1 <- NA
dist[which((dist$ID1 %in% WA[,1] & dist$ID2 %in% North2[,1])) ,8] <- "Yes"
dist[which((dist$ID1 %in% North2[,1] & dist$ID2 %in% WA[,1])) ,8] <- "Yes"

p <- ggplot(dist, aes(geopolitical_dis,PREC_dis, color=Type1)) +
  geom_point( size=3,alpha=0.2)+
  theme_classic()+
  ggtitle("WA VS North2")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

print(p)

#North1 VS North2
dist$Type1 <- NA
dist[which((dist$ID1 %in% North2[,1] & dist$ID2 %in% North1[,1])) ,8] <- "Yes"
dist[which((dist$ID1 %in% North1[,1] & dist$ID2 %in% North2[,1])) ,8] <- "Yes"

p <- ggplot(dist, aes(geopolitical_dis,PREC_dis, color=Type1)) +
  geom_point( size=3,alpha=0.2)+
  theme_classic()+
  ggtitle("North1 VS North2")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
print(p)
dev.off()
