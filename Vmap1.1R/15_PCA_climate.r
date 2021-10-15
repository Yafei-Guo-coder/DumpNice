#working directory
#/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/NegativeContral
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/NegativeContral")
library(RColorBrewer)
location <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/04_795taxaIBS/795_Location.txt",header=T,stringsAsFactors = F)

Elevation <- read.table("Elevation.txt",header=T,stringsAsFactors = F)
ele <- merge(Elevation,location,by="ID")
colnames(ele)[2] <-"PC1"
Temp <- read.table("temp_PC.txt",header=T,stringsAsFactors = F)
Temp$ID <- rownames(Temp)
tem <- merge(Temp,location,by="ID")
Prec <- read.table("Prec_PC.txt",header=T,stringsAsFactors = F)
Prec$ID <- rownames(Prec)
pre <- merge(Prec,location,by="ID")

bio <- read.table("20bio_PC.txt",header=T,row.names = 1)
bio$ID <- rownames(bio)
bio <- merge(bio,location,by="ID")

data <- list(ele,tem,pre)
pdf("20bio_PC1.pdf",width = 12,height = 7.5)
for(i in c(1:2)){
mp <- NULL
mapworld <- borders("world",colour = "grey90",fill="white") 
mp <- ggplot()+mapworld+ylim(-60,90)+xlim(-25,200)
mp2 <- mp+geom_point(aes(x=max$Logititude,y=max$Latitude),size=2.6)+scale_size(range=c(1,1))+ 
  scale_colour_gradientn(colours = brewer.pal(11, "RdYlBu"))+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  theme(panel.grid = element_blank()) + 
  #theme(panel.grid =element_blank()) +   ## 删去网格线
  #theme(axis.text = element_blank()) +   ## 删去刻度标签
  #theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank())
print(mp2)
}
dev.off()

#看筛选样本的分布位点
#taxa <- read.table("group_V2/max_temp.txt",header=F,stringsAsFactors = F)
#colnames(taxa)[1] <-"ID"
#max <- merge(taxa,location,by="ID")

