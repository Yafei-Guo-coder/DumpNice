library(readxl)
library(ggmap)
library(RColorBrewer)
a <- read_excel("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx", sheet=1, na='NA')
taxa <- a[,c(4,7,8,10,14)]

#绘图
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="white") 
mp<-ggplot()+mapworld+ylim(-60,90) +theme_classic()
color <- brewer.pal(8, "Dark2")[c(1,2,3,4,6)]
#AABBDD
mp2<-mp+geom_point(aes(x=taxa$Logitude, y=taxa$Latitude, color=taxa$Genome),size=0.5)+
  scale_size(range=c(1,1))+ 
  scale_color_manual(values = color) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) 
  theme(axis.ticks = element_blank()) 
