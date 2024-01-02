library(readxl)
library(ggmap)
library(RColorBrewer)
a <- read_excel("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx", sheet=1, na='NA')
taxa <- a[,c(1,7,9,11,12)]

#cultivar
taxa <- taxa[which(taxa$`type(USDA_8)`=="Cultivar"),]
#landrace
taxa <- taxa[which(taxa$`type(USDA_8)`=="Landrace"),]
#otherHexaploids
taxa <- taxa[which(taxa$`type(USDA_8)`=="OtherHexaploids"),]
#otherTetraploids
taxa <- taxa[which(taxa$`type(USDA_8)`=="OtherTetraploids"),]
#free_threshing_tetraploids
taxa <- taxa[which(taxa$`type(USDA_8)`=="Free_threshing_tetraploids"),]
#strangulata
taxa <- taxa[which(taxa$`type(USDA_8)`=="Strangulata"),]
#wild_emmer
taxa <- taxa[which(taxa$`type(USDA_8)`=="Wild_emmer"),]
#Domesticated_emmer
taxa <- taxa[which(taxa$`type(USDA_8)`=="Domesticated_emmer"),]

#绘图:全部的
mp <- NULL
mapworld <- borders("world",colour = "gray90",fill="gray90") 
mp <- ggplot() + 
  mapworld + 
  #xlim(0,90) +
  ylim(-60,90) + 
  theme_classic()

color <- brewer.pal(8, "Dark2")[c(1,2,3,4,6)]
color <- brewer.pal(8, "Dark2")[c(6)]

#AABBDD
mp+geom_point(aes(x=taxa$Logitude, y=taxa$Latitude, color=taxa$Genome),size=1)+
  scale_size(range=c(1,1)) + 
  scale_color_manual(values = color) +
  #theme_classic() +
  #theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "null",
        panel.border = element_blank())

mp+geom_point(aes(x=taxa$Logitude, y=taxa$Latitude,color = color), size=1.5)+
  scale_size(range=c(1,1)) + 
  scale_color_manual(values = color) +
  #theme_classic() +
  #theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "null",
        panel.border = element_blank())  

mp+geom_point(aes(x=p$Longitude, y=p$Latitude,color = p$outcome), size=1)+
  scale_size(range=c(1,1)) + 
  scale_color_manual(values = c("#215085","#D53B3D")) 
