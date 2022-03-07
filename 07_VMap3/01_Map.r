library(readxl)
library(ggmap)
library(RColorBrewer)
a <- read_excel("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx", sheet=1, na='NA')
taxa <- a[,c(4,7,8,10,14)]

#cultivar
taxa <- taxa[which(taxa$`type(8)`=="Cultivar"),]
#landrace
taxa <- taxa[which(taxa$`type(8)`=="Landrace"),]
#otherHexaploids
taxa <- taxa[which(taxa$`type(8)`=="OtherHexaploids"),]
#otherTetraploids
taxa <- taxa[which(taxa$`type(8)`=="OtherTetraploids"),]
#free_threshing_tetraploids
taxa <- taxa[which(taxa$`type(8)`=="Free_threshing_tetraploids"),]
#strangulata
taxa <- taxa[which(taxa$`type(8)`=="Strangulata"),]
#wild_emmer
taxa <- taxa[which(taxa$`type(8)`=="Wild_emmer"),]

#绘图:全部的
mp <- NULL
mapworld <- borders("world",colour = "gray70",fill="white") 
mp <- ggplot() + 
  mapworld + 
  #xlim(0,90) +
  ylim(-60,90) + 
  theme_classic()
color <- brewer.pal(8, "Dark2")[c(1,2,3,4,6)]
color <- brewer.pal(8, "Dark2")[c(7)]
#AABBDD
mp+geom_point(aes(x=taxa$Logitude, y=taxa$Latitude, color=taxa$Genome),size=0.5)+
  scale_size(range=c(1,1))+ 
  scale_color_manual(values = color) +
  #theme_classic()+
  #theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "null",
        panel.border = element_blank()) 
  
#Check Sample
library(dplyr)
a <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global.xlsx", sheet=1, na='NA')
b <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (1).xlsx", sheet=1, na='NA')  
c <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (2).xlsx", sheet=1, na='NA')  
d <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (3).xlsx", sheet=1, na='NA')  
e <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (4).xlsx", sheet=1, na='NA')  
f <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (5).xlsx", sheet=1, na='NA')  
g <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (6).xlsx", sheet=1, na='NA')  
h <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (7).xlsx", sheet=1, na='NA')  
i <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (8).xlsx", sheet=1, na='NA')  
j <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (9).xlsx", sheet=1, na='NA')  
k <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (10).xlsx", sheet=1, na='NA')  
l <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (11).xlsx", sheet=1, na='NA')  
m <- read_excel("/Users/guoyafei/Downloads/Search Accessions GRIN-Global (12).xlsx", sheet=1, na='NA') 

all <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m)
all2 <- distinct(all)
write.table(all2, "USDA_sample.txt", quote=F, sep="\t", row.names = F)
