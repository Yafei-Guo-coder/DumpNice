library(readxl)
library(ggmap)
library(RColorBrewer)
a <- read_excel("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx", sheet=1, na='NA')
taxa <- as.data.frame(a[which(a$Latitude != "-" & a$Longitude != "-" ),c(10,14,17,18)])

taxa_jm <- as.data.frame(a[which(a$Latitude != "-" & a$Longitude != "-" & a$DataSource != "XD" & a$Genome == "AABBDD"),c(14,17,18)])

taxa$Latitude <- as.numeric(taxa$Latitude)
taxa$Longitude <- as.numeric(taxa$Longitude)

taxa_jm$Latitude <- as.numeric(taxa_jm$Latitude)
taxa_jm$Longitude <- as.numeric(taxa_jm$Longitude)

#绘图:全部的
mp <- NULL
mapworld <- borders("world",colour = "gray90",fill="gray90") 
mp <- ggplot() + 
  mapworld + 
  #xlim(0,90) +
  ylim(-60,90) + 
  theme_classic()

color <- brewer.pal(8, "Dark2")[c(4,6,1,2,3,5)]
color <- brewer.pal(8, "Dark2")[c(6)]

mp+geom_point(aes(x=taxa$Longitude, y=taxa$Latitude,color = taxa$Genome), size=1.5)+
  scale_size(range=c(1,1))+
  scale_color_manual(values = color) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "null",
        panel.border = element_blank())  

#AABBDD
sub <- taxa[which(taxa$type_final == "Wild_emmer"),]
sub <- taxa[which(taxa$type_final == "Domesticated_emmer"),]
sub <- taxa[which(taxa$type_final == "Free_threshing_tetraploids"),]
sub <- taxa[which(taxa$type_final == "Landrace"),]
sub <- taxa[which(taxa$type_final == "Cultivar"),]
sub <- taxa[which(taxa$type_final == "Strangulata"),]

mp+geom_point(aes(x=sub$Longitude, y=sub$Latitude), color=color[6],size=1.5)+
  scale_size(range=c(1,1))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "null",
        panel.border = element_blank())

mp+geom_point(aes(x=taxa_jm$Longitude, y=taxa_jm$Latitude, color = taxa_jm$type_final),size=1.5,alpha=0.8,stroke = 0.08)+
  scale_size(range=c(1,1))+
  scale_color_manual(values = c("#E69F00","#56B4E9","#A07D35")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank())




