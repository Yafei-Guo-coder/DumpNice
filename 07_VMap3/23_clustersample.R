#分群
#各地区的landrace: WA, EU, IA, SH, EA, AF, AM
setwd("/Users/guoyafei/Documents/02_Vmap3/16_clustersample")
library(readxl)
library(ggmap)
library(RColorBrewer)
Lulab_germplasm_Info <- read_excel("/Users/guoyafei/Documents/02_Vmap3/16_clustersample/Lulab_germplasm_Info.xlsx")
landrace_WA <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`type(USDA_8)` == "Landrace" & Lulab_germplasm_Info$`Ctnt(9)` == "WA" & as.numeric(Lulab_germplasm_Info$Longitude) < 50 & as.numeric(Lulab_germplasm_Info$Longitude) > 35 & as.numeric(Lulab_germplasm_Info$Latitude) > 35 & as.numeric(Lulab_germplasm_Info$Latitude) < 40),]
landrace_WA$xpclr <- "WA"
dim(landrace_WA)
[1] 48 21
landrace_EU <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`type(USDA_8)` == "Landrace" & Lulab_germplasm_Info$`Ctnt(9)` == "EU" & as.numeric(Lulab_germplasm_Info$Longitude) < 30 & as.numeric(Lulab_germplasm_Info$Latitude) < 50) ,]
landrace_EU$xpclr <- "EU"
dim(landrace_EU)
[1] 44 21
landrace_SH <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`type(USDA_8)` == "Landrace" & Lulab_germplasm_Info$`Ctnt(9)` == "SA" & as.numeric(Lulab_germplasm_Info$Latitude) < 35) ,]
landrace_SH <- landrace_SH[c(6:30,32:47,91:93),]
landrace_SH$xpclr <- "SH"
dim(landrace_SH)
[1] 44 21
landrace_CA_1 <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`type(USDA_8)` == "Landrace" & Lulab_germplasm_Info$`Ctnt(9)` == "WA" & as.numeric(Lulab_germplasm_Info$Longitude) > 55) ,]
landrace_CA_2 <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`type(USDA_8)` == "Landrace" & Lulab_germplasm_Info$`Ctnt(9)` == "CA" & as.numeric(Lulab_germplasm_Info$Latitude) > 35),]
landrace_CA <- rbind(landrace_CA_1,landrace_CA_2)
landrace_IA <- landrace_CA[which(as.numeric(landrace_CA$Latitude) > 25 & as.numeric(landrace_CA$Latitude) < 46) ,]
landrace_IA$xpclr <- "IA"
dim(landrace_IA)
[1] 46 21
landrace_EA <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`type(USDA_8)` == "Landrace" & Lulab_germplasm_Info$`Ctnt(9)` == "EA" & as.numeric(Lulab_germplasm_Info$Longitude) > 98 & as.numeric(Lulab_germplasm_Info$Longitude) < 120 & as.numeric(Lulab_germplasm_Info$Latitude) > 31 & as.numeric(Lulab_germplasm_Info$Latitude) < 41) ,]
landrace_EA$xpclr <- "EA"
dim(landrace_EA)
[1] 53 21
landrace_AM <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`type(USDA_8)` == "Landrace" & Lulab_germplasm_Info$`Ctnt(9)` == "Sth_AM") ,]
landrace_AM$xpclr <- "AM"
dim(landrace_AM)
[1] 20 21
landrace_AF <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`type(USDA_8)` == "Landrace" & Lulab_germplasm_Info$`Ctnt(9)` == "AF") ,]
landrace_AF$xpclr <- "AF"
dim(landrace_AF)
[1] 20 21

mp <- NULL
mapworld <- borders("world",colour = "gray70",fill="white") 
mp <- ggplot() + 
  mapworld + 
  #xlim(0,90) +
  ylim(-60,90) + 
  theme_classic()
color <- brewer.pal(8, "Dark2")[c(1,2,3,4,6)]
color1 <- brewer.pal(8, "Dark2")
color2 <- brewer.pal(8, "Set1")
color3 <- brewer.pal(8, "Set2")

mp+geom_point(aes(x=as.numeric(landrace_WA$Longitude), y=as.numeric(landrace_WA$Latitude),color=landrace_WA$`Ctnt(9)`),size=1)+
  scale_fill_manual(values = c(color1,color2,color3)) 
mp+geom_point(aes(x=as.numeric(landrace_EU$Longitude), y=as.numeric(landrace_EU$Latitude),color=landrace_EU$`Ctnt(9)`),size=1)+
  scale_fill_manual(values = c(color1,color2,color3)) 
mp+geom_point(aes(x=as.numeric(landrace_SA$Longitude), y=as.numeric(landrace_SA$Latitude),color=landrace_SA$`Ctnt(9)`),size=1)+
  scale_fill_manual(values = c(color1,color2,color3)) 
mp+geom_point(aes(x=as.numeric(landrace_CA$Longitude), y=as.numeric(landrace_CA$Latitude),color=landrace_CA$`Ctnt(9)`),size=1)+
  scale_fill_manual(values = c(color1,color2,color3)) 
mp+geom_point(aes(x=as.numeric(landrace_EA$Longitude), y=as.numeric(landrace_EA$Latitude),color=landrace_EA$`Ctnt(9)`),size=1)+
  scale_fill_manual(values = c(color1,color2,color3))
mp+geom_point(aes(x=as.numeric(landrace_AF$Longitude), y=as.numeric(landrace_AF$Latitude),color=landrace_AF$`Ctnt(9)`),size=1)+
  scale_fill_manual(values = c(color1,color2,color3)) 
mp+geom_point(aes(x=as.numeric(landrace_AM$Longitude), y=as.numeric(landrace_AM$Latitude),color=landrace_AM$`Ctnt(9)`),size=1)+
  scale_fill_manual(values = c(color1,color2,color3))

mp+geom_point(aes(x=as.numeric(landrace$Longitude), y=as.numeric(landrace$Latitude),color=landrace$`Ctnt(9)`),size=1)+
  scale_fill_manual(values = c(color1,color2,color3)) 

xpclr <- rbind(landrace_WA,landrace_EU,landrace_IA,landrace_SH,landrace_EA,landrace_AF,landrace_AM)
mp+geom_point(aes(x=as.numeric(xpclr$Longitude), y=as.numeric(xpclr$Latitude),color=xpclr$xpclr),size=1)+
  scale_fill_manual(values = c(color1,color2,color3)) 
write.table(xpclr,"xpclr-landrace.txt",sep="\t",quote=F,row.names = F)

Lulab_germplasm_Info <- read_excel("/Users/guoyafei/Documents/02_VmapIII/16_clustersample/Lulab_germplasm_Info.xlsx")
ABD <- Lulab_germplasm_Info[which((Lulab_germplasm_Info$Genome == "AABBDD" )),]
landrace_1 <- ABD[which((ABD$DataSource == "LuLab")),]
landrace_2 <- ABD[which((ABD$DataSource == "WEGA")),]
landrace_3 <- ABD[which((ABD$DataSource == "XD")),]
landrace_4 <- ABD[which((ABD$DataSource == "JiaoLab")),]
all <- rbind(landrace_1,landrace_2,landrace_3,landrace_4)
mp+geom_point(aes(x=as.numeric(all$Longitude), y=as.numeric(all$Latitude)),color= "orange", size=1)+
  scale_fill_manual(values = "orange") 

#不同年份的cultivar:
setwd("/Users/guoyafei/Documents/02_VmapIII/16_clustersample")
library(readxl)
library(ggmap)
library(RColorBrewer)
Lulab_germplasm_Info <- read_excel("/Users/guoyafei/Documents/02_VmapIII/16_clustersample/Lulab_germplasm_Info.xlsx")
cultivar <- Lulab_germplasm_Info[which(Lulab_germplasm_Info$`Common_name(21)` == "Cultivar" & Lulab_germplasm_Info$Introducedyear != "NA"),]
mp+geom_point(aes(x=as.numeric(cultivar$Longitude), y=as.numeric(cultivar$Latitude),color=as.numeric(cultivar$Introducedyear)),size=1)+
  scale_colour_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 2000)

cultivar_50 <- cultivar[cultivar$Introducedyear <= 1955,]
cultivar_50$year <- "1950"
cultivar_60 <- cultivar[cultivar$Introducedyear <= 1965 & cultivar$Introducedyear > 1955,]
cultivar_60$year <- "1960"
cultivar_70 <- cultivar[cultivar$Introducedyear <= 1975 & cultivar$Introducedyear > 1965,]
cultivar_70$year <- "1970"
cultivar_80 <- cultivar[cultivar$Introducedyear <= 1985 & cultivar$Introducedyear > 1975,]
cultivar_80$year <- "1980"
cultivar_90 <- cultivar[cultivar$Introducedyear <= 1995 & cultivar$Introducedyear > 1985,]
cultivar_90$year <- "1990"
cultivar_00 <- cultivar[cultivar$Introducedyear <= 2005 & cultivar$Introducedyear > 1995,]
cultivar_00$year <- "2000"
cultivar_10 <- cultivar[cultivar$Introducedyear <= 2015 & cultivar$Introducedyear > 2005,]
cultivar_10$year <- "2010"
cultivar_20 <- cultivar[cultivar$Introducedyear <= 2025 & cultivar$Introducedyear > 2015,]
cultivar_20$year <- "2020"

output <- rbind(cultivar_50,cultivar_60,cultivar_70,cultivar_80,cultivar_90,cultivar_00,cultivar_10,cultivar_20)
write.table(output,"cultivar_year.txt",quote=F, row.names = F)

ggplot(cultivar, aes(x = as.numeric(Introducedyear))) +
  geom_histogram(binwidth = 1, fill = "lightblue", colour = "black")+
  theme_classic()


