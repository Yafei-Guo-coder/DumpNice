#Ppd-1和VRN基因随纬度变化的structure分布
#Working diirectory: 203:yafei:/data1/home/yafei/008_Software/snpEff/Xp-clr_6VIP
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)
require(gridExtra)
library(RColorBrewer)
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/")
annotation_col <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/Taxa_Region_412.txt",header=T,stringsAsFactors = F)
rownames(annotation_col) <- annotation_col[,1]
#load Qmatrix files
#Ppd-1
#sfiles <- list.files(path="/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/Ppd-1", full.names=T)
#slist <- readQ(files=sfiles)
#tabulateQ(qlist=readQ(sfiles))
#summariseQ(tabulateQ(qlist=readQ(sfiles)))
#k <- read.table("Qmatrix.txt",header=T,stringsAsFactors = F,fill = NA)
#Ppd_k2 <- cbind(name,k[,2])
#Ppd_k2$value <- 1
#colnames(Ppd_k2) <- c("ID","Logititude", "Latitude", "qp_gragh", "Xp_clr","Ancestor","value")
#Ppd_k2 = Ppd_k2[!is.na(Ppd_k2$Latitude),]
#Ppd_k2 = Ppd_k2[!is.na(Ppd_k2$Ancestor),]
#Ppd_k2 <- cast(Ppd_k2,Latitude+Logititude~Ancestor) 
#mapPies(Ppd_k2,xlim=c(-120,140),ylim=c(0,40),nameX="Logititude",nameY="Latitude",nameZs=c("1","2","3"),symbolSize=1,
#        zColours=brewer.pal(8, "Set2")[c(2,3,4)],barOrient='vert',oceanCol="white",borderCol = "black", landCol="grey",main="i", lwd = 0.3)

ppd <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/TXT/5k.Ppd-1.11.33947048-33961269.pos.recode.txt",header=F,stringsAsFactors = F)
name <- read.table("Ppd-1.txt",header=F,stringsAsFactors = F)
name <- annotation_col[name[,1],]
ppd2 <- cbind(name,t(ppd[28,]))
ppd3 =ppd2[!is.na(ppd2$Latitude),]
ppd3$R <- NA
ppd3[which(ppd3$Latitude >20 & ppd3$Latitude <=30),7] <- 1
ppd3[which(ppd3$Latitude >30 & ppd3$Latitude <=40),7] <- 2
ppd3[which(ppd3$Latitude >40 & ppd3$Latitude <=50),7] <- 3
ppd3[which(ppd3$Latitude >50),7] <- 4
ppd3[which(ppd3$Latitude == "DD"),7] <- 5
data <- ppd3[,c(6,7)]
colnames(data) <- c("Type", "Region")
aggregate(data$Type, by=list(type=data$Region),mean)

VRN <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/TXT/5k.VRN2-2.24.58272633-58284728.pos.recode.txt",header=F,stringsAsFactors = F)
name <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/VRN-2_250.txt",header=F,stringsAsFactors = F)
name <- annotation_col[name[,1],]
VRN2 <- cbind(name,t(VRN[12,]))
VRN3 =VRN2[!is.na(VRN2$Latitude),]
VRN3$R <- NA
VRN3[which(VRN3$Latitude >20 & VRN3$Latitude <=30),7] <- 1
VRN3[which(VRN3$Latitude >30 & VRN3$Latitude <=40),7] <- 2
VRN3[which(VRN3$Latitude >40 & VRN3$Latitude <=50),7] <- 3
VRN3[which(VRN3$Latitude >50), 7] <- 4
VRN3[which(VRN3$Latitude == "DD"),7] <- 5
data <- VRN3[,c(6,7)]
colnames(data) <- c("Type", "Region")
aggregate(data$Type, by=list(type=data$Region),mean)

ppd_A1 <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/TXT/5k.Ppd-A1.7.36928684-36943202.pos.recode.txt",header=F,stringsAsFactors = F)
name <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/AABBgenome.txt",header=F,stringsAsFactors = F)
name <- annotation_col[name[,1],]
ppd_A12 <- cbind(name,t(ppd_A1[6,]))
ppd_A13 =ppd_A12[!is.na(ppd_A12$Latitude),]
ppd_A13$R <- NA
ppd_A13[which(ppd_A13$Latitude >20 & ppd_A13$Latitude <=30),7] <- 1
ppd_A13[which(ppd_A13$Latitude >30 & ppd_A13$Latitude <=40),7] <- 2
ppd_A13[which(ppd_A13$Latitude >40 & ppd_A13$Latitude <=50),7] <- 3
ppd_A13[which(ppd_A13$Latitude >50), 7] <- 4
ppd_A13[which(ppd_A13$Latitude == "AABB"),7] <- 5
data <- ppd_A13[,c(6,7)]
colnames(data) <- c("Type", "Region")
aggregate(data$Type, by=list(type=data$Region),mean)

Rht_B1 <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/TXT/5k.Rht-B1.21.30856268-30868723.pos.recode.txt",header=F,stringsAsFactors = F)
name <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/AABBgenome.txt",header=F,stringsAsFactors = F)
name <- annotation_col[name[,1],]
Rht_B12 <- cbind(name,t(Rht_B1[6,]))
Rht_B13 =Rht_B12[!is.na(Rht_B12$Latitude),]
Rht_B13$R <- NA
Rht_B13[which(Rht_B13$Latitude >20 & Rht_B13$Latitude <=30),7] <- 1
Rht_B13[which(Rht_B13$Latitude >30 & Rht_B13$Latitude <=40),7] <- 2
Rht_B13[which(Rht_B13$Latitude >40 & Rht_B13$Latitude <=50),7] <- 3
Rht_B13[which(Rht_B13$Latitude >50), 7] <- 4
Rht_B13[which(Rht_B13$Latitude == "AABB"),7] <- 5
data <- Rht_B13[,c(6,7)]
colnames(data) <- c("Type", "Region")
aggregate(data$Type, by=list(type=data$Region),mean)

Sr45 <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/TXT/5k.Sr45.5.19336296-19351065.pos.recode.txt",header=F,stringsAsFactors = F)
name <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/DDgenome.txt",header=F,stringsAsFactors = F)
name <- annotation_col[name[,1],]
Sr452 <- cbind(name,t(Sr45[6,]))
Sr453 =Sr452[!is.na(Sr452$Latitude),]
Sr453$R <- NA
Sr453[which(Sr453$Latitude >20 & Sr453$Latitude <=30),7] <- 1
Sr453[which(Sr453$Latitude >30 & Sr453$Latitude <=40),7] <- 2
Sr453[which(Sr453$Latitude >40 & Sr453$Latitude <=50),7] <- 3
Sr453[which(Sr453$Latitude >50), 7] <- 4
Sr453[which(Sr453$Latitude == "DD"),7] <- 5
data <- Sr453[,c(6,7)]
colnames(data) <- c("Type", "Region")
aggregate(data$Type, by=list(type=data$Region),mean)


Sr33 <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/TXT/5k.Sr33.5.11446423-11464353.pos.recode.txt",header=F,stringsAsFactors = F)
name <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/DDgenome.txt",header=F,stringsAsFactors = F)
name <- annotation_col[name[,1],]
Sr332 <- cbind(name,t(Sr33[6,]))
Sr333 =Sr332[!is.na(Sr332$Latitude),]
Sr333$R <- NA
Sr333[which(Sr333$Latitude >20 & Sr333$Latitude <=30),7] <- 1
Sr333[which(Sr333$Latitude >30 & Sr333$Latitude <=40),7] <- 2
Sr333[which(Sr333$Latitude >40 & Sr333$Latitude <=50),7] <- 3
Sr333[which(Sr333$Latitude >50), 7] <- 4
Sr333[which(Sr333$Latitude == "DD"),7] <- 5
data <- Sr333[,c(6,7)]
colnames(data) <- c("Type", "Region")
aggregate(data$Type, by=list(type=data$Region),mean)