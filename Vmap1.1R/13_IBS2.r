#Server Working directory
#204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6/VCF
library(ggplot2)
library(RColorBrewer)
library(ggmap)
library(sp)
library(maptools)
library(maps)
library(psych)
library(reshape)
library("corrplot")

setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/04_795taxaIBS/V2_noAFAM/")
path <- "/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/04_795taxaIBS/V2_noAFAM/group" 
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
group <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})

#> names(group)
#[1] "sub_225Landrace.txt"        "sub_Domesticated_emmer.txt" "sub_Macha.txt"             
#[4] "sub_Persian_wheat.txt"      "sub_Polish_wheat.txt"       "sub_Rivet_wheat.txt"       
#[7] "sub_Spelt.txt"              "sub_Xinjiang_wheat.txt"   

path <- "/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/04_795taxaIBS/ibs" 
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
ibs <- lapply(filePath, function(x){
  read.table(x, row.names = 1, header=T,stringsAsFactors = F)})

#> names(ibs)
#[1] "Alineage.ibs.txt" "Blineage.ibs.txt" "Dlineage.ibs.txt"

location <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/04_795taxaIBS/795_Location.txt",header=T,stringsAsFactors = F)

Spelt <- data.frame(ID=group[[7]][,1],IBS="NA",Type="Spelt")
Macha <- data.frame(ID=group[[3]][,1],IBS="NA",Type="Macha")
Xinjiang <- data.frame(ID=group[[8]][,1],IBS="NA",Type="Xinjiang")
Persian <- data.frame(ID=group[[4]][,1],IBS="NA",Type="Persian")

#1. 提取IBS(行是祖先，列是杂交种)
  #Spelt
    Landrace_Spelt_A <- ibs[[1]][group[[1]][,1],group[[7]][,1]]
    Landrace_Spelt_B <- ibs[[2]][group[[1]][,1],group[[7]][,1]]
    Landrace_Spelt_AB <- cbind(Landrace_Spelt_A,Landrace_Spelt_B)
    Landrace_Spelt_AB_mean <- data.frame(ID=names(apply(Landrace_Spelt_AB,1,mean)),IBS=as.numeric(apply(Landrace_Spelt_AB,1,mean)))
    Landrace_Spelt_AB_mean$Type = "Landrace_Spelt_AB"
    A <- rbind(Landrace_Spelt_AB_mean,Spelt)
    A_addLoc <- merge(A,location,by="ID")[,2:5]
    A_addLoc$IBS <- as.numeric(A_addLoc$IBS)
    melt <- melt(A_addLoc,id=c("Type","Latitude","Logititude"))
    A_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    A_od <- A_cast[order(A_cast$IBS),]
    A_od[1:5,1] <- "Min"
    
    DomeEmmer_Spelt_A <- ibs[[1]][group[[2]][,1],group[[7]][,1]]
    DomeEmmer_Spelt_B <- ibs[[2]][group[[2]][,1],group[[7]][,1]]
    DomeEmmer_Spelt_AB <- cbind(DomeEmmer_Spelt_A,DomeEmmer_Spelt_B)
    DomeEmmer_Spelt_AB_mean <- data.frame(ID=names(apply(DomeEmmer_Spelt_AB,1,mean)),IBS=as.numeric(apply(DomeEmmer_Spelt_AB,1,mean)))
    DomeEmmer_Spelt_AB_mean$Type = "DomeEmmer_Spelt_AB"
    B <- rbind(DomeEmmer_Spelt_AB_mean,Spelt)
    B_addLoc <- merge(B,location,by="ID")[,2:5]
    B_addLoc$IBS <- as.numeric(B_addLoc$IBS)
    melt <- melt(B_addLoc,id=c("Type","Latitude","Logititude"))
    B_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    B_od <- B_cast[order(B_cast$IBS),]
    B_od[1:5,1] <- "Min"

    Landrace_Spelt_D <- ibs[[3]][group[[1]][,1],group[[7]][,1]]
    Landrace_Spelt_D_mean <- data.frame(ID=names(apply(Landrace_Spelt_D,1,mean)),IBS=as.numeric(apply(Landrace_Spelt_D,1,mean)))
    Landrace_Spelt_D_mean$Type = "Landrace_Spelt_D"
    C <- rbind(Landrace_Spelt_D_mean,Spelt)
    C_addLoc <- merge(C,location,by="ID")[,2:5]
    C_addLoc$IBS <- as.numeric(C_addLoc$IBS)
    melt <- melt(C_addLoc,id=c("Type","Latitude","Logititude"))
    C_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    C_od <- C_cast[order(C_cast$IBS),]
    C_od[1:5,1] <- "Min"

  #Macha
    Landrace_Macha_A <- ibs[[1]][group[[1]][,1],group[[3]][,1]]
    Landrace_Macha_B <- ibs[[2]][group[[1]][,1],group[[3]][,1]]
    Landrace_Macha_AB <- cbind(Landrace_Macha_A,Landrace_Macha_B)
    Landrace_Macha_AB_mean <- data.frame(ID=names(apply(Landrace_Macha_AB,1,mean)),IBS=as.numeric(apply(Landrace_Macha_AB,1,mean)))
    Landrace_Macha_AB_mean$Type = "Landrace_Macha_AB"
    D <- rbind(Landrace_Macha_AB_mean,Macha)
    D_addLoc <- merge(D,location,by="ID")[,2:5]
    D_addLoc$IBS <- as.numeric(D_addLoc$IBS)
    melt <- melt(D_addLoc,id=c("Type","Latitude","Logititude"))
    D_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    D_od <- D_cast[order(D_cast$IBS),]
    D_od[1:5,1] <- "Min"
    
    Landrace_Macha_D <- ibs[[3]][group[[1]][,1],group[[3]][,1]]
    Landrace_Macha_D_mean <- data.frame(ID=names(apply(Landrace_Macha_D,1,mean)),IBS=as.numeric(apply(Landrace_Macha_D,1,mean)))
    Landrace_Macha_D_mean$Type = "Landrace_Macha_D"
    E <- rbind(Landrace_Macha_D_mean,Macha)
    E_addLoc <- merge(E,location,by="ID")[,2:5]
    E_addLoc$IBS <- as.numeric(E_addLoc$IBS)
    melt <- melt(E_addLoc,id=c("Type","Latitude","Logititude"))
    E_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    E_od <- E_cast[order(E_cast$IBS),]
    E_od[1:5,1] <- "Min"
    
    DomeEmmer_Macha_A <- ibs[[1]][group[[2]][,1],group[[3]][,1]]
    DomeEmmer_Macha_B <- ibs[[2]][group[[2]][,1],group[[3]][,1]]
    DomeEmmer_Macha_AB <- cbind(DomeEmmer_Macha_A,DomeEmmer_Macha_B)
    DomeEmmer_Macha_AB_mean <- data.frame(ID=names(apply(DomeEmmer_Macha_AB,1,mean)),IBS=as.numeric(apply(DomeEmmer_Macha_AB,1,mean)))
    DomeEmmer_Macha_AB_mean$Type = "DomeEmmer_Macha_AB"
    F <- rbind(DomeEmmer_Macha_AB_mean,Macha)
    F_addLoc <- merge(F,location,by="ID")[,2:5]
    F_addLoc$IBS <- as.numeric(F_addLoc$IBS)
    melt <- melt(F_addLoc,id=c("Type","Latitude","Logititude"))
    F_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    F_od <- F_cast[order(F_cast$IBS),]
    F_od[1:5,1] <- "Min"
    
  #Xinjiang
    Landrace_Xinjiang_A <- ibs[[1]][group[[1]][,1],group[[8]][,1]]
    Landrace_Xinjiang_B <- ibs[[2]][group[[1]][,1],group[[8]][,1]]
    Landrace_Xinjiang_AB <- cbind(Landrace_Xinjiang_A,Landrace_Xinjiang_B)
    Landrace_Xinjiang_AB_mean <- data.frame(ID=names(apply(Landrace_Xinjiang_AB,1,mean)),IBS=as.numeric(apply(Landrace_Xinjiang_AB,1,mean)))
    Landrace_Xinjiang_AB_mean$Type = "Landrace_Xinjiang_AB"
    G <- rbind(Landrace_Xinjiang_AB_mean,Xinjiang)
    G_addLoc <- merge(G,location,by="ID")[,2:5]
    G_addLoc$IBS <- as.numeric(G_addLoc$IBS)
    melt <- melt(G_addLoc,id=c("Type","Latitude","Logititude"))
    G_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    G_od <- G_cast[order(G_cast$IBS),]
    G_od[1:5,1] <- "Min"
    
    Landrace_Xinjiang_D <- ibs[[3]][group[[1]][,1],group[[8]][,1]]
    Landrace_Xinjiang_D_mean <- data.frame(ID=names(apply(Landrace_Xinjiang_D,1,mean)),IBS=as.numeric(apply(Landrace_Xinjiang_D,1,mean)))
    Landrace_Xinjiang_D_mean$Type = "Landrace_Xinjiang_D"
    H <- rbind(Landrace_Xinjiang_D_mean,Xinjiang)
    H_addLoc <- merge(H,location,by="ID")[,2:5]
    H_addLoc$IBS <- as.numeric(H_addLoc$IBS)
    melt <- melt(H_addLoc,id=c("Type","Latitude","Logititude"))
    H_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    H_od <- H_cast[order(H_cast$IBS),]
    H_od[1:5,1] <- "Min"
    
    Polish_Xinjiang_A <- ibs[[1]][group[[5]][,1],group[[8]][,1]]
    Polish_Xinjiang_B <- ibs[[2]][group[[5]][,1],group[[8]][,1]]
    Polish_Xinjiang_AB <- cbind(Polish_Xinjiang_A,Polish_Xinjiang_B)
    Polish_Xinjiang_AB_mean <- data.frame(ID=names(apply(Polish_Xinjiang_AB,1,mean)),IBS=as.numeric(apply(Polish_Xinjiang_AB,1,mean)))
    Polish_Xinjiang_AB_mean$Type = "Polish_Xinjiang_AB"
    I <- rbind(Polish_Xinjiang_AB_mean,Xinjiang)
    I_addLoc <- merge(I,location,by="ID")[,2:5]
    I_addLoc$IBS <- as.numeric(I_addLoc$IBS)
    melt <- melt(I_addLoc,id=c("Type","Latitude","Logititude"))
    I_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    I_od <- I_cast[order(I_cast$IBS),]
    I_od[1:5,1] <- "Min"
    
  #Persian
    Landrace_Persian_A <- ibs[[1]][group[[1]][,1],group[[4]][,1]]
    Landrace_Persian_B <- ibs[[2]][group[[1]][,1],group[[4]][,1]]
    Landrace_Persian_AB <- cbind(Landrace_Persian_A,Landrace_Persian_B)
    Landrace_Persian_AB_mean <- data.frame(ID=names(apply(Landrace_Persian_AB,1,mean)),IBS=as.numeric(apply(Landrace_Persian_AB,1,mean)))
    Landrace_Persian_AB_mean$Type = "Landrace_Persian_AB"
    J <- rbind(Landrace_Persian_AB_mean,Persian)
    J_addLoc <- merge(J,location,by="ID")[,2:5]
    J_addLoc$IBS <- as.numeric(J_addLoc$IBS)
    melt <- melt(J_addLoc,id=c("Type","Latitude","Logititude"))
    J_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    J_od <- J_cast[order(J_cast$IBS),]
    J_od[1:5,1] <- "Min"
    
    Rivet_Persian_A <- ibs[[1]][group[[6]][,1],group[[4]][,1]]
    Rivet_Persian_B <- ibs[[2]][group[[6]][,1],group[[4]][,1]]
    Rivet_Persian_AB <- cbind(Rivet_Persian_A,Rivet_Persian_B)
    Rivet_Persian_AB_mean <- data.frame(ID=names(apply(Rivet_Persian_AB,1,mean)),IBS=as.numeric(apply(Rivet_Persian_AB,1,mean)))
    Rivet_Persian_AB_mean$Type = "Rivet_Persian_AB"
    K <- rbind(Rivet_Persian_AB_mean,Persian)
    K_addLoc <- merge(K,location,by="ID")[,2:5]
    K_addLoc$IBS <- as.numeric(K_addLoc$IBS)
    melt <- melt(K_addLoc,id=c("Type","Latitude","Logititude"))
    K_cast <- cast(melt,Type+Latitude+Logititude~variable,mean)
    K_od <- K_cast[order(K_cast$IBS),]
    K_od[1:5,1] <- "Min"

    
#2. 画map图
data <- list(A_od,B_od,C_od,D_od,E_od,F_od,G_od,H_od,I_od,J_od,K_od)
#类型1
pdf("map.pdf",width = 12,height = 7.5)
for(i in c(1:11)){
  mp <- NULL
  mapworld <- borders("world",colour = "grey90",fill="white") 
  mp <- ggplot()+mapworld+ylim(-60,90)+xlim(-25,200)
  mp2 <- mp+geom_point(aes(x=data[[i]]$Logititude,y=data[[i]]$Latitude,shape= data[[i]]$Type, color=data[[i]]$IBS),size=2.6)+scale_size(range=c(1,1))+ 
    scale_colour_gradientn(colours = brewer.pal(11, "RdYlBu"))+
    theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
    theme(panel.grid = element_blank()) + 
    #theme(panel.grid =element_blank()) +   ## 删去网格线
    #theme(axis.text = element_blank()) +   ## 删去刻度标签
    #theme(axis.ticks = element_blank()) +   ## 删去刻度线
    theme(panel.border = element_blank())
    
  print(mp2+guides(fill=guide_legend(title=NULL)))
}
dev.off()

#类型2
ggplot(data[[1]], aes(Logititude, Latitude))+
  borders("world", xlim = c(-10,150), ylim = c(0, 90))+
  geom_point(aes(color=IBS,shape=Type),size=1.5)+
  #geom_point(aes(color=mean,shape=Value),size=3)+
  scale_colour_gradient(low = "yellow",high = "black")+
  scale_size_area()+
  coord_quickmap()+
  theme_classic()+
  theme(legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  guides(fill=guide_legend(title=NULL))


