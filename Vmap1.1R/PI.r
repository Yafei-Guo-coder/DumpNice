#PI
#画图文件:yafei@159.226.116.204
#/data1/home/yafei/Project3/Pi/A/by_chr
#/data1/home/yafei/Project3/Pi/B/by_chr
#/data1/home/yafei/Project3/Pi/D/by_chr
#/data1/home/yafei/Project3/Pi/AB/by_chr
#/data1/home/yafei/Project3/Pi/ABD/by_chr
#画图代码
library(RColorBrewer)
library(ggplot2)
setwd("/data1/home/yafei/Project3/Pi/")
nameA <- c("Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian wheat","Khorasan wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic")
nameB <- c("Speltoides","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian wheat","Khorasan wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic")
nameD <- c("Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic","Strangulata", "Meyeri", "Anathera")
nameAB <- c("Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian_wheat","Khorasan_wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic")
nameABD <- c("Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic")

#A lineage
path <- "/data1/home/yafei/Project3/Pi/A/by_chr"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)
for(i in 1:length(data)){
  data[[i]]$Name <- nameA[i]
  if (i<4) {
    data[[i]]$Fill <- "pink"
  } else if (i<13) {
    data[[i]]$Fill <- "lightblue"
  } else {
    data[[i]]$Fill <- "orange"
  }
}
allA <- data.frame()
for(i in 1:length(data)){
  allA <- rbind(allA, data[[i]])
}
allA$lineage <- "A lineage"

#B lineage 
path <- "/data1/home/yafei/Project3/Pi/B/by_chr"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
for(i in 1:length(data)){
  data[[i]]$Name <- nameB[i]
  if (i<2) {
    data[[i]]$Fill <- "red"
  } else if (i<11) {
    data[[i]]$Fill <- "lightblue"
  } else {
    data[[i]]$Fill <- "orange"
  }
}
allB <- data.frame()
for(i in 1:length(data)){
  allB <- rbind(allB, data[[i]])
}
allB$lineage <- "B lineage"

#D lineage
path <- "/data1/home/yafei/Project3/Pi/D/by_chr"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)
for(i in 1:length(data)){
  data[[i]]$Name <- nameD[i]
  if (i<11) {
    data[[i]]$Fill <- "orange"
  } else {
    data[[i]]$Fill <- "yellow"
  }
}
allD <- data.frame()
for(i in 1:length(data)){
  allD <- rbind(allD, data[[i]])
}
allD$lineage <- "D lineage"

all <- rbind(allA,allB,allD)

P2<- ggplot(all, aes(x=Name, y=V4,fill=Fill),na.rm=TRUE) + 
  facet_grid(lineage ~ .) +
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  #  scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"))+
  #  scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"), breaks=c(breaks=c("red", "pink", "lightblue","orange","yellow")),labels = c(" SS", " AA", " AABB"," AABBDD"," DD"))+
  scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"),labels = c(" SS", " AA", " AABB"," AABBDD"," DD"))+
  scale_x_discrete(limits= c("Speltoides","Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Rivet wheat","Polish wheat","Durum","Khorasan wheat","Persian wheat","Ispahanicum","Georgian wheat","Indian dwarf wheat","Yunan wheat","Vavilovii","Macha","Spelt","Tibetan semi-wild","Xinjiang wheat","Club wheat","Cultivar","Landrace","Strangulata", "Meyeri", "Anathera"))+
  theme_bw()+ #背景变为白色
  #B:0.01 D:0.004 AB:0.01
  ylab("PI")+xlab("Lineages")+scale_y_continuous(limits=c(0,0.02)) + 
  theme_classic()+
  theme(legend.text = element_text(size = 20),
        legend.key.size=unit(1,'cm'),
        legend.key.width=unit(0.8,'cm'),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 30,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_text(size=20))

pdf(file = "Pi.pdf",width=20,height=10) #结果保存
print(P2)
dev.off()

#画每组比较的boxplot图

P3 <- ggplot(all, aes(Name, V4)) + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA, aes(fill=factor(lineage)))+
  scale_fill_manual(name = "Lineages", values = c("#66C2A5","#FC8D62","#FFD92F"), breaks=c(breaks=c("A lineage", "B lineage", "D lineage")),labels = c("A lineage", " B lineage", " D lineage"))+
  guides(fill=guide_legend(title=NULL)) +
  scale_x_discrete(limits= c("Wild einkorn","Domesticated einkorn","Urartu", "Speltoides","Wild emmer","Domesticated emmer","Rivet wheat","Polish wheat","Durum","Khorasan wheat","Persian wheat","Ispahanicum","Georgian wheat","Indian dwarf wheat","Yunan wheat","Vavilovii","Macha","Spelt","Tibetan semi-wild","Xinjiang wheat","Club wheat","Cultivar","Landrace","Strangulata", "Meyeri", "Anathera"))+
  theme_bw()+ #背景变为白色
  ylab("PI")+xlab("Lineages")+scale_y_continuous(limits=c(0,0.02)) + 
  theme_classic()+
  theme(legend.text = element_text(size = 20),
        legend.key.size=unit(1,'cm'),
        legend.key.width=unit(0.8,'cm'),
        legend.background = element_rect(),
        legend.position=c(0.85,0.85),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 30,face="plain"))

pdf(file = "Pi2.pdf",width=20,height=10) #结果保存
print(P3)
dev.off()

