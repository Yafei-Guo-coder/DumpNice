#调色板：display.brewer.pal(n = 8, name = 'Dark2')
#工作目录：203:/data2/yafei/Project3/TajimasD/VmapData
#绘图：修改输出文件名
library(ggplot2)
#改名字
nameA <- c("Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian wheat","Khorasan wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar")
nameB <- c("Speltoides","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian wheat","Khorasan wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar")
nameD <- c("Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Strangulata", "Meyeri", "Anathera")
nameAB <- c("Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian_wheat","Khorasan_wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar")
nameABD <- c("Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar")

#修改 路径和name以及输出文件名
#A lineage
path <- "/data2/yafei/Project3/TajimasD/VmapData/A"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 

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
#Blineage
path <- "/data2/yafei/Project3/TajimasD/VmapData/B"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 

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
#Dlineage
path <- "/data2/yafei/Project3/TajimasD/VmapData/D"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 

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

P2<- ggplot(all, aes(x=Name, y=TajimaD,fill=Fill),na.rm=TRUE) + 
  facet_grid(lineage ~ .) +
  #geom_violin() +
  scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"), breaks=c(breaks=c("red", "pink", "lightblue","orange","yellow")),labels = c(" SS", " AA", " AABB"," AABBDD"," DD"))+
  #guides(fill=guide_legend(title=NULL)) +
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  scale_x_discrete(limits= c("Speltoides","Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Rivet wheat","Polish wheat","Durum","Khorasan wheat","Persian wheat","Ispahanicum","Georgian wheat","Indian dwarf wheat","Yunan wheat","Vavilovii","Macha","Spelt","Tibetan semi-wild","Xinjiang wheat","Club wheat","Cultivar","Landrace","Strangulata", "Meyeri", "Anathera"))+
  theme_bw()+ #背景变为白色
  #B:0.01 D:0.004 AB:0.01
  ylab("TajimaD")+xlab("Lineages")+scale_y_continuous() + 
  theme_classic()+
  theme(legend.text = element_text(size = 200),
        legend.key.size=unit(10,'cm'),
        legend.key.width=unit(8,'cm'),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=200),
        axis.text.y=element_text(size=200),
        axis.title.y=element_text(size = 300,face="plain"),
        strip.text.y = element_text(size=200),
        strip.text.y = element_text(margin = margin(.1, 0, .1, 0, "cm")))+
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")

pdf(file = "TajimaD.pdf",width=200,height=100) #结果保存
print(P2)
dev.off()

#画每组比较的boxplot图

P3 <- ggplot(all, aes(x=Name, y=TajimaD)) + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA, aes(fill=factor(lineage)))+
  scale_fill_manual(name = "Lineages", values = c("#66C2A5","#FC8D62","#FFD92F"), breaks=c(breaks=c("A lineage", "B lineage", "D lineage")),labels = c(" A lineage", " B lineage", " D lineage"))+
  guides(fill=guide_legend(title=NULL)) +
  scale_x_discrete(limits= c("Wild einkorn","Domesticated einkorn","Urartu", "Speltoides","Wild emmer","Domesticated emmer","Rivet wheat","Polish wheat","Durum","Khorasan wheat","Persian wheat","Ispahanicum","Georgian wheat","Indian dwarf wheat","Yunan wheat","Vavilovii","Macha","Spelt","Tibetan semi-wild","Xinjiang wheat","Club wheat","Cultivar","Landrace","Strangulata", "Meyeri", "Anathera"))+
  theme_bw()+ #背景变为白色
  ylab("TajimaD")+xlab("Lineages")+scale_y_continuous() + 
  theme_classic()+
  theme(legend.text = element_text(size = 200),
        legend.key.size=unit(10,'cm'),
        legend.key.width=unit(8,'cm'),
        legend.background = element_rect(),
        legend.position=c(0.5,0.85),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=200),
        axis.text.y=element_text(size=200),
        axis.title.y=element_text(size = 300,face="plain")) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")

pdf(file = "TajimaD2.pdf",width=200,height=100) #结果保存
print(P3)
dev.off()