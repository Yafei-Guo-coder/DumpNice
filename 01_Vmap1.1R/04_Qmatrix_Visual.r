#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)
require(gridExtra)
setwd("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/01_BasicStatistic/06_Structure/V1/Select")
#load Qmatrix files
sfiles <- list.files(path="/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/01_BasicStatistic/06_Structure/V1/Select/QMatrix", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))

pdf("Q.pdf",width = 24,height = 8)
p1 <- plotQ(slist[1:7],returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="all",showindlab=F,showyaxis=T,showticks=T,sharedindlab=T,clustercol=brewer.pal(7, "Set2"))
grid.arrange(p1$plot[[1]],p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],nrow=5)
dev.off()

library(readxl)
library(ggmap)
library(RColorBrewer)
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/Fst/")
group <- read.table("group.txt",header=F,stringsAsFactors = F)
a <- read_excel("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx", sheet=1, na='NA')
b<-a[which(a$`Taxa(vmap3)` %in% group$V3),c(1,11,12)]
rownames(group)<-group$V3
c <- cbind(b,group[b$`Taxa(vmap3)`,])

mp <- NULL
mapworld <- borders("world",colour = "gray70",fill="white") 
mp <- ggplot() + 
  mapworld + 
  #xlim(0,90) +
  #ylim(10,70) + 
  theme_classic()
color <- brewer.pal(8, "Dark2")[c(1,2,3,4,5,6,7)]
color <- brewer.pal(8, "Dark2")[1:7]
data <- c[which(c$V1=="domemmer"),]
data <- c[which(c$V1=="freethresh"),]
data <- c[which(c$V1=="wildemmer"),]
#AABBDD
sub <- data[which(data$Region %in% c("EU","WA","EA2","SH","IA")),]
mp+geom_point(aes(x=sub$Longitude, y=sub$Latitude, color=sub$Region,alpha=0.3),size=1)+
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
        panel.border = element_blank())







