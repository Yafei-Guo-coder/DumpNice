#调色板：display.brewer.pal(n = 8, name = 'Set2')
#yafei@204:/data1/home/yafei/Project3/Pi/
#ggplot2画图
#改名字
library(ggplot2)
nameA <- c("Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian wheat","Khorasan wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic")
nameB <- c("Speltoides","Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian wheat","Khorasan wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic")
nameD <- c("Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic","Strangulata", "Meyeri", "Anathera")
nameAB <- c("Wild emmer","Domesticated emmer","Georgian wheat","Ispahanicum","Rivet wheat","Polish wheat","Persian_wheat","Khorasan_wheat","Durum","Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic")
nameABD <- c("Spelt","Macha","Club wheat","Indian dwarf wheat","Yunan wheat","Xinjiang wheat","Tibetan semi-wild","Vavilovii","Landrace","Cultivar","Synthetic")

#修改 路径和name以及输出文件名
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
colnames(all)[1:6] <- c("Chr","Start","End","Pi","Value","Text")
#A: 4, 13
#B: 2, 11
#D: 11
#AB: 10
#ggplot(data, aes(x=Maf))+ geom_histogram() + facet_grid(. ~ Chr) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
library(RColorBrewer)
#绘图：修改横坐标顺序和输出文件名
library(ggplot2)

p <- data %>%
  #mutate(text = fct_reorder(Text, Value)) %>% # Reorder data
  ggplot(aes(x=Value, y=Pi, fill=lineage, color=lineage)) +
  geom_violin(outlier.shape = NA, outlier.color = NA,trim=TRUE,width=1) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  #theme_ipsum() +
  theme(
    legend.position="none"
  ) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  xlab("") + ylim(0,0.02)+ facet_grid(.~lineage)+
  ylab("Assigned Probability (%)")

P2<- ggplot(all, aes(x=Text, y=Value,fill=Fill),na.rm=TRUE) + 
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
  theme(legend.text = element_text(size = 200),
        legend.key.size=unit(10,'cm'),
        legend.key.width=unit(8,'cm'),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=200),
        axis.text.y=element_text(size=200),
        axis.title.y=element_text(size = 300,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_text(size=200))

pdf(file = "Pi.pdf",width=200,height=100) #结果保存
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
  theme(legend.text = element_text(size = 200),
        legend.key.size=unit(10,'cm'),
        legend.key.width=unit(8,'cm'),
        legend.background = element_rect(),
        legend.position=c(0.85,0.85),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=200),
        axis.text.y=element_text(size=200),
        axis.title.y=element_text(size = 300,face="plain"))

pdf(file = "Pi2.pdf",width=200,height=100) #结果保存
print(P3)
dev.off()

ggplot(data, aes(x =Pi, y=Value,fill=lineage)) +
  #scale_y_discrete(limits= c("Speltoides","Wild einkorn","Domesticated einkorn","Urartu","Wild emmer","Domesticated emmer","Rivet wheat","Polish wheat","Durum","Khorasan wheat","Persian wheat","Ispahanicum","Georgian wheat","Indian dwarf wheat","Yunan wheat","Vavilovii","Macha","Spelt","Tibetan semi-wild","Xinjiang wheat","Club wheat","Cultivar","Landrace","Strangulata", "Meyeri", "Anathera"))+
  geom_density_ridges() +
  theme_ridges() + xlim(0,0.005)+
  theme(legend.position = "none")


data %>%
  mutate(text = fct_reorder(Value, Pi)) %>%
  ggplot( aes(y=Value, x=Pi,  fill=Value)) +
  geom_density_ridges(alpha=0.6, stat="binline", bins=200) +
  theme_ridges() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.001, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") + xlim(0,0.02)+
  ylab("Assigned Probability (%)")

library(tidyverse)

# Create dataset
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
  value=sample( seq(10,100), 60, replace=T)
)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(all, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  
  geom_segment(data=grid_data, aes(x = end, y = 0.15, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.1, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.05, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),3), y = c(0.05,0.1,0.15), label = c("0.05", "0.1", "0.15") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(0,1) +
  theme_minimal()+ 
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
p

setwd("/Users/guoyafei/Documents/Lulab/Project-2-Migration/基本统计/Plots/Drift/")
dataA <- read.table("A_drift.txt",header=F,stringsAsFactors = F)
dataA$group <- "A"
dataB <- read.table("B_drift.txt",header=F,stringsAsFactors = F)
dataB$group <- "B"
dataD <- read.table("D_drift.txt",header=F,stringsAsFactors = F)
dataD$group <- "D"
all$id <- c(1:64)
colnames(all) <- c("individual","value","group","id")
data<-all
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.12, xend = start, yend = 0.12), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.09, xend = start, yend = 0.09), colour = "grey", alpha=1, size=0.5, inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.06, xend = start, yend = 0.06), colour = "grey", alpha=1, size=0.5, inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.03, xend = start, yend = 0.03), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0.03, 0.06, 0.09, 0.12), label = c("0.03", "0.06", "0.09", "0.12") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) + 
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-0.1,0.2) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.01, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.01, xend = end, yend = -0.01), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -0.02, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
p



#计算pi的平均值和方差
library(tidyverse)
library(reshape2)
setwd("/Users/guoyafei/Documents/Lulab/Project-2-Migration/基本统计/Plots/PI/")
data <- read.table("Pi.txt",header=T,stringsAsFactors = F,sep="\t")
pi <- data[which(data$Pi>0),c(4:5,7)]
mydata<-melt(pi,id=c("Value","lineage"))
mean <- dcast(mydata,Value+lineage~variable,mean)
sd <- dcast(mydata,Value+lineage~variable,sd)
colnames(mean)[3] <- "mean"
colnames(sd)[3] <- "sd"
all <- as.data.frame(cbind(mean,sd[,3]))
colnames(all)[4] <- "sd"
all$id <- c(1:58)
colnames(all) <- c("individual","group","value","value2","id")
data <- all
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))


data <- read.table("Pi_stat.txt",header=T,stringsAsFactors = F,sep="\t")
data <- read.table("PI_mean_sd.txt",header=T,stringsAsFactors = F,sep="\t")
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))
# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]



pd <- position_dodge(0.1)
p <- ggplot(data2, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.12, xend = start, yend = 0.12), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.09, xend = start, yend = 0.09), colour = "grey", alpha=1, size=0.5, inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.06, xend = start, yend = 0.06), colour = "grey", alpha=1, size=0.5, inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.03, xend = start, yend = 0.03), colour = "grey", alpha=1, size=0.5 , inherit.aes = FALSE ) +
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0.03, 0.06, 0.09, 0.12), label = c("0.003", "0.006", "0.009", "0.012") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) + 
  #geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1, position=pd) +
  geom_point(position=pd)+
  ylim(-0.03,0.03) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.01, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.01, xend = end, yend = -0.01), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -0.02, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
p