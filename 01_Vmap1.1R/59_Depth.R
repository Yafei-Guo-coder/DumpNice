#depth.r
library(cowplot) 
library(ggExtra)
library(forcats)
library(ggplot2)
library(RColorBrewer)
library(aplot)
setwd("/Users/guoyafei/Desktop/NP/reviewer/depth")
# method 1 ------------
dataFile <- c("tetra_depth.txt.gz")
df1 <- read.delim(dataFile) %>% mutate(GenomeType = "AABB")
dataFile <- c("hexa_depth.txt.gz")
df2 <- read.delim(dataFile)%>% mutate(GenomeType = "AABBDD")
dataFile <- c("diploid_depth.txt.gz")
df3 <- read.delim(dataFile)%>% mutate(GenomeType = "DD")

df <- rbind(df1,df2,df3) %>% 
  group_by(GenomeType) %>% 
  slice_sample(.,n = 3000) %>% 
  ungroup()

data <- read.table("all_100k.depth",header=F,stringsAsFactors = F)
colnames(data) <- c("Depth","SD")

colB <- c('#ffd702','#fc6e6e',"#87cef9")
df$GenomeType <- factor(df$GenomeType,levels = c("AABB","AABBDD","DD"))

p <- ggplot(data,aes(x=AverageDepth,y=SD))+
  geom_point(alpha=0.5,shape=1,size=0.8)+
  xlab("Mean of depth")+ylab("SD")+
  #scale_color_manual(values = colB)+
  xlim(0,15)+ylim(0,10)+
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.7)
    ,legend.position = 'none'
    ,text = element_text(size = 14))

p
p2 <- ggMarginal(p, type="density",groupColour = F, groupFill = F)
p2
ggsave(p2,filename = "~/Documents/test.pdf",height = 4,width = 4)

# method 2 ------------
data <- read.table("all_100k.depth",header=F,stringsAsFactors = F)
colnames(data) <- c("Depth","SD")
data <- data[which(data$Depth > 0 ),]
data <- data[which(data$Depth < 15 & data$SD < 10),]
data <- data[sample(1:85255, 30000),]
  p1<-ggplot(data,aes(x = Depth, y = SD))+
  geom_point(alpha= 0.1)+
  scale_color_brewer(palette = "Set2")+
  #scale_color_manual(values=c("green","blue","grey"))+
  theme_bw()
#geom_hline(yintercept = 3.86971,lty="dashed")+
#geom_hline(yintercept = 6.15844,lty="dashed")+
#geom_vline(xintercept = 3.91843,lty="dashed")+
#geom_vline(xintercept = 9.36970,lty="dashed")
p2<-ggplot(data,aes(Depth))+
  geom_density(fill="grey",alpha=0.5)+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p3<-ggplot(data,aes(SD))+
  geom_density(fill="grey",alpha=0.5)+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  coord_flip()
#累计分布函数
p4<-ggplot(data, aes(Depth)) + 
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_color_brewer(palette = "Set2")+
  stat_ecdf()+
  scale_y_reverse()
p5<-ggplot(data, aes(SD)) + 
  scale_color_brewer(palette = "Set2")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  stat_ecdf()+
  coord_flip()+
  scale_y_reverse()
p1%>%
  insert_top(p2,height = 0.3)%>%
  insert_right(p3,0.3)

