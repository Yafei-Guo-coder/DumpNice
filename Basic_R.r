#palette using black
cbbPalette <- c(
  
branch: 乌拉尔图: "#999999"（灰色）, 野生一粒:"#E69F00"(橘黄色), 栽培一粒:"#56B4E9"(明兰色), AABBDD:"#009E73"(橄榄绿), AABB:"#F0E442"（明黄）
label: WA: "#000000"（黑色）, EU："#0072B2"(天蓝), EA："#D55E00"(橘红色), SCA："#CC79A7"(皮粉色)),AF: "#E78AC3"(粉色) AM:"#66C2A5"(浅绿) "#A07D35"(浅棕色)

#循环：for
v <- LETTERS[1:6]
for ( i in v) {
  if (i == "D") { 
    next
  }
  print(i)
}
#连接两个字符串：paste
paste (a, sep = " ", collapse = NULL)
paste0(a, collapse = NULL)
#定义一个空向量
x=vector()
x<-numeirc(0) #长度可变的存储数字的向量
x=character() #创建出来的为字符串向量
x[1]=1
x[2]=3
x<-NULL; x[1]<-2
vector(mode="numeric",length=0) #定义一个空向量，往里面添加元素即可
x=matrix(nrow = 2,ncol=3) #创建空矩阵
#添加元素
c1 <- c(1,2,3,4,5) #创建一个向量
c1 <- c(c1,5) #追加一个元素
c1 <- c(c1,c(5,6)) #追加一个向量
c1 <- c(c1[1:2],c(5,6),c[3:5]) #指定位置来添加的元素
c1 <- append(c1,8) #在向量最后追加一个元素8
c1 <- append(c1,c(11,22)) #在向量后追加向量
c1 <- append(c1,35,3) #在第3个元素后插入新元素，也可以插入向量

#删除元素
c1 <- c1[-1] #从向量中指定位置为1的元素
c1 <- c1[-c(2:3)] #可以给定一个位置向量来删除多个元素
c1 <- c1[c(3:5)] #与上面的方式相反，保留想要的元素
#展示调色板
library(RColorBrewer)
display.brewer.all()
brewer.pal(9, "Set1")
scale_fill_brewer(palette = "RdYiBu")
scale_fill_manual(values = c("#FDC086","#BEAED4")) 
scale_color_manual(values = color) 
scale_colour_discrete(breaks = c("#838B8B","#FFD700", "#97FFFF", "#D8BFD8", "#FF6349"), labels = c('EU','WA','SCA','EA-N','EA-S'))

#批量读取文件
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/VIP_gene/V2/snpEff/snpEff" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep="\t")})

library(reshape)
melt <- melt(A_addLoc,id=c("ID","Type","Latitude","Logititude"))

p <- ggplot(all, aes(RDA1, RDA2,color=RDA_Region)) +
  geom_point( size=3) +
  #geom_point( size=3) +
  #stat_ellipse(aes(group=label$cols), level = 0.95, show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  #scale_color_brewer(brewer.pal(6, "Set2")[c(1,2,3,4,6)])+
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
  #theme(panel.grid =element_blank()) +   ## 删去网格线
  #theme(axis.text = element_blank()) +   ## 删去刻度标签
  #theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank()) +
  #theme(panel.grid.major = element_line(color = 'gray', size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'))+
  #legend.title = (element_blank(), legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'RDA1 (34.96%)', y = 'RDA2 (11.23%)') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb_forward_r.env, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.4, 'cm')), size = 1, color = 'brown',alpha=0.5) +
  geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'brown', size = 6)+
  #scale_colour_discrete(breaks = c("#838B8B","#FFD700", "#97FFFF", "#D8BFD8", "#FF6349"), labels = c('EU','WA','SCA','EA-N','EA-S'))+
  guides(fill=guide_legend(title=NULL))
  #geom_label_repel(aes(label =sample, color = group), size = 3, box.padding = unit(0, 'lines'), show.legend = FALSE)

geom_smooth(method = "lm", color = "black", fill = "lightgray") 
ggtitle("North1 VS North2")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

out <- strsplit(sub('_',':', dist[,1]) , ":")
data$tissue <- unlist(lapply(X = data$V6, FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))
R CMD INSTALL clusterProfiler_4.0.5.tar.gz 

for(i in c(1:length(data))){
  name <- out1[[i]][1]
}

openstraightmap
dfe

#定义一个空的data.frame
all <- data.frame(CHROM="", BIN_START="", BIN_END="", MEAN_FST="", Gene_start="", Gene_end="",Gene_id="",Pop="", site="", lineage="", stringsAsFactors=FALSE)
p = ggplot(data, aes(x = taxaNum,y = mean)) +
  geom_point() + 
  geom_errorbar(aes(ymin=(mean-sd),ymax=(mean+sd)),width=2,size=2) +
  scale_fill_manual(values = c("red","blue")) +  
  labs(x = "Taxa Number",y = "Snp proportion") +
  facet_wrap(Type1~Type2,scales="free") +
  theme_bw() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  )+
  theme(legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 10))
  #geom_text(aes(label = dat$Num),position=position_dodge(width = 0.5),size = 5,vjust = -0.25)+ ##########

#聚类方法：
ward最小方差法

#读Excel文件
library(gdata)
read.xls("/home/slave/test.xls",sheet=1,na.strings=c("NA","#DIV/0!"))

result <- substring("Extract", 5, 7)
print(result)

#画直方图
data <- read.table("~/Desktop/all.ibs",header=F,stringsAsFactors = F)
ggplot(dataD, aes(x = V2)) +
  geom_histogram(binwidth = 1, fill = "lightblue", colour = "black")+
  theme_classic()
N <- data[which(data$V6 == "Wild_emmer_N"),]
S <- data[which(data$V6 == "Wild_emmer_S"),]
ggplot(N, aes(x = V2,group = V6)) +
  geom_histogram(binwidth = 0.001, fill = "lightblue", colour = "black")+
  theme_classic()

#箱线图
ggplot(Bdata, aes(x = B))+
    geom_boxplot(fill = '#f8766d', notch = TRUE) +
    theme_classic() +  
    theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
) 

ggplot() + geom_bar(data = data_order, aes(x = x, y=y),stat = "identity")

data_order %>%
  mutate(name = fct_reorder(name, desc(y))) %>%
  ggplot(aes(x=name, y=y)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") + ylab("")
theme_bw()

data_sorted <- data %>% 
  mutate(class = fct_reorder(x, y))

data_order$name <- factor(data_order$x,levels = unique(data_order$x))

colnames(data) <- c("taxa","elevation",paste("temp",seq(1,11),sep=""),paste("prec",seq(1,8),sep=""),c("latitude","longitude","population1","population"))

data <- read.table("mergeall_env.txt",header=F,stringsAsFactors = F)
p <- list()
p[[1]] <- ggscatter(data, x = "longitude", y = "elevation", add = "reg.line", color = "population") + stat_cor(label.x = 10, aes(color = population,label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) + theme( legend.position="none")
library(gridExtra)
pdf("longidute3.pdf",width = 22,height = 22)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],nrow=4)
dev.off()

ggbuild + ggplot (抠图上的点)

#画地图
library("sf")
ggplot() +
  geom_sf(data = world,color = "dark grey", fill = "white",alpha=1) +
  geom_point(data = LR_clust,aes(x=Longitude,y=Latitude,color=Group_BayPass)) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("442 landrace distribution for GEA analysis") +
  coord_sf(xlim = c(0, 140),
           ylim = c(10, 60),
           expand = T)

a <- vector()
for ( i in c(1:30)) {
  a <- append(a,as.numeric(dpois(0, i, log = FALSE)))
}
b <- c(1:30)
data <- as.data.frame(cbind(b,a))

psmc <- read.table("/Users/guoyafei/Desktop/NP/PSMC_miss.txt",header=F,stringsAsFactors = F)
vmap <- read.table("/Users/guoyafei/Desktop/NP/Chr1A_snpcount_missing.txt", header=F,stringsAsFactors = F)

colnames(data) <- c("coverage","missing")
data$type <- "possion"
colnames(vmap) <- c("coverage","missing")
vmap$type <- "vmap"
colnames(psmc) <- c("coverage","missing")
psmc$type <- "psmc"
all <- rbind(data,vmap,psmc)

brewer.pal(9, "Set11")
scale_fill_brewer(palette = "RdYiBu")
scale_fill_manual(values = c("#FDC086","#BEAED4")) 
scale_color_manual(values = color) 

ggplot(data = all, 
       aes(x = coverage, y = missing, color = type)) +
  geom_point() +
  geom_smooth(se = TRUE, method = "gam", formula = y ~ s(log(x)))+
  theme_bw()+
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB"))
  
####wanghao
setwd("/Users/guoyafei/Documents/02_VmapIII/13_Haplotype/wanghao")
data <- read.table("Chr_A_dom.xpclr.TE.bed",header=F,stringsAsFactors = F)
ggplot(data, aes(x=V2, y=V4)) +
  geom_point(size=0.5)+
  geom_line()+
  theme_bw()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("XP-CLR") +
  geom_hline(yintercept=3.37,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759451668,color='red')
  
data <- read.table("Chr_AB_dom.xpclr.TE.bed",header=F,stringsAsFactors = F)
ggplot(data, aes(x=V2, y=V4)) +
  geom_point(size=0.5)+
  geom_line()+
  theme_bw()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("XP-CLR") +
  geom_hline(yintercept=1.908,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759451668,color='red')

data <- read.table("Chr_AB_imp.xpclr.TE.bed",header=F,stringsAsFactors = F)
ggplot(data, aes(x=V2, y=V4)) +
  geom_point(size=0.5)+
  geom_line()+
  theme_bw()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("XP-CLR") +
  geom_hline(yintercept=3.09,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759451668,color='red')

3.63

data <- read.table("Chr_ABD_imp.xpclr.TE.bed",header=F,stringsAsFactors = F)
ggplot(data, aes(x=V2, y=V4)) +
  geom_point(size=0.5)+
  geom_line()+
  theme_bw()+
  #theme(axis.title.y = element_blank()) 
  xlab("Position") + ylab("XP-CLR") +
  geom_hline(yintercept=3.63,color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=759451668,color='red')
防止列名改变
check.names=F,

setwd("/Users/guoyafei/Desktop/")
data <- read.table("/Users/guoyafei/Desktop/AB.k10.cultivar.salti.lassip.mlg.txt.shuf",header=F,stringsAsFactors = F)
ggplot(data,aes(V10))+
  geom_density(fill="grey",alpha=0.5)+
  #scale_y_continuous(expand = c(0,0))+
  theme_bw()
ggplot(data,aes(V11))+
  geom_density(fill="grey",alpha=0.5)+
  #scale_y_continuous(expand = c(0,0))+
  theme_bw()
ggplot(data,aes(V12))+
  geom_density(fill="grey",alpha=0.5)+
  #scale_y_continuous(expand = c(0,0))+
  theme_bw()

data2 <- read.table("/Users/guoyafei/Desktop/AB.k10.cultivar.genes.txt",header=T,stringsAsFactors = F)
ggplot(data2,aes(m))+
  geom_density(fill="grey",alpha=0.5)+
  #scale_y_continuous(expand = c(0,0))+
  theme_bw()
ggplot(data2,aes(A))+
  geom_density(fill="grey",alpha=0.5)+
  #scale_y_continuous(expand = c(0,0))+
  theme_bw()
ggplot(data2,aes(T))+
  geom_density(fill="grey",alpha=0.5)+
  #scale_y_continuous(expand = c(0,0))+
  theme_bw()

data <- read.table("/Users/guoyafei/Desktop/adapgene/prec7.alt.raw.out", header=F,stringsAsFactors = F)
data$maf <- NA
data[which(data$V4 < 0.5),5] <- data[which(data$V4 < 0.5),4]
data[which(data$V4 >= 0.5),5] <- 1-data[which(data$V4 >= 0.5),4]
data$effect <- NA
data[,6] <- abs(data[,3])
ggplot(data, aes(maf,effect)) +
  geom_point( size=3,alpha=0.3,color="blue") +
  theme_classic()

ggplot(data, aes(V7,V8)) +
  geom_point( size=3,alpha=0.3,color="#A65628") +
  ylim(-150,150)+
  theme_classic()

data <- read.table("/Users/guoyafei/Desktop/adapgene/shuf1k.gcta.out", header=F,stringsAsFactors = F)
data <- read.table("/Users/guoyafei/Desktop/adapgene/sig.shuf1k.gcta.out", header=F,stringsAsFactors = F)
data <- read.table("/Users/guoyafei/Desktop/adapgene/prec7.sig.gcta.out", header=F,stringsAsFactors = F)

#一页多图
library(ggplot2)
library(gplots)
library(RColorBrewer)
library("corrplot")
library("cowplot")
library("gridExtra")
library(ggpubr)
setwd("/Users/guoyafei/Desktop/coFu/effect")
data <- read.table("/Users/guoyafei/Desktop/coFu/effect/TraesCS5A02G391300.5k.all.txt",header=T,stringsAsFactors = F)
all <- read.table("/Users/guoyafei/Desktop/coFu/effect/TraesCS5A02G391300.5k.crosstab.txt",check.names=F,header=T,row.names=1, stringsAsFactors = F)
ld <- as.matrix(all)
p <- list()
p[[1]]<-ggplot(data = data, aes(x = pos, y = -log10(gwasp), color = snpeff)) +
  geom_point() +
  theme_bw()
p[[2]]<-ggplot(data = data, aes(x = pos, y = r2, color = -log10(p))) +
  geom_point() +
  theme_bw()
pdf("test.pdf",width=10,height=5)
ggarrange(p[[1]], p[[2]], heights = c(1, 1),
          ncol = 1, nrow = 2, align = "v")
corrplot(ld,method = "color",tl.col="black",tl.srt = 45, addrect = 2,addCoef.col= NULL, type = "upper")
dev.off()

setwd("~/Desktop/maf/")
library(ggplot2)
data <- read.table("shuf800k_1A.frq.all.txt",header=F,stringsAsFactors = F)
ggplot(data,aes(V2,fill=V3))+
     geom_histogram(bins=30,aes(x=V2, y=..density..),position="identity",alpha = 0.5)

#数据整合
M <- aggregate(data, by=list(type),FUN=mean)
#长格式变宽格式
cats <- cast(all,id1~variable,value.var="value",fun.aggregate = mean)
#宽格式变长格式
cats <- melt(data,id="type")
#相关性分析
corrgram(cats, order=F, lower.panel=panel.shade, upper.panel=NULL, text.panel=panel.txt, main="Car Milage Data in PC2/PC1 Order")
heatmap(cats, Colv = NA, Rowv = NA, scale="column")
pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, cluster_col = F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"),cluster_row = FALSE, annotation_col = anno2, annotation_names_col = F)

corrplot(cats,method = "color",col.lim = c(20, 30),type = 'lower',tl.col="black",tl.srt = 45,addrect=1,addCoef.col = "grey",number.cex=0.5,number.digits=2,tl.cex=1,cl.cex=1,cl.lim = c(0, 1))

ggplot(data,aes(a,b))+
  theme_bw()+
  geom_bar()

ggplot(all,aes(beta)) +
  theme_bw() +
  geom_histogram(binwidth = 50,aes(x=beta, y=..density..),position="identity",alpha = 0.5)

a <- data[,4,drop=F]
colnames(a) <- "beta"
b <- data[,5,drop=F]
colnames(b) <- "beta"

all <- as.data.frame(cbind(pop,p005,p001))

pdf("sites.pdf")

p <- ggplot(data, aes(a,b)) +
  geom_bar(stat="identity") +
  theme_classic() + 
  geom_text(aes(label=b)) +
  labs(x = 'Number of groups in which a site under selection', y = 'Number of sites')

data <- read.table("~/Desktop/test3.txt", header=T,stringsAsFactors = F)
sub <- data[which(data$RECEIVED.USDA. != "NA"),]

sub <- data[which(data$SOURCE_DATE.USDA. != "NA"),]

ggplot(data,aes(y=V2,group=V1))+
  theme_bw()+
  geom_histogram()

data <- read.table("/Users/guoyafei/Desktop/type1.mode.salti.out", header=F, stringsAsFactors = F)
data$V1 <- as.factor(data$V1)
data[which(data$V2 > 1),2 ] <- 2
data$type<- 0
data[which(data$V1 > 1 & data$V1 < 7 ),4] <- 1
data[which(data$V1 > 6 & data$V1 < 14 ),4 ] <- 2
data[which(data$V1 > 13 & data$V1 < 21 ),4 ] <- 3
data$V1 <- data$type

setwd("/Users/guoyafei/Desktop/3type/选择/")
library(ggplot2)
library(reshape)
data <- read.table("pop1.type1.mode.salti.out", header = F, stringsAsFactors = F)
data2 <- cast(data,V2~V3)
data3 <- melt(data2,id="V1")
data4 <- data3[which(data3$V2 != 0),]
data5 <- data4[which(data3$V1 != 0),]
pdf("test.pdf")
p <- ggplot(data=data5, mapping=aes(x = V1, y = value,fill=V2))+
  geom_bar(stat="identity",position=position_dodge(0.75))+
  theme_classic()
print(p)
dev.off()
library(corrgram)
ggcorr(clim, method = c("everything", "pearson"), nbreaks = 5,hjust = 1, size = 4, color = "grey50",layout.exp = 1, legend.position = "bottom", legend.size = 12) +
  guides(fill = guide_colorbar(barwidth = 18, title.vjust = 0.75)) +
  theme(legend.title = element_text(size = 14))







