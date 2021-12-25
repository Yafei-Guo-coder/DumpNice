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
brewer.pal(8, "Set1")
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

R CMD INSTALL clusterProfiler_4.0.5.tar.gz 

for(i in c(1:length(data))){
  name <- out1[[i]][1]
}

openstraightmap
dfe

#定义一个空的data.frame
all <- data.frame(CHROM="", BIN_START="", BIN_END="", MEAN_FST="", Gene_start="", Gene_end="",Gene_id="",Pop="", site="", lineage="", stringsAsFactors=FALSE)
p = ggplot(data, aes(x = taxaNum,y = mean))+
  geom_point()+ 
  geom_errorbar(aes(ymin=(mean-sd),ymax=(mean+sd)),width=2,size=2)+
  scale_fill_manual(values = c("red","blue"))+  
  labs(x = "Taxa Number",y = "Snp proportion") +
  facet_wrap(Type1~Type2,scales="free")+
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
ggplot(dataD, aes(x = V2)) +
  geom_histogram(binwidth = 1, fill = "lightblue", colour = "black")+
  theme_classic()
#箱线图
ggplot(Bdata, aes(x = B))+
    geom_boxplot(fill = '#f8766d', notch = TRUE)+theme_classic()+  theme(
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
