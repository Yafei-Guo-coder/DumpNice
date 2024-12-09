#################################################################适应性位点在全基因组上的物理分布################################################################
setwd("/Users/guoyafei/Desktop/3type/location")
library(CMplot)
end <- read.table("chrend_21chr.txt",header=F,stringsAsFactors = F)
colnames(end) <- c("snp","chr","pos")
setwd("/Users/guoyafei/Desktop/3type/location")
mydata <- read.table("type1.21chr.bed", header=F, stringsAsFactors = F)
mydata <- read.table("type2.21chr.bed", header=F, stringsAsFactors = F)
mydata <- read.table("type3.21chr.bed", header=F, stringsAsFactors = F)

mydata <- read.table("rank04_pop15_pop7_2.745_0.709_2.036.pos.21chr.txt", header=F, stringsAsFactors = F)
colnames(mydata) <- c("chr","start","pos","snp")
data <- mydata[,c(4,1,3)]
all <- rbind(data,end)
head(data)
# snp         chr       pos
# snp1_1    1        2041
CMplot(data,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="pdf",memo="snp_density",dpi=300)
#all position
mydata1 <- read.table("type1.21chr.bed", header=F, stringsAsFactors = F)
mydata2 <- read.table("type2.21chr.bed", header=F, stringsAsFactors = F)
mydata3 <- read.table("type3.21chr.bed", header=F, stringsAsFactors = F)
colnames(mydata1) <- c("chr","start","pos","snp")
colnames(mydata2) <- c("chr","start","pos","snp")
colnames(mydata3) <- c("chr","start","pos","snp")
my2 <- mydata2[,c(4,1,3)]
my2$type <- "type2"
end$type <- "type2"
type2 <- rbind(my2,end)
my3 <- mydata3[,c(4,1,3)]
my3$type <- "type3"
end$type <- "type3"
type3 <- rbind(my3,end)
all <- rbind(type1,type2,type3)
#但是我不会把他们画在同一条染色体上

#位点全基因组密度分布
setwd("/Users/guoyafei/Documents/Vmap3/density")
require(RIdeogram)
gene_density <- read.table("all.count.txt", header=T, stringsAsFactors = F)
wheat_karyotype <- read.table("wheat_karyotype.txt", header=T, stringsAsFactors = F)
ideogram(karyotype = wheat_karyotype, overlaid = gene_density)
convertSVG("chromosome.svg", device = "pdf")
#wheat_karyotype
#Chr Start       End  CE_start    CE_end
#1     0 248956422 122026459 124932724
#2     0 242193529  92188145  94090557
#gene_density
#Chr   Start     End Value
#1       1 1000000    65
#1 1000001 2000000    76

#亚基因组位点数量分布
library(ggplot2)
library(reshape)
library(UpSetR)
setwd("/Users/guoyafei/Desktop/RAiSD/gene")
#part1

expressionInput <- c(
  AABB = 82420294, AABBDDAB = 46757947, AABBDDD = 22849127, DD = 11670269,
  `AABB&AABBDDAB` = 162162813, `AABBDDD&DD` = 18599770)

expressionInput <- c(
  AABB = 82, AABBDDAB = 47, AABBDDD = 23, DD = 12,
  `AABB&AABBDDAB` = 162, `AABBDDD&DD` = 18)

upset(fromExpression(expressionInput), keep.order = TRUE,  text.scale = c(2),point.size = 2.5, line.size = 1.5)

setwd("/Users/guoyafei/Desktop/3type/location")
library(RIdeogram)
library(RColorBrewer)
library(ggplot2)
library(CMplot)
end <- read.table("chrend_21chr.txt",header=F,stringsAsFactors = F)
end$start <- 0
DmChromosome <- end[,c(2,4,3)]
colnames(DmChromosome) <- c("Chr","Start","End")

pos <- read.table("type3.verystrong.21chr.10k.7k.pos", header=F,stringsAsFactors = F)
data <- pos[,c(3,1,2)]
pos$type <- "decisive"
pos$shape <- "box"
pos$color <- "FFFFB3"
pos$name <- "region"
colnames(pos) <- c("Chr","Start","End","Type","Shape","color","name")
gene_list <- pos[,c(4,5,1,2,3,6,7)]
sub <- gene_list[1:10,]
ideogram(karyotype=DmChromosome, #染色体信息
                  #overlaid=Dmdensity, #基因密度
                   label=sub, #目标基因标签
                   label_type="marker",width=135)

####################################################################曼哈顿图(lfmm & baypass)并标注基因###########################################################
library(qqman)
library(tidyverse)
library(ggrepel)
setwd("/Users/guoyafei/Desktop/3type/曼哈顿")

input <- c("baypass.solar1_21chr.txt","lfmm.solar1_21chr.txt")
output <- c("baypass.solar1_21chr.pdf","lfmm.solar1_21chr.pdf")

point <- read.table("select.30k.21chr.bed", header=F,stringsAsFactors = F)
point$SNP <- paste(point$V1,point$V2,sep="-")
snpsOfInterest <- point$SNP
anno <- point[,c(1,2,7)]
anno$P <- 30
colnames(anno) <- c("CHR","BP","SNP","P")

for (i in 1){
  gwasResults2 <- read.table(input[i], header=T, stringsAsFactors = F)
  gwasResults3 <- gwasResults2[sample(nrow(gwasResults2), 25000), ]
  gwasResults <- rbind(gwasResults3,anno)
  don <- gwasResults %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) 

  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  
  data <- subset(don, is_highlight=="yes")
  all <- merge(data,point,by="SNP")
  p <- ggplot(don, aes(x=BPcum, y=P)) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    #scale_color_manual(values = rep(c("grey","skyblue","#E69F00"), 7)) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +   
    #Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    theme_bw() +
    geom_vline(xintercept=all$BPcum) +
    scale_color_manual(values = rep(c("grey","#CC79A7","#A07D35"), 7)) +
    #scale_color_manual(values = c("grey","skyblue","#E69F00")) +
    #geom_point(data=all, aes(x=BPcum, y=P,color=as.factor(V6), size=2)) +
    new_scale_color() +
    geom_label_repel(data=all, aes(label=V5,color=as.factor(V6)), size=2,max.overlaps=20) +
    scale_color_manual(values = c("#009E73","#0072B2","#D55E00")) +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=15),
      axis.text.y=element_text(size=15),
      axis.title.y=element_text(size = 15),
      axis.title.x=element_text(size = 15),
    )+
    xlab("CHR")+
    ylab("Bayes Factor")
  #scale_y_continuous(limits = c(0,7))+
  #geom_point(data=point,aes(x=BPcum,y=-log10(P)),color="red")
  pdf(output[i],height = 2.5,width = 15)
  print(p)
  dev.off()
}

##########################################################位点受选择情况#######################################################
library(ggplot2)
library(reshape)
setwd("/Users/guoyafei/Desktop/3type/选择/")

input <- paste("pop",c(1,2,3,7,9,11,15,16,17,18,20,22,23),".type1.mode.salti.out",sep="")
output <- paste("pop",c(1,2,3,7,9,11,15,16,17,18,20,22,23),".type1.mode.salti2.pdf",sep="")

for(i in c(1:length(input))){
  data <- read.table(input[i], header = F, stringsAsFactors = F)
  data$value <- 0
  data[which(data$V3 == 1),5] <- 1
  data[which(data$V3 > 1),5] <- 2
  data2 <- as.data.frame(cast(data,V2~value))
  data22 <- data2[2:17,c(1,3,4)] 
  data3 <- melt(data22,id="V2",value.name = c("1","2"))
  data3$type <- NA
  data3[which(data3$V2 > 1 & data3$V2 < 8), 4] <- "Class1"
  data3[which(data3$V2 >=8 & data3$V2 < 14), 4] <- "Class2"
  data3[which(data3$V2 >=14), 4] <- "Class3"
  pdf(output[i])
  p <- ggplot(data=data3, mapping=aes(x = type, y = value,fill=variable))+
    geom_bar(stat="identity",position="fill")+
    theme_classic()
  print(p)
  dev.off()
}

#########################################位点BF和beta值的分布情况##############################
setwd("/Users/guoyafei/Desktop/beta/")
library(ggplot2)
type1 <- read.table("elevation.inter.txt", header=F, stringsAsFactors = F)

type1$color <- NA
type1[which(type1$Fold >= 3), 4] <- 1
type1[which(type1$Fold < 3), 4] <- 0
data1 <- type1[which(type1$color ==1),]
data2 <- type1[which(type1$color ==0),]

pdf("elevation.inter.abs.beta.pdf")
ggplot()+
  geom_point(data=type1,aes(x = abs(V8),y = abs(V5),alpha=0.2), color="#000000")+ 
  #geom_point(data=data1,aes(x = BF,y = Beta, alpha=0.2),color="red")+ 
  #scale_fill_manual(values = c("red","blue"))+  
  #labs(x = "Bayes Factor",y = "Beta") +
  labs(x = "lfmm",y = "baypass") +
  #facet_wrap(Type1~Type2,scales="free")+
  theme_classic() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=30),
    axis.text.y=element_text(size=30),
    axis.title.y=element_text(size = 30),
    axis.title.x=element_text(size = 30)
  )+
  theme(legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))
dev.off()

#################################################################### 位点合并片段长度统计 ###########################################################
setwd("/Users/guoyafei/Desktop/baypass/size")
library(ggplot2)
type1 <- read.table("type1.size.txt",header=F,stringsAsFactors = F)  
type1$type <- "type1"
type2 <- read.table("type2.size.txt",header=F,stringsAsFactors = F)  
type2$type <- "type2"
type3 <- read.table("type3.size.txt",header=F,stringsAsFactors = F)  
type3$type <- "type3"
all <- rbind(type1,type2,type3)
sub <- all[which(all$V1 > 3.8),]

ggplot(sub, aes(V1, log10(V2), color = type,shape=type)) +
  geom_point(size=3,alpha=1) +
  scale_color_manual(values = c("#009E73","#56B4E9","#E69F00")) +
  theme_classic()+
  xlab("Region size, bp (log10)")+
  ylab("No. of Features (log10)")+
  #title("Bayenv feature sizes")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

############################################# 位点效应值统计 #######################################
setwd("/Users/guoyafei/Desktop/baypass/size")

data1 <- read.table("type1.baypass.lfmm.10k.bed.region7k.shuf10k", header=F,stringsAsFactors = F)
data1$type <- "type1"
data2 <- read.table("type2.baypass.lfmm.10k.bed.region7k.shuf10k", header=F,stringsAsFactors = F)
data3 <- read.table("type3.baypass.lfmm.10k.bed.region7k.shuf10k", header=F,stringsAsFactors = F)
data2$type <- "type2"
data3$type <- "type3"

all <- rbind(data1,data2,data3)
sub <- all[which(all$V4 !=0),]
a <- lm(sub$V4~sub$V5)
summary(a)

ggplot(sub, aes(V4,V5,color=type)) +
  geom_point(size=3,alpha=0.1) +
  scale_color_manual(values = c("#009E73","#56B4E9","#E69F00")) +
  theme_classic()+
  xlab("bayes factor")+
  ylab("beta")+
  #title("Bayenv feature sizes")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

sub <- all[which(all$V7 !=0),]
b <- -log10(sub$V7)
a <- lm(b~sub$V8)
summary(a)
ggplot(sub, aes(-log10(V7), V8,color=type)) +
  geom_point(size=3,alpha=0.05) +
  scale_color_manual(values = c("#009E73","#56B4E9","#E69F00")) +
  theme_classic()+
  xlab("lfmm -log10(P) value")+
  ylab("beta")+
  #title("Bayenv feature sizes")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

############################################# 三种类型位点等位基因频率变化的分布 #######################################
setwd("/Users/guoyafei/Desktop/baypass/class/freq/")
library(reshape)
library(ggplot2)

data <- read.table("elevation-bottom.env.frq", header=T,stringsAsFactors = F)
data <- read.table("elevation-top.env.frq", header=T,stringsAsFactors = F)

data <- read.table("elevation.20.env.frq", header=T,stringsAsFactors = F)
data <- read.table("elevation.15.env.frq", header=T,stringsAsFactors = F)
data <- read.table("elevation.10.env.frq", header=T,stringsAsFactors = F)

data <- read.table("temp9.20.env.frq", header=T,stringsAsFactors = F)
data <- read.table("temp9.15.env.frq", header=T,stringsAsFactors = F)
data <- read.table("temp9.10.env.frq", header=T,stringsAsFactors = F)

data <- read.table("type3.decisive.env.frq", header=T,stringsAsFactors = F)
data <- read.table("type3.verystrong.env.frq", header=T,stringsAsFactors = F)
data <- read.table("type3.strong.env.frq", header=T,stringsAsFactors = F)


data <- read.table("type1.decisive.temp9.env.frq", header=T,stringsAsFactors = F)
data <- read.table("type1.verystrong.temp9.env.frq", header=T,stringsAsFactors = F)
data <- read.table("type1.strong.temp9.env.frq", header=T,stringsAsFactors = F)

data <- read.table("type2.decisive.solar1.env.frq", header=T,stringsAsFactors = F)
data <- read.table("type2.verystrong.solar1.env.frq", header=T,stringsAsFactors = F)
data <- read.table("type2.strong.solar1.env.frq", header=T,stringsAsFactors = F)

data <- read.table("type3.decisive.soil12.env.frq", header=T,stringsAsFactors = F)
data <- read.table("type3.verystrong.soil12.env.frq", header=T,stringsAsFactors = F)
data <- read.table("type3.strong.soil12.env.frq", header=T,stringsAsFactors = F)
#temp1
#sub <- data[,c(1,2,3,29,40:1039)]
#temp5
#sub <- data[,c(1,2,3,35,40:1039)]
#temp10
#sub <- data[,c(1,2,3,30,40:1039)]
#temp9
sub <- data[,c(1,2,3,39,40:1039)]
#elevation
sub <- data[,c(1,2,3,4,40:1039)]
#solar1
#sub <- data[,c(1,2,3,26,40:1039)]
#prec7
#sub <- data[,c(1,2,3,11,40:1039)]
#prec6
#sub <- data[,c(1,2,3,10,40:1039)]
#soil12
#sub <- data[,c(1,2,3,16,40:1039)]
#sub2 <- melt(sub, id=c("pop","X","Y","bio_temp1"))
sub2 <- melt(sub, id=c("pop","X","Y","bio_temp9"))
sub2 <- melt(sub, id=c("pop","X","Y","bio_Elevation"))
#sub2 <- melt(sub, id=c("pop","X","Y","bio_solar1"))
#sub2 <- melt(sub, id=c("pop","X","Y","bio_prec6"))
#sub2 <- melt(sub, id=c("pop","X","Y","bio_soil12"))
sub3 <- sub2[c(3000:4500),]
#sub3 <- sub2[c(501:1000),]
ggplot(sub3, aes(x=bio_Elevation, y=value, color = variable)) +
  geom_point(size=0.5)+
  geom_line()+
  theme_classic()+
  ylab("Allele Frequency")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

############################################# 三种类型位点受选择以及lassip的结果分析 #######################################
setwd("/Users/guoyafei/Desktop/lassip")
library(ggplot2)
data <- read.table("withoutcentromere.top001.envlassip.5pop.txt",header=T,sep="\t",stringsAsFactors = F)
data <- read.table("withoutcentromere.top01.envlassip.txt",header=T,sep="\t",stringsAsFactors = F)
data <- read.table("withoutcentromere.top005.envlassip.test.nocoding.txt",header=T,sep="\t",stringsAsFactors = F)
data <- read.table("5pop_envlassip.30.6.3type.txt",header=T,stringsAsFactors = F)
data[which(data$pop1_m > 2), 6] <- 2
ggplot(data=data, mapping=aes(x = pop1_m,fill=type))+
  facet_grid(name~.)+
  geom_bar()+
  theme_classic()

library(reshape)
data2 <- data[,c(9,10,6,11)]
data2$select <- NA
data2[which(data2$pop1_m > 1), 5] <- "soft"
data2[which(data2$pop1_m == 1), 5] <- "hard"
data2[which(data2$pop1_m < 1), 5] <- "none"
data2 <- data2[which(data2$select != "none"),]
data2 <- data2[which(data2$pop != "AF" & data2$pop != "AM" ),]

data2 <- data2[which(data2$name != "AF" & data2$name != "AM" ),]
data2$value <- 1

data3 <- cast(data2,name+type+pop~select)
data3 <- cast(data2,name+type~select)
data4 <- melt(data3,id=c("name"))
data4$type <- factor(data4$type,levels = c("decisive","verystrong","strong"))
data4$type <- factor(data4$type,levels = c("top30","mid30","bottom30"))
data4$pop <- factor(data4$pop,levels = c("EU","WA","IA","SH","EA"))
data4$name <- factor(data4$name,levels = c("EU","WA","IA","SH","EA"))

data3 <- cast(data2,pop~select)
data4 <- melt(data3,id=c("pop"))
data4$pop <- factor(data4$pop,levels = c("EU","WA","IA","SH","EA"))

data4$type <- "non-coding"
all <- data4
out <- rbind(all,data4)
out$type2 <- paste(out$select,out$type,sep="-")
out$type2 <- factor(out$type2,levels = c("hard-non-coding","hard-coding" ,"soft-non-coding","soft-coding"))

ggplot(data=data4, mapping=aes(x = name, y=value, fill=select))+
  facet_grid(pop~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  #theme(strip.text.x = element_text(size = 30, colour = "#FC0043")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.text.y = element_text(size = 15)) +
  #theme(strip.background.x = element_rect(fill = "#00FCE6", colour = "#00FCE6")) +
  theme(strip.background.y = element_rect(fill = "#FC4E07")) +
  #theme(strip.placement = "outside") + # 分面条带放外面
  theme(strip.switch.pad.grid = unit(1, "inch")) + # 设置分面条带与坐标轴的距离
  ylab("Counts")

ggplot(data=out, mapping=aes(x = pop, y=value, fill=type2))+
  #facet_grid(type~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="fill")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  ylab("Frequency")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
  
ggplot(data=data4, mapping=aes(x = type, y=value, fill=select))+
  facet_grid(pop~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  #theme(strip.text.x = element_text(size = 30, colour = "#FC0043")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.text.y = element_text(size = 15)) +
  #theme(strip.background.x = element_rect(fill = "#00FCE6", colour = "#00FCE6")) +
  theme(strip.background.y = element_rect(fill = "#FC4E07")) +
  #theme(strip.placement = "outside") + # 分面条带放外面
  theme(strip.switch.pad.grid = unit(1, "inch")) + # 设置分面条带与坐标轴的距离
  ylab("Counts")

ggplot(data=data4, mapping=aes(x = name, y=value, fill=select))+
  facet_grid(type~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  #theme(strip.text.x = element_text(size = 30, colour = "#FC0043")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.text.y = element_text(size = 15)) +
  #theme(strip.background.x = element_rect(fill = "#00FCE6", colour = "#00FCE6")) +
  theme(strip.background.y = element_rect(fill = "#FC4E07")) +
  #theme(strip.placement = "outside") + # 分面条带放外面
  theme(strip.switch.pad.grid = unit(1, "inch")) + # 设置分面条带与坐标轴的距离
  ylab("Counts")

ggplot(data=data4, mapping=aes(x = pop, y=value, fill=select))+
  facet_grid(type~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  ylab("Frequency")+ 
  theme(strip.text.y = element_text(size = 15)) +
  #theme(strip.background.x = element_rect(fill = "#00FCE6", colour = "#00FCE6")) +
  theme(strip.background.y = element_rect(fill = "#FC4E07")) +
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 25))

################################################### 位点的BF和GWAS beta的相关性分析 ##############################################
setwd("/Users/guoyafei/Desktop/baypass/class/")
library(ggplot2)
data <- read.table("elevation_BF_beta_shuf2k.txt",header=T,sep="\t",stringsAsFactors = F)

data <- read.table("temp1.baypass.lfmm.GF.GWAS.both.bed",header=T,sep="\t",stringsAsFactors = F)
ggplot(data, aes(x=abs(data$lfmm_beta), y=abs(data$gwas_BETA))) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  #xlab("BETA")+
  #ylab("GF")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none", plot.title = element_text(color="red", size=20, face="bold.italic"), legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

ggplot(data, aes(x=BF, y=abs(BETA))) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("Elevation Bayes Factor")+
  ylab("Elevation GWAS BETA")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

################################################## 位点的BF和GWAS beta和群体间fst的相关性分析 ##############################################
setwd("/Users/guoyafei/Desktop/baypass/class/fst")
data <- read.table("elevation_BF_beta_fst_top1_max_rmAFAM.txt", header = F,stringsAsFactors = F)
data <- read.table("temp1.baypass.lfmm.GF.GWAS.lfmm.fst.merge.bed", header = T,stringsAsFactors = F)
sub <- data[sample(1:dim(data)[1],5000),]
sub <- data[which(data$GF_rsq !=0),]
sub2 <- sub[which(sub$gwas_P < 0.05),]
ggplot(sub2, aes(x=effect, y=gwas_BETA)) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("Elevation GWAS BETA")+
  ylab("Population Mean Fst")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

sub <- data[which(data$effect != 0),]
sub2 <- sub[which(sub$gwas_P < 0.05),]
ggplot(sub2, aes(x=BF, y=abs(effect))) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("Elevation GWAS BETA")+
  ylab("Population Mean Fst")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

ggplot(sub2, aes(x=-log(lfmm_P), y=abs(lfmm_beta))) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("Elevation GWAS BETA")+
  ylab("Population Mean Fst")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

sub2 <- sub[which(sub$GF_rsq != 0),]
ggplot(sub2, aes(x=GF_rsq, y=abs(lfmm_beta))) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("GF_rsq")+
  ylab("gwas_BETA")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

ggplot(sub2, aes(x=max.HUDSON_FST., y=abs(gwas_BETA))) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("Elevation GWAS BETA")+
  ylab("Population Mean Fst")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

ggplot(data, aes(x=abs(BETA), y=median.HUDSON_FST.)) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("Elevation GWAS BETA")+
  ylab("Population Median Fst")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

ggplot(data, aes(x=abs(BETA), y=min.HUDSON_FST.)) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("Elevation GWAS BETA")+
  ylab("Population Min Fst")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))

ggplot(data, aes(x=abs(BETA), y=max.HUDSON_FST.)) +
  geom_point(size=0.5)+
  #geom_line()+
  theme_classic()+
  xlab("Elevation GWAS BETA")+
  ylab("Population Max Fst")+
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))


################################################## 适应性位点，选择分析，基因联合分析 ##############################################
library(ggplot2)
library(reshape)
library(UpSetR)
setwd("/Users/guoyafei/Desktop/RAiSD/gene")
#part1

data <- read.table("withoutcentromere.top05.envlassip.addgene.txt",header=T,sep="\t",stringsAsFactors = F)
data[which(data$pop1_m > 2), 6] <- 2
ggplot(data=data, mapping=aes(x = pop1_m,fill=type))+
  facet_grid(name~.)+
  geom_bar()+
  theme_classic()
data2 <- data[,c(9,10,6,11,15)]
data2$select <- NA
data2[which(data2$pop1_m > 1), 6] <- "soft"
data2[which(data2$pop1_m == 1), 6] <- "hard"
data2[which(data2$pop1_m < 1), 6] <- "none"

data2 <- data2[which(data2$select != "none"),]
data2 <- data2[which(data2$pop != "AF" & data2$pop != "AM" ),]

EA <- data2[which(data2$pop == "EA"),5,drop=F]
EA2 <-  EA[!duplicated(EA, fromLast=TRUE), ] 
EU <- data2[which(data2$pop == "EU"),5,drop=F]
EU2 <-  EU[!duplicated(EU, fromLast=TRUE), ] 
WA <- data2[which(data2$pop == "WA"),5,drop=F]
WA2 <-  WA[!duplicated(WA, fromLast=TRUE), ] 
SH <- data2[which(data2$pop == "SH"),5,drop=F]
SH2 <-  SH[!duplicated(SH, fromLast=TRUE), ] 
IA <- data2[which(data2$pop == "IA"),5,drop=F]
IA2 <-  IA[!duplicated(IA, fromLast=TRUE), ] 

a <- list(EA<-EA2,EU<-EU2,WA<-WA2,SH<-SH2,IA<-IA2)
names(a) <- c("EA","EU","WA","SH","IA")
tmp <- names(a)
p <- upset(fromList(a), keep.order = TRUE, sets=tmp, text.scale = c(2),point.size = 2.5, line.size = 1.5)

#part2
data1 <- read.table("withoutcentromere.top005.envlassip.addgene.uniq.gene.out.merge.txt",header=F,sep="\t",stringsAsFactors = F)
data2 <- read.table("withoutcentromere.top005.envlassip.addgene.repeat.gene.out.merge.txt",header=F,sep="\t",stringsAsFactors = F)
data1$type <- "uniq"
data2$type <- "repeat"
all <- rbind(data1,data2)
all$select <- NA
all[which(all$V1 > 2), 6] <- "soft"
all[which(all$V1 == 1), 6] <- "hard"
all[which(all$V1 == 0), 6] <- "none"
all2 <- all[which(all$select != "none"),]
colnames(all2) <- c("M","T","pop","gene","type","select")
all2$value <- 1
data3 <- cast(all2,pop+type+select~value)



data3$pop <- factor(data3$pop,levels = c("EU","WA","IA","SH","EA"))
colnames(data3)[4] <- "value"

ggplot(data=data3, mapping=aes(x = pop, y=value, fill=select))+
  facet_grid(type~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="fill")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  #theme(strip.text.x = element_text(size = 30, colour = "#FC0043")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.text.y = element_text(size = 15)) +
  #theme(strip.background.x = element_rect(fill = "#00FCE6", colour = "#00FCE6")) +
  theme(strip.background.y = element_rect(fill = "#FC4E07")) +
  #theme(strip.placement = "outside") + # 分面条带放外面
  theme(strip.switch.pad.grid = unit(1, "inch")) + # 设置分面条带与坐标轴的距离
  ylab("Counts")


ggplot(data=data4, mapping=aes(x = pop, y=value, fill=select))+
  #facet_grid(pop~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="fill")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  ylab("Frequency")+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))


ggplot(data=data4, mapping=aes(x = type, y=value, fill=select))+
  facet_grid(pop~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  #theme(strip.text.x = element_text(size = 30, colour = "#FC0043")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.text.y = element_text(size = 15)) +
  #theme(strip.background.x = element_rect(fill = "#00FCE6", colour = "#00FCE6")) +
  theme(strip.background.y = element_rect(fill = "#FC4E07")) +
  #theme(strip.placement = "outside") + # 分面条带放外面
  theme(strip.switch.pad.grid = unit(1, "inch")) + # 设置分面条带与坐标轴的距离
  ylab("Counts")

ggplot(data=data4, mapping=aes(x = name, y=value, fill=select))+
  facet_grid(type~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  #theme(strip.text.x = element_text(size = 30, colour = "#FC0043")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.text.y = element_text(size = 15)) +
  #theme(strip.background.x = element_rect(fill = "#00FCE6", colour = "#00FCE6")) +
  theme(strip.background.y = element_rect(fill = "#FC4E07")) +
  #theme(strip.placement = "outside") + # 分面条带放外面
  theme(strip.switch.pad.grid = unit(1, "inch")) + # 设置分面条带与坐标轴的距离
  ylab("Counts")

ggplot(data=data4, mapping=aes(x = pop, y=value, fill=select))+
  facet_grid(type~.)+
  #geom_bar(position="dodge")+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  scale_fill_manual(values =c("#00AFBB", "#E7B800","#4682B4","#FC4E07" ))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_discrete(expand = c(0, 0))+
  ylab("Frequency")+ 
  theme(strip.text.y = element_text(size = 15)) +
  #theme(strip.background.x = element_rect(fill = "#00FCE6", colour = "#00FCE6")) +
  theme(strip.background.y = element_rect(fill = "#FC4E07")) +
  theme(plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=20),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 25))

### ibs distance ####
setwd("/Users/guoyafei/Documents/Vmap3/ibs")
data <- read.table("AB_CS.ibs.txt", header=T,stringsAsFactors = F)
a <- aggregate(data, by=list(data$type1),FUN=mean)[,c(1,3)]
b <- aggregate(data, by=list(data$type1),FUN=sd)[,c(1,3)]
all1 <- cbind(a,b)
colnames(all1) <- c("type","mean","type2","sd")

data <- read.table("D_CS.ibs.txt", header=T,stringsAsFactors = F)
c <- aggregate(data, by=list(data$type1),FUN=mean)[,c(1,3)]
d <- aggregate(data, by=list(data$type1),FUN=sd)[,c(1,3)]
all2 <- cbind(c,d)
colnames(all2) <- c("type","mean","type2","sd")

all1$type2 <- "AB"
all2$type2 <- "D"
all <- rbind(all1,all2)


ggplot(data=all, mapping=aes(x = type, y = mean, group= type2, fill=type2)) +
  #scale_fill_manual(values = c("#A07D35","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00"))+ 
  geom_bar(stat="identity",position=position_dodge(0.85),width=0.75) +
  geom_errorbar(aes(ymax = mean+sd, ymin = mean-sd),
                position = position_dodge(0.85), width = 0.1)+
  #facet_grid(tissue~.)+
  theme_classic()

AB <- read.table("AB_CS.ibs.txt", header=T,stringsAsFactors = F)
AB$type3 <- "AB"

D <- read.table("D_CS.ibs.txt", header=T,stringsAsFactors = F)
D$type3 <- "D"
all <- rbind(AB,D)
all$type1 <- factor(all$type1,levels = c("Wild_emmer", "Domesticated_emmer","Ispahanicum","Georgian_wheat","Durum","Persian_wheat","Polish_wheat","Khorasan_wheat","Rivet_wheat","Spelt","Macha","Xinjiang_wheat","Vavilovii","Tibetan_semi_wild","Club_wheat","Indian_dwarf_wheat","Yunan_wheat","Landrace","Cultivar","OtherHexaploids","Strangulata"))
all$type2 <- factor(all$type2,levels = c("Wild_emmer", "Domesticated_emmer","Free_threshing_tetraploids","OtherTetraploids","Landrace","Cultivar","OtherHexaploids","Strangulata"))
ggplot(all, mapping=aes(x = type2, y = ABD_1143,fill=type3))+
  geom_jitter(aes(color=type3),shape=16,alpha = 0.3, position = position_jitter(0.2))+
  geom_boxplot() +
  #geom_dotplot(binaxis = "y", stackdir = "center",
  #             position = position_dodge(1),dotsize=0.02)+
  theme_classic() +  
  theme(
    #legend.position="none",
    #panel.border = element_blank(),
    #axis.line.y = element_line(),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(angle = 45,hjust = 1,size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
  )




### 每个环境鉴定到的基因上面显著位点的分布 ####
setwd("/Users/guoyafei/Desktop/baypass/Num")
data <- read.table("type1_Num.txt", header=T, stringsAsFactors = F)
ggplot(data, aes(x = CDS))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE) +
  geom_density()+
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
