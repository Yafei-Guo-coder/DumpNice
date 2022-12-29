#h12----
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
library(reshape2)
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/H12")
data <- read.table("h12.selection.pro.txt", header=F,stringsAsFactors = F)
data$V1 <- factor(data$V1,levels = c("hap25","hap50","hap75","hap100"))
data2 <- as.data.frame(data[which(data$V2=="landrace" | data$V2=="cultivar" ),])
#整体柱状图----
ggplot(data2, aes(x = V1,y=V3))+
  facet_grid(.~V2)+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_bar(stat = 'identity') +
  theme_bw() +
  ylim(0,1) +
  scale_color_manual(values = brewer.pal(12, "Paired"))+ 
  xlab("H12 haplotype") +
  ylab("Hard sweep proportion")

#散点图----
data <- read.table("flower.selection.txt", header=T,stringsAsFactors = F)
data$hap <- factor(data$hap,levels = c("hap25","hap50","hap75","hap100"))
data2 <- as.data.frame(data[which(data$hap=="hap75"),])
wildemmer <- data2[which(data2$name == "wildemmer"),4]
domemmer <- data2[which(data2$name == "domemmer"),4]
freethresh <- data2[which(data2$name == "freethresh"),4]
landrace <- data2[which(data2$name == "landrace"),4]
cultivar <- data2[which(data2$name == "cultivar"),4]
strangulata <- data2[which(data2$name == "strangulata"),4]

library(UpSetR)
listInput <- list(wildemmer = wildemmer, domemmer = domemmer, freethresh = freethresh, landrace=landrace,cultivar=cultivar,strangulata=strangulata)
upset(fromList(listInput), keep.order = TRUE,sets = c("strangulata","cultivar","landrace","freethresh","domemmer","wildemmer"), order.by = "freq",text.scale = c(2),point.size = 2.5, line.size = 1.5)

cast <- dcast(data2,hap+name~type)
a <- cast[,1:3]
a$type <- "hard"
colnames(a)[1:3] <- c("hap","name","count")
b<- cast[,c(1,2,4)]
b$type <- "soft"
colnames(b)[1:3] <- c("hap","name","count")
all <- rbind(a,b)

all$name <- factor(all$name,levels = c("wildemmer","domemmer","freethresh","landrace","cultivar","strangulata"))

ggplot(all,mapping=aes(x = name, y = count,fill=type))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_bar(stat="identity",position="stack")+
  theme_bw()+
  scale_fill_brewer(palette="Set2")

ggplot(cast, aes(x = wildemmer,y = domemmer))+
  geom_point()+ 
  labs(x = "Taxa Number",y = "Snp proportion") +
  theme_bw()


################################################ hard versus soft figure 4 #############################################
library(ggplot2)
#library(gridExtra)
library(RColorBrewer)
#library(fitdistrplus)
library(ggrepel)

setwd("/Users/guoyafei/Documents/02_VmapIII/18_h12/hap100")
name <- read.table("/Users/guoyafei/Documents/02_VmapIII/11_piratio/pi-ratio/gene.txt", header=F, stringsAsFactors = F)
rownames(name) <- name$V1

#thresh
#hap75
#cultivar_t <- 0.06642
#landrace_t <- 0.06306

#hap100
cultivar_t <- 0.041964165079476
landrace_t <- 0.033131204459876

#landrace
all1 <- read.table("landrace.AB.h12.shuf5k", header=F,stringsAsFactors = F)
gene1 <- read.table("landrace.gene.AB.h12",header=F,stringsAsFactors = F)
gene1$name <- name[gene1$V1,2]
gene1$type <- "NA"
gene1[which(gene1$V2 > landrace_t),5] <- "yes"
sub_gene1 <- gene1[which(gene1$type  == "yes"),]

#cultivar
all1 <- read.table("cultivar.AB.h12.shuf5k", header=F,stringsAsFactors = F)
gene1 <- read.table("cultivar.gene.AB.h12",header=F,stringsAsFactors = F)
gene1$name <- name[gene1$V1,2]
gene1$type <- "NA"
gene1[which(gene1$V2 > landrace_t),5] <- "yes"
sub_gene1 <- gene1[which(gene1$type  == "yes"),]


ggplot(all1, aes(V2, V3)) +
  geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
  geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_vline(xintercept=landrace_t,color='red',linetype = "dotted")+
  geom_hline(yintercept=0.059,color='red',linetype = "dotted")+
  geom_text_repel(data = sub_gene1,aes(V2, V3, label = name),max.overlaps = 100) +
  xlab("H12")+
  ylab("H2/H1")+
  theme_bw()+
  theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

################################# consensus haplotype frequency ##################################
library(fitdistrplus)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
setwd("/Users/guoyafei/Documents/02_VmapIII/18_h12/")
data <- read.table("haplotype.compare.sort.common.txt", header=T, stringsAsFactors = F, sep="\t")
all <- read.table("haplotype.compare.shuf5k.txt", header=T, stringsAsFactors = F, sep="\t")

path <- "/Users/guoyafei/Documents/02_VmapIII/22_xpclr/module" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data1 <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep="\t")})

for (i in c(1:length(data1))) {
  name <- as.data.frame(data1[[i]])
  rownames(name) <- name$V1
  sub <- data[which(data$gene %in% name$V1),]
  out <- paste(names(data1)[i],"common.pdf",sep=".")
  line <- data$landrace_num/sqrt(data$landrace_num+data$cultivar_num)
  fit <- fitdist(line, "norm")
  a <- summary(fit)
  mean <- mean(sub$landrace_num/sqrt(sub$landrace_num+sub$cultivar_num))
  z <- round(abs(mean-a$estimate[[1]])/a$estimate[[2]],2)
  pdf(out,width=6,height = 6)
  p <- ggplot(all, aes(landrace_num, cultivar_num)) +
    geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
    geom_point(data = sub, aes(landrace_num, cultivar_num, size=2, alpha = 0.5,colour = "red")) +
    geom_abline(intercept=0,slope=1,color='red',linetype = "dotted")+
    xlim(0,1)+
    ylim(0,1)+
    xlab("landrace_freq")+
    ylab("cultivar_freq")+
    theme_bw()+
    ggtitle(paste("z-score =",z))+
    theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  print(p)
  dev.off()
}

