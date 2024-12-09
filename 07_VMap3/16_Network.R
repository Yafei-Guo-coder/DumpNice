setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/")
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(fitdistrplus)
display.brewer.all()

#All R2 distribution-----
setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/all")
path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/all/shuf5k" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=",")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  p[[i]] <- ggplot(data=name, aes(V4))+geom_histogram(binwidth = 0.01,aes(y=..density..))+
    theme_bw()+
    #geom_vline(xintercept=flower_mean[i],color='red')+
    xlab("R2")+
    ggtitle(names[i])
}
pdf("all_R2.pdf",width = 10,height = 7)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()

#expression Module------
names <- c("S1coleoptile","S1root","S2leaf","S2root","S3leaf","S4Awn","S4leaf","S4spike","S4stem","S5Anthers","S5leaf","S6grain","S6leaf","S7grain","S7leaf","S8grain","S8leaf","S9grain")
#r03
#seed_degree <- c(1048,640,482,609,381,534,513,1595,779,545,361,666,381,3463,468,1450,954,707)
#seed_mean <- c(0.38746660305344,0.38995875,0.38906286307054,0.37375697865353,0.37878005249344,0.38261235955056,0.38886900584795,0.39889937304075,0.37766379974326,0.40305119266055,0.39592936288089,0.3913990990991,0.40296719160105,0.42257938203869,0.39139102564103,0.39571020689655,0.40205377358491,0.39270311173975)
#flower_degree <- c(3587,2322,2605,2783,1464,1840,2069,6453,2792,2230,1309,2441,1529,15982,1597,6682,2685,4492)
#flower_mean <- c(0.38503788681349,0.37501042204996,0.41120671785029,0.38280783327345,0.3798131147541,0.37531902173913,0.39039574673755,0.39073706803037,0.37848896848138,0.41091049327354,0.37982123758594,0.38136239246211,0.39326108567691,0.42212898886247,0.38440682529743,0.40172698293924,0.39550283054004,0.40384808548531)

#r01
seed_mean <- c(0.19360946752223,0.17978099141796,0.1818948787062,0.18066811814149,0.17666317490494,0.1817860858465,0.18148121623976,0.2070258840028,0.18960978580991,0.18521523826458,0.17388497119916,0.18560951145038,0.17658887775551,0.24460937871581,0.1782179683512,0.2025521420301,0.19895087968218,0.18338525582007)
seed_degree <- c(8883,7807,5565,8566,5260,6267,5739,9983,7470,5624,5729,6550,4990,13456,5877,10364,7048,7603)
flower_degree <- c(37369,32774,26497,39885,22516,25965,24015,44613,31396,23569,23470,29696,21666,63902,23592,45579,27011,36798)
flower_mean <- c(0.18823638845032,0.17687784829438,0.18615279088199,0.17914648113326,0.17395727926808,0.17524488349702,0.17882783676869,0.20190767489297,0.18394685947254,0.18303373923374,0.17029267575628,0.18018855401401,0.17445518785193,0.24182734656192,0.17366643777552,0.20645209636017,0.18651013290881,0.19476062557748)

#PPI
mean <- c(0.19157261572783,0.16553474853037,0.19440211002261,0.1828303451582,0.16642771381579,0.18085498676081,0.17802098939929,0.21247850241546,0.18794293369056,0.18016781954887,0.17991066846725,0.18150084210526,0.17022178861789,0.25068429130304,0.1708906274206,0.20651469754253,0.1891889375,0.18416850728155)
degree <- c(1793,1531,1327,2086,1216,1133,1415,2484,1493,1330,1481,1425,1230,3323,1291,2116,1600,1648)

#flower------
setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/shufflower/")
path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/shufflower/mean" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 0.001,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=flower_mean[i],color='red')+
    xlab("mean R2")+
    ggtitle(names[i])
}
pdf("flower_mean.pdf",width = 10,height = 7)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()

path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/shufflower/degree" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 40,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=flower_degree[i],color='red')+
    xlab("degree")+
    ggtitle(names[i])
}

pdf("flower_degree.pdf",width = 10,height = 7)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()

#seed------
setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/shufseed")
path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/shufseed/mean" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 0.001,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=seed_mean[i],color='red')+
    xlab("mean R2")+
    ggtitle(names[i])
}
pdf("seed_mean.pdf",width = 10,height = 7)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()

path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/shufseed/degree" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 20,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=seed_degree[i],color='red')+
    xlab("degree")+
    ggtitle(names[i])
}
pdf("seed_degree.pdf",width = 10,height = 7)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()

#PPI
brewer.pal(8, "Set1")
setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/PPI/")
data <- read.table("count.txt",header=T,stringsAsFactors = F)
a <- data[,c(1:2)]
a$type <- "Score>0"
colnames(a) <- c("proteinN", "count", "type")
b <- data[,c(1,3)]
b$type <- "Score>0.5"
colnames(b) <- c("proteinN", "count", "type")
c <- rbind(a,b)
ggplot(c, aes(x=proteinN, y=count,fill=type)) + geom_bar(stat="identity") +theme_classic()+xlim(0,50)+
  scale_fill_manual(values = c("#FC8D62","#8DA0CB","#E78AC3","#FFD92F"))+
  theme(legend.text = element_text(size=10),legend.title=element_blank(),legend.position=c(0.8,0.8))

data <- read.table("gene-OGs-count.txt", header=F,stringsAsFactors = F)
data[1,1] <- 1500
data[5,1] <- 1200
ggplot(data, aes(x=V2, y=V1, fill=V3)) + geom_bar(stat="identity") + theme_classic() + xlab("OGs counts") + ylab("Gene Number")+
  geom_text(aes(label=V1), position = position_dodge2(width = 0.9, preserve = 'single'), vjust = -0.2, hjust = 0.5) +
  scale_fill_manual(values = c("#FC8D62","#8DA0CB","#E78AC3","#FFD92F"))+
  theme(legend.text = element_text(size=10),legend.title=element_blank(),legend.position=c(0.8,0.8))

data <- read.table("Network.node.0.5cfms.degree.txt", header=F,stringsAsFactors = F)
ggplot(data, aes(x=V1)) + geom_density(aes(y = ..density..),color = "#FC8D62",fill="#FC8D62") + 
  theme_classic() + xlab("Degree") + ylab("Density")


#ppi_expr------
setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/PPI/degree30shuf")
path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/PPI/degree30shuf/mean" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 0.001,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=mean[i],color='red')+
    xlab("mean R2")+
    ggtitle(names[i])
}
pdf("seed_mean.pdf",width = 10,height = 7)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()

path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/PPI/degree30shuf/degree" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 20,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=degree[i],color='red')+
    xlab("degree")+
    ggtitle(names[i])
}
pdf("seed_degree.pdf",width = 10,height = 7)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()




#Correlation strength
names <- c("S1coleoptile","S1root","S2leaf","S2root","S3leaf","S4Awn","S4leaf","S4spike","S4stem","S5Anthers","S5leaf","S6grain","S6leaf","S7grain","S7leaf","S8grain","S8leaf","S9grain")
flower <- c(7034.21,5796.99,4932.49,7145.26,3916.82,4550.23,4294.55,9007.71,5775.2,4313.92,3996.77,5350.88,3779.75,15453.3,4097.14,9409.88,5037.83,7166.8)
seed <- c(1719.83,1403.55,1012.25,1547.6,929.248,1139.25,1041.52,2066.74,1416.39,1041.65,996.187,1215.74,881.178,3291.46,1047.39,2099.25,1402.21,1394.28)
PPI <- c(343.49,253.434,257.972,381.384,202.376,204.909,251.9,527.797,280.599,239.623,266.448,258.639,209.373,833.024,220.62,436.985,302.702,303.51)
#flower------
setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/out")
path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/out/shufflower" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  fit <- fitdist(data[[i]]$V1, "norm")
  a <- summary(fit)
  z <- round(abs(flower[i]-a$estimate[[1]])/a$estimate[[2]],2)
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 50,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=flower[i],color='red')+
    xlab("Correlation strength")+
    ggtitle(paste(names[i],"z-score =",z))
}
pdf("flower.pdf",width = 15,height = 10)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()

path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/r01/out/shufseed" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  fit <- fitdist(data[[i]]$V1, "norm")
  a <- summary(fit)
  z <- round(abs(seed[i]-a$estimate[[1]])/a$estimate[[2]],2)
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 30,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=seed[i],color='red')+
    xlab("Correlation strength")+
    ggtitle(paste(names[i],"z-score =",z))
}

pdf("seedsize.pdf",width = 15,height = 10)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()

#ppi_expr------
setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/PPI/")
path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/PPI/out" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F,sep=" ")})
p <-list ()
for(i in c(1:length(data))){
  name <- as.data.frame(data[[i]])
  fit <- fitdist(data[[i]]$V1, "norm")
  a <- summary(fit)
  z <- round(abs(PPI[i]-a$estimate[[1]])/a$estimate[[2]],2)
  p[[i]] <- ggplot(data=name, aes(V1))+geom_histogram(binwidth = 5,aes(y=..density..))+
    theme_bw()+
    geom_vline(xintercept=PPI[i],color='red')+
    xlab("Correlation strength")+
    ggtitle(paste(names[i],"z-score =",z))
}
pdf("PPI.pdf",width = 15,height = 10)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],nrow=4)
dev.off()


#### phyperTest co-expression network#####
setwd("/Users/guoyafei/Desktop/network/2018science/phyperTest")
shuf_integra <- read.table("shuf/integra_abiotic.txt", header=T, stringsAsFactors = F)
shuf_phyper <- read.table("shuf/phyper_abiotic.txt", header=T, stringsAsFactors = F)

#### z-score Test PPI network#####
setwd("/Users/guoyafei/Desktop/network/PPI")
data <- read.table("z.score.txt", header=T, stringsAsFactors = F)
#(215-200.52)/3.63388
#3.984722
#(202-175.98)/3.41317
#7.623412
#(239-229.57)/3.76148
#2.506992
#(200-179.08)/3.69226
#5.665907
#(214-187.37)/3.69606
#7.20497

data$type <- rownames(data)
pdf("enrich.pdf",width=4, height=2)
ggplot(data, aes(y = type, x = mean)) +
  geom_point(aes(y = type, x = point)) + 
  geom_errorbar(aes(xmin=(mean-sd),xmax=(mean+sd)),width=1,size=1) +
  scale_fill_manual(values = c("red","blue")) +  
  labs(x = "",y = "") +
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
dev.off()


