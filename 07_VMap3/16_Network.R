setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/shufgene/New")
library(ggplot2)
library(gridExtra)
#filenames <- c("shuf.127.txt","shuf.176.txt","shuf.245.txt","shuf.279.txt","shuf.287.txt","shuf.471.txt","shuf.712.txt","shuf.78.txt")
#outnames <- c("shuf.127.pdf","shuf.176.pdf","shuf.245.pdf","shuf.279.pdf","shuf.287.pdf","shuf.471.pdf","shuf.712.pdf","shuf.78.pdf")
filenames <- c("shuf.127.txt","shuf.176.txt","shuf.279.txt","shuf.78.txt")
outnames <- c("shuf.127.pdf","shuf.176.pdf","shuf.279.pdf","shuf.78.pdf")

data <- read.table(filenames[1], header=T,stringsAsFactors = F)
p <- list()
p[[1]] <- ggplot(data, aes(x=genenum))+geom_histogram(bins=5,aes(y=..density..))+xlab("Gene number")+theme_bw()
p[[2]] <- ggplot(data, aes(x=degree))+geom_histogram(bins=50,aes(y=..density..))+xlab("Degree")+theme_bw()
p[[3]] <- ggplot(data, aes(x=sumR))+geom_histogram(bins=50,aes(y=..density..))+xlab("sumR")+theme_bw()
pdf(outnames[1],width = 10.7,height = 3.7)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=1)
dev.off()

data1 <- read.table(filenames[1], header=T,stringsAsFactors = F)
data2 <- read.table(filenames[2], header=T,stringsAsFactors = F)
data3 <- read.table(filenames[3], header=T,stringsAsFactors = F)
data4 <- read.table(filenames[4], header=T,stringsAsFactors = F)

######
data1 <- read.table("shuf.217.degree",header=F,stringsAsFactors = F)
ggplot(data1, aes(x=V1))+geom_histogram(bins=15,aes(y=..density..))+xlab("Degree")+theme_bw()


data2 <- read.table("shuf.217.degree.sum",header=F,stringsAsFactors = F)
ggplot(data2, aes(x=V2))+geom_histogram(bins=15,aes(y=..density..))+xlab("sumR")+theme_bw()


pdf(outnames[1],width = 10.7,height = 3.7)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=1)
dev.off()


#PPI
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
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












