library(ggplot2)
genome <- read.table("/Users/guoyafei/公共资源/01_数据资源/chrConvertionRule.txt", header=T,stringsAsFactors = F)
ge <- genome[c(2:43),c(1,3)]

sum <- read.table("/Users/guoyafei/Desktop/Va.txt", header=F,stringsAsFactors = F)
mean <- read.table("/Users/guoyafei/Desktop/mean2.txt", header=F,stringsAsFactors = F)
snpnum <- read.table("/Users/guoyafei/Desktop/snpnum.txt", header=F,stringsAsFactors = F)

all <- cbind(ge,sum)[,c(2,4)]

colnames(all) <- c("ID","length","sum","mean","snpnum")
all$V3 <- all$V2/10000000

a <- summary(lm(as.numeric(sum[,2])~ as.numeric(snpnum[,2])))
a <- summary(lm(as.numeric(ge[,2])~ as.numeric(sum[,2])))
a <- summary(lm(as.numeric(ge[,2])~ as.numeric(snpnum[,2])))

x <- append(x,a$adj.r.squared)

ggplot(data=all, aes(x=EndIndex.exclusive., y=V3))+
  geom_point(size=3,aes(x=EndIndex.exclusive., y=V3))+
  stat_smooth(method="lm")+
  xlab("length")+
  ylab("variance explained")+
  theme_classic()