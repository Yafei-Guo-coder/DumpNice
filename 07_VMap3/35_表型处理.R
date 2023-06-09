#表型数据处理
setwd("/Users/guoyafei/Desktop/表型")
library(ggplot2)
data <- read.table("/Users/guoyafei/Desktop/表型/heading.day.txt", header=T,stringsAsFactors = F)
ggplot(data, aes(x = Q)) +
  geom_histogram(binwidth = 1, fill = "lightblue", colour = "black")+
  theme_classic()
shapiro.test(data$S)

qqnorm(data$P, ylab = 'TaxiIn') # 做滑行时间的QQ图，将y轴命名为TaxiIn
qqline(data$P)

qqnorm(test, ylab = 'TaxiIn') # 做滑行时间的QQ图，将y轴命名为TaxiIn
qqline(test)

test <- runif(100, min = 2, max = 4)
test <- rnorm(100, mean = 5, sd = 3)
shapiro.test(rnorm(100, mean = 5, sd = 3))
shapiro.test(runif(100, min = 2, max = 4))


#BLUE
library(lme4)
library(emmeans)
library(data.table)
library(tidyverse)

dat = fread("flowering.day.txt",data.table = F)
#col = 1:5
#dat[,col] = dat %>% select(all_of(col)) %>% map_df(as.factor)
#str(dat)
m1 = lmer(day ~ line + (1|env:line) + (1|env) + (1|env:rep), data=dat)
re1 = as.data.frame(emmeans(m1,"line"))
write.table(re1,"flowering.day.blue.txt", quote=F,sep="\t",row.names=F)

dat = fread("green.day.txt",data.table = F)
#col = 1:5
#dat[,col] = dat %>% select(all_of(col)) %>% map_df(as.factor)
#str(dat)
m1 = lmer(day ~ line + (1|env:line) + (1|env) + (1|env:rep), data=dat)
re1 = as.data.frame(emmeans(m1,"line"))
write.table(re1,"green.day.blue.txt", quote=F,sep="\t",row.names=F)

dat = fread("heading.day.txt",data.table = F)
#col = 1:5
#dat[,col] = dat %>% select(all_of(col)) %>% map_df(as.factor)
#str(dat)
m1 = lmer(day ~ line + (1|env:line) + (1|env) + (1|env:rep), data=dat)
re1 = as.data.frame(emmeans(m1,"line"))
write.table(re1,"heading.day.blue.txt", quote=F,sep="\t",row.names=F)


#blue 值和实际值的关系
setwd("/Users/guoyafei/Desktop/表型/trait")
data1 <- read.table("heading.blue.txt", header=F, stringsAsFactors = F)
data2 <- read.table("/Users/guoyafei/Desktop/表型/heading.day.txt", header=T,stringsAsFactors = F)

ggplot(data1,aes(x = V2))+
  geom_density()+
  theme_classic()
ggplot(data2,aes(x = S))+
  geom_density()+
  theme_classic()
#flowering是个双峰分布

setwd("/Users/guoyafei/Desktop/表型/group")

data <- read.table("/Users/guoyafei/Desktop/表型/group/green.day.blue.txt",header=F,stringsAsFactors = F)

land <- data[which(data$V3 == "Landrace"),]
cul <- data[which(data$V3 == "Cultivar"),]

ggplot(land,aes(x = V2))+
  geom_density()+
  theme_classic()
ggplot(cul,aes(x = V2))+
  geom_density()+
  theme_classic()
ggplot(data,aes(x = V2))+
  geom_density()+
  theme_classic()


#----------------------------------check trait-----------------------
#table
setwd("/Users/guoyafei/Desktop/表型/trait/raw")
head <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/heading.blue.txt", header=F, stringsAsFactors = F)
a <- (head$P + head$Q)/2
b <- (head$R + head$S)/2
head$mean1 <- (head$P + head$Q)/2
head$mean2 <- (head$R + head$S)/2
sub <- head[which(head$mean1 != "NA" & head$mean2 != "NA"),]
cor(sub$mean1, sub$mean2)
head_blue <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/yellow.blue.txt",header=F,stringsAsFactors = F)
a <- head_blue$V2
mean(na.omit(a))
mean(na.omit(b))
sd(na.omit(a))
sd(na.omit(b))
max(na.omit(a))
max(na.omit(b))
min(na.omit(a))
min(na.omit(b))

#plot
library(ggplot2)
library(dplyr)
library(hrbrthemes)
setwd("/Users/guoyafei/Desktop/表型/trait/plot")

cultivar <- read.table("/Users/guoyafei/Desktop/表型/group/cultivar.txt", header= F, stringsAsFactors = F)
landrace <- read.table("/Users/guoyafei/Desktop/表型/group/landrace.txt", header= F, stringsAsFactors = F)

head <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/heading-flowering.blue.txt", header=F, stringsAsFactors = F)
rownames(head) <- head$V1
head$type <- "All"
land_head <- head[landrace$V1,]
land_head$type  <- "Landrace"
cul_head <- head[cultivar$V1,]
cul_head$type <- "Cultivar"

all <- rbind(head,land_head,cul_head)
all <- rbind(land_head,cul_head)

p <- head %>%
  ggplot( ) +
  #geom_histogram(data=head, aes(x=V2,fill=type), alpha=0.8) +
  geom_density(data=head, color="#999999",aes(x=V2,fill=type), alpha=0.6)+
  scale_fill_manual(values=c("#999999")) +
  #theme_ipsum() +
  labs(fill="")+
  theme_bw()
pdf("heading-flowering1.pdf",width=4.63, height=3.13)
print(p)
dev.off()
p <- all %>%
  ggplot( aes(x=V2, fill=type)) +
  geom_histogram(color="#e9ecef",  alpha=0.8, position = 'identity') +
  #geom_density(data=head, aes(x=V2), alpha=0.8)+
  scale_fill_manual(values=c("#0072B2", "#E69F00")) +
  #theme_ipsum() +
  labs(fill="")+
  theme_bw()
pdf("heading-flowering2.pdf",width=4.63, height=3.13)
print(p)
dev.off()

#表型分布规律
setwd("/Users/guoyafei/Desktop/表型/trait/分布差异")
cultivar <- read.table("/Users/guoyafei/Desktop/表型/group/cultivar.txt", header= F, stringsAsFactors = F)
landrace <- read.table("/Users/guoyafei/Desktop/表型/group/landrace.txt", header= F, stringsAsFactors = F)

head <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/yellow.blue.txt", header=F, stringsAsFactors = F)
rownames(head) <- head$V1
head$type <- "All"
land_head <- head[landrace$V1,]
land_head$type  <- "Landrace"
cul_head <- head[cultivar$V1,]
cul_head$type <- "Cultivar"
#p.value
t.test(land_head$V2,cul_head$V2)
#heading: < 2.2e-16
#heading-flowering: 4.159e-15
#flowering: < 2.2e-16
#flowering-green: < 2.2e-16
#green: < 2.2e-16
#green-yellow: 0.4543
#yellow: < 2.2e-16

#mean.value (landrace VS cultivar)(相对于landrace, cultivar发生了什么变化)
mean(land_head$V2, na.rm = T)
mean(cul_head$V2, na.rm = T)
heading: 214.7068 208.1617（抽穗提前6天）
heading-flowering: 3.912943 4.663681（生育期阶段一延长1天：）
flowering: 218.5399 212.8254（开花提前）
flowering-green: 29.29549 31.64418（生育期阶段二延长2天：灌浆时间变长）
green: 247.8354 244.4695（灌浆完成提前3天）
green-yellow: 2.657951 2.538757（硬化时间一致）
yellow: 249.9636 246.823（成熟时间提前3天）


