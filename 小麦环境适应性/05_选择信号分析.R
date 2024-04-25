setwd("/Users/guoyafei/Desktop/RAiSD/poppairSelect/")
#5条传播路线的选择基因数量upset图 #####
#temp selection regions
library(UpSetR)

route01 <- read.table("poppair01.3k.gene", header=F,stringsAsFactors = F)
route01 <- route01[,1]
route02 <- read.table("poppair02.3k.gene", header=F,stringsAsFactors = F)
route02 <- route02[,1]
route03 <- read.table("poppair03.3k.gene", header=F,stringsAsFactors = F)
route03 <- route03[,1]
route04 <- read.table("poppair04.3k.gene", header=F,stringsAsFactors = F)
route04 <- route04[,1]
route05 <- read.table("poppair05.3k.gene", header=F,stringsAsFactors = F)
route05 <- route05[,1]

all <- list(
  route01 <- route01,
  route02 <- route02,
  route03 <- route03,
  route04 <- route04,
  route05 <- route05
)

names(all) <- c("route01","route02","route03","route04","route05")
upset(fromList(all))

upset(fromList(all), keep.order = TRUE,text.scale = c(1),point.size = 1.5, line.size = 1)


#5条传播路线重复选择的基因 #####
a <- c(1:5)
b <- c(1371,420,126,23,4)
data <- as.data.frame(cbind(a,b))
ggplot() + 
  geom_bar(data = data, aes(x = a, y = b),stat = "identity")+
  theme_classic()+
  xlab("Number of groups in which a gene under selection")+
  ylab("Number of genes")

random <- read.table("random.3k.select.gene", header=F,stringsAsFactors = F)
both <- random[which(random$V2 ==2),]
t.test(both$V1,mu=420)

#重复选择基因的RAiSD信号分布 #####
library(ggplot2)
#RAiSA delta u
data <- read.table("6pop.21chr.B.txt",header=F,stringsAsFactors = F, sep="\t")
#16	350186873	350193206	TraesCS3B02G566000
#16	350233874	350240793	TraesCS3B02G566100
#16	350238925	350242868	TraesCS3B02G566200
#16	350245397	350249402	TraesCS3B02G566300

#threshold
#24
#18
#24
#37
#54
sub <- data[which(data$V1 == 8),]

ggplot(sub, aes(x=V2, y=abs(V4))) +
  geom_point(size=0.8,color="#807dba")+
  theme_classic()+
  geom_hline(yintercept=24, color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=798342142,color='#E69F00')+
  geom_vline(xintercept=798404671,color='#E69F00')+
  ylim(0,300)+
  xlab("Positon(Mb)")+
  ylab("delta u")

ggplot(sub, aes(x=V2, y=abs(V5))) +
  geom_point(size=0.8,color="#807dba")+
  theme_classic()+
  geom_hline(yintercept=18, color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=798342142,color='#E69F00')+
  geom_vline(xintercept=798404671,color='#E69F00')+
  ylim(0,100)+
  xlab("Positon(Mb)")+
  ylab("delta u")

ggplot(sub, aes(x=V2, y=abs(V6))) +
  geom_point(size=0.8,color="#807dba")+
  theme_classic()+
  geom_hline(yintercept=24, color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=798342142,color='#E69F00')+
  geom_vline(xintercept=798404671,color='#E69F00')+
  ylim(0,200)+
  xlab("Positon(Mb)")+
  ylab("delta u")

ggplot(sub, aes(x=V2, y=abs(V7))) +
  geom_point(size=0.8,color="#807dba")+
  theme_classic()+
  geom_hline(yintercept=37, color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=798342142,color='#E69F00')+
  geom_vline(xintercept=798404671,color='#E69F00')+
  ylim(0,400)+
  xlab("Positon(Mb)")+
  ylab("delta u")

ggplot(sub, aes(x=V2, y=abs(V8))) +
  geom_point(size=0.8,color="#807dba")+
  theme_classic()+
  geom_hline(yintercept=54, color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=798342142,color='#E69F00')+
  geom_vline(xintercept=798404671,color='#E69F00')+
  ylim(0,400)+
  xlab("Positon(Mb)")+
  ylab("delta u")


#基因lassip信号分析#####
library(ggplot2)
setwd("/Users/guoyafei/Desktop/RAiSD/poppairSelect/lassip")
pop01 <- read.table("step5_21chr.gene.pop06.txt", header=T, stringsAsFactors = F)
sub <- pop01[which(pop01$pos > 797842142 & pop01$pos < 798904671),]
ggplot(sub, aes(x=pos, y=pop01_T)) +
    geom_line(size=0.8,color="#807dba")+
    theme_classic()+
    geom_vline(xintercept=798342142,color='#E69F00')+
    geom_vline(xintercept=798348475,color='#E69F00')+
    geom_vline(xintercept=798389143,color='#fc9272')+
    geom_vline(xintercept=798396062,color='#fc9272')+
    geom_vline(xintercept=798394194,color='#6baed6')+
    geom_vline(xintercept=798398137,color='#6baed6')+
    geom_vline(xintercept=798400666,color='#969696')+
    geom_vline(xintercept=798404671,color='#969696')+
  #ylim(0,300)+
  xlab("Positon(Mb)")+
    ylab("T")
ggplot(sub, aes(x=pos, y=pop01_m)) +
  geom_line(size=0.8,color="#807dba")+
  theme_classic()+
  geom_hline(yintercept=1, color='#E69F00',linetype = "dashed")+
  geom_vline(xintercept=798342142,color='#E69F00')+
  geom_vline(xintercept=798348475,color='#E69F00')+
  geom_vline(xintercept=798389143,color='#fc9272')+
  geom_vline(xintercept=798396062,color='#fc9272')+
  geom_vline(xintercept=798394194,color='#6baed6')+
  geom_vline(xintercept=798398137,color='#6baed6')+
  geom_vline(xintercept=798400666,color='#969696')+
  geom_vline(xintercept=798404671,color='#969696')+
  #ylim(0,300)+
  xlab("Positon(Mb)")+
  ylab("m")


#计算不同类型的趋同选择信号水平差异####
setwd("/Users/guoyafei/Desktop/RAiSD/poppairSelect/shuf")
two <- read.table("random.type1.3k.2sstdev.select.gene", header=F,stringsAsFactors = F)

#type1
num <- c(313,229,175,128,108)
fold <- c(3.279064,4.610619,6.230577,7.588421,9.712778)
x <- factor(c("2.5sstdev","3sstdev","3.5sstdev","4sstdev","4.5sstdev"),levels=(c("2.5sstdev","3sstdev","3.5sstdev","4sstdev","4.5sstdev")))
all1 <- as.data.frame(cbind(as.numeric(num),as.numeric(fold)))
all1$name <- x
all1$type <- "type1"

#type2
num <- c(312,237,179,141,116)
fold <- c(2.86381,4.18109,5.55934,7.3297,9.07484)
x <- factor(c("2.5sstdev","3sstdev","3.5sstdev","4sstdev","4.5sstdev"),levels=(c("2.5sstdev","3sstdev","3.5sstdev","4sstdev","4.5sstdev")))
all2 <- as.data.frame(cbind(as.numeric(num),as.numeric(fold)))
all2$name <- x
all2$type <- "type2"

#type3
num <- c(331,234,173,131,103)
fold <- c(3.15967,4.27888,5.62254,7.10545,8.33059)
x <- factor(c("2.5sstdev","3sstdev","3.5sstdev","4sstdev","4.5sstdev"),levels=(c("2.5sstdev","3sstdev","3.5sstdev","4sstdev","4.5sstdev")))
all3 <- as.data.frame(cbind(as.numeric(num),as.numeric(fold)))
all3$name <- x
all3$type <- "type3"


#type4
num <- c(270,187,147,107,78)
fold <- c(2.84015,3.76681,5.23835,6.33057,7.01404)
x <- factor(c("2.5sstdev","3sstdev","3.5sstdev","4sstdev","4.5sstdev"),levels=(c("2.5sstdev","3sstdev","3.5sstdev","4sstdev","4.5sstdev")))
all4 <- as.data.frame(cbind(as.numeric(num),as.numeric(fold)))
all4$name <- x
all4$type <- "type4"

all <- rbind(all1,all2,all3,all4)
ggplot(data=all,aes(x=name, y=V2, size = V1,color=type)) +
  #facet_grid(type~.)+
  geom_point(alpha=0.7) +
  scale_size(range = c(2, 10), name="Population (M)")+
  theme_classic()+
  ylim(2,10)+
  xlab("Threshold of delta u score")+
  ylab("Fold change")+
  scale_fill_manual(values = c("#ffff33","#ff7f00","#377eb8","#a65628")) +
  scale_color_manual(values = c("#ffff33","#ff7f00","#377eb8","#a65628") )

#受选择基因的表达差异分析####
setwd("/Users/guoyafei/Desktop/RAiSD/expression")
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
data <- read.table("type1.3k.3sstdev.geneExpressionSummary", header = T, stringsAsFactors = F)
data$selectTime <- as.factor(data$selectTime)
ggplot(data=data, aes(x=selectTime, y=TPM_mean,group=selectTime, fill=selectTime)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")


data <- read.table("type1.3k.3sstdev.geneExpressionTissue", sep="\t", check.names = F,header = T, stringsAsFactors = F)
data$selectTime <- as.factor(data$selectTime)
ggplot(data=data, aes(x=selectTime, y=data[,38],group=selectTime, fill=data[,38])) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  #theme_ipsum() +
  theme(legend.position="none", plot.title = element_text(size=11)) +
  ggtitle("A boxplot with jitter") +
  xlab("")

sub <- data[,c(1,2,5:38)]
cats <- melt(sub,id=c("Transcript", "selectTime"))

tissue <- cats[which(cats$variable != "Transcripts" & cats$variable != "Internode2" ),]
#tissue$variable <- as.character(tissue$variable)


sub <- tissue[which(tissue$selectTime > 1),]
ggplot(data=sub, aes(x=variable, y=value,group=Transcript, color=Transcript, fill=Transcript)) +
  geom_line() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  #theme_ipsum() +
  theme(legend.position="none", plot.title = element_text(size=11)) +
  ggtitle("A boxplot with jitter") +
  xlab("")

#一组基因在所有组织中的表达网络连接度分析
setwd("/Users/guoyafei/Desktop/network/tissue")
#读文件type1####
type1_data1 <- read.table("S1coleoptile.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data1$type <- "data1"
type1_data1$select <- NA
type1_data1[which(type1_data1$V2 == 1),6] <- "unique"
type1_data1[which(type1_data1$V2  > 1),6] <- "repeat"
type1_data1$four <- "type1"
type1_data2 <- read.table("S1root.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data2$type <- "data2"
type1_data2$select <- NA
type1_data2[which(type1_data2$V2 == 1),6] <- "unique"
type1_data2[which(type1_data2$V2  > 1),6] <- "repeat"
type1_data2$four <- "type1"
type1_data3 <- read.table("S2leaf.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data3$type <- "data3"
type1_data3$select <- NA
type1_data3[which(type1_data3$V2 == 1),6] <- "unique"
type1_data3[which(type1_data3$V2  > 1),6] <- "repeat"
type1_data3$four <- "type1"
type1_data4 <- read.table("S2root.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data4$type <- "data4"
type1_data4$select <- NA
type1_data4[which(type1_data4$V2 == 1),6] <- "unique"
type1_data4[which(type1_data4$V2  > 1),6] <- "repeat"
type1_data4$four <- "type1"
type1_data5 <- read.table("S3leaf.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data5$type <- "data5"
type1_data5$select <- NA
type1_data5[which(type1_data5$V2 == 1),6] <- "unique"
type1_data5[which(type1_data5$V2  > 1),6] <- "repeat"
type1_data5$four <- "type1"
type1_data6 <- read.table("S4Awn.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data6$type <- "data6"
type1_data6$select <- NA
type1_data6[which(type1_data6$V2 == 1),6] <- "unique"
type1_data6[which(type1_data6$V2  > 1),6] <- "repeat"
type1_data6$four <- "type1"
type1_data7 <- read.table("S4leaf.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data7$type <- "data7"
type1_data7$select <- NA
type1_data7[which(type1_data7$V2 == 1),6] <- "unique"
type1_data7[which(type1_data7$V2  > 1),6] <- "repeat"
type1_data7$four <- "type1"
type1_data8 <- read.table("S4spike.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data8$type <- "data8"
type1_data8$select <- NA
type1_data8[which(type1_data8$V2 == 1),6] <- "unique"
type1_data8[which(type1_data8$V2  > 1),6] <- "repeat"
type1_data8$four <- "type1"
type1_data9 <- read.table("S4stem.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data9$type <- "data9"
type1_data9$select <- NA
type1_data9[which(type1_data9$V2 == 1),6] <- "unique"
type1_data9[which(type1_data9$V2  > 1),6] <- "repeat"
type1_data9$four <- "type1"
type1_data10 <- read.table("S5Anthers.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data10$type <- "data10"
type1_data10$select <- NA
type1_data10[which(type1_data10$V2 == 1),6] <- "unique"
type1_data10[which(type1_data10$V2  > 1),6] <- "repeat"
type1_data10$four <- "type1"
type1_data11 <- read.table("S5leaf.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data11$type <- "data11"
type1_data11$select <- NA
type1_data11[which(type1_data11$V2 == 1),6] <- "unique"
type1_data11[which(type1_data11$V2  > 1),6] <- "repeat"
type1_data11$four <- "type1"
type1_data12 <- read.table("S6grain.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data12$type <- "data12"
type1_data12$select <- NA
type1_data12[which(type1_data12$V2 == 1),6] <- "unique"
type1_data12[which(type1_data12$V2  > 1),6] <- "repeat"
type1_data12$four <- "type1"
type1_data13 <- read.table("S6leaf.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data13$type <- "data13"
type1_data13$select <- NA
type1_data13[which(type1_data13$V2 == 1),6] <- "unique"
type1_data13[which(type1_data13$V2  > 1),6] <- "repeat"
type1_data13$four <- "type1"
type1_data14 <- read.table("S7grain.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data14$type <- "data14"
type1_data14$select <- NA
type1_data14[which(type1_data14$V2 == 1),6] <- "unique"
type1_data14[which(type1_data14$V2  > 1),6] <- "repeat"
type1_data14$four <- "type1"
type1_data15 <- read.table("S7leaf.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data15$type <- "data15"
type1_data15$select <- NA
type1_data15[which(type1_data15$V2 == 1),6] <- "unique"
type1_data15[which(type1_data15$V2  > 1),6] <- "repeat"
type1_data15$four <- "type1"
type1_data16 <- read.table("S8grain.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data16$type <- "data16"
type1_data16$select <- NA
type1_data16[which(type1_data16$V2 == 1),6] <- "unique"
type1_data16[which(type1_data16$V2  > 1),6] <- "repeat"
type1_data16$four <- "type1"
type1_data17 <- read.table("S8leaf.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data17$type <- "data17"
type1_data17$select <- NA
type1_data17[which(type1_data17$V2 == 1),6] <- "unique"
type1_data17[which(type1_data17$V2  > 1),6] <- "repeat"
type1_data17$four <- "type1"
type1_data18 <- read.table("S9grain.type1.both.allLine",header=F,stringsAsFactors = F)
type1_data18$type <- "data18"
type1_data18$select <- NA
type1_data18[which(type1_data18$V2 == 1),6] <- "unique"
type1_data18[which(type1_data18$V2  > 1),6] <- "repeat"
type1_data18$four <- "type1"
#读文件type2####
type2_data1 <- read.table("S1coleoptile.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data1$type <- "data1"
type2_data1$select <- NA
type2_data1[which(type2_data1$V2 == 1),6] <- "unique"
type2_data1[which(type2_data1$V2  > 1),6] <- "repeat"
type2_data1$four <- "type2"
type2_data2 <- read.table("S1root.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data2$type <- "data2"
type2_data2$select <- NA
type2_data2[which(type2_data2$V2 == 1),6] <- "unique"
type2_data2[which(type2_data2$V2  > 1),6] <- "repeat"
type2_data2$four <- "type2"
type2_data3 <- read.table("S2leaf.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data3$type <- "data3"
type2_data3$select <- NA
type2_data3[which(type2_data3$V2 == 1),6] <- "unique"
type2_data3[which(type2_data3$V2  > 1),6] <- "repeat"
type2_data3$four <- "type2"
type2_data4 <- read.table("S2root.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data4$type <- "data4"
type2_data4$select <- NA
type2_data4[which(type2_data4$V2 == 1),6] <- "unique"
type2_data4[which(type2_data4$V2  > 1),6] <- "repeat"
type2_data4$four <- "type2"
type2_data5 <- read.table("S3leaf.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data5$type <- "data5"
type2_data5$select <- NA
type2_data5[which(type2_data5$V2 == 1),6] <- "unique"
type2_data5[which(type2_data5$V2  > 1),6] <- "repeat"
type2_data5$four <- "type2"
type2_data6 <- read.table("S4Awn.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data6$type <- "data6"
type2_data6$select <- NA
type2_data6[which(type2_data6$V2 == 1),6] <- "unique"
type2_data6[which(type2_data6$V2  > 1),6] <- "repeat"
type2_data6$four <- "type2"
type2_data7 <- read.table("S4leaf.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data7$type <- "data7"
type2_data7$select <- NA
type2_data7[which(type2_data7$V2 == 1),6] <- "unique"
type2_data7[which(type2_data7$V2  > 1),6] <- "repeat"
type2_data7$four <- "type2"
type2_data8 <- read.table("S4spike.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data8$type <- "data8"
type2_data8$select <- NA
type2_data8[which(type2_data8$V2 == 1),6] <- "unique"
type2_data8[which(type2_data8$V2  > 1),6] <- "repeat"
type2_data8$four <- "type2"
type2_data9 <- read.table("S4stem.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data9$type <- "data9"
type2_data9$select <- NA
type2_data9[which(type2_data9$V2 == 1),6] <- "unique"
type2_data9[which(type2_data9$V2  > 1),6] <- "repeat"
type2_data9$four <- "type2"
type2_data10 <- read.table("S5Anthers.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data10$type <- "data10"
type2_data10$select <- NA
type2_data10[which(type2_data10$V2 == 1),6] <- "unique"
type2_data10[which(type2_data10$V2  > 1),6] <- "repeat"
type2_data10$four <- "type2"
type2_data11 <- read.table("S5leaf.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data11$type <- "data11"
type2_data11$select <- NA
type2_data11[which(type2_data11$V2 == 1),6] <- "unique"
type2_data11[which(type2_data11$V2  > 1),6] <- "repeat"
type2_data11$four <- "type2"
type2_data12 <- read.table("S6grain.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data12$type <- "data12"
type2_data12$select <- NA
type2_data12[which(type2_data12$V2 == 1),6] <- "unique"
type2_data12[which(type2_data12$V2  > 1),6] <- "repeat"
type2_data12$four <- "type2"
type2_data13 <- read.table("S6leaf.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data13$type <- "data13"
type2_data13$select <- NA
type2_data13[which(type2_data13$V2 == 1),6] <- "unique"
type2_data13[which(type2_data13$V2  > 1),6] <- "repeat"
type2_data13$four <- "type2"
type2_data14 <- read.table("S7grain.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data14$type <- "data14"
type2_data14$select <- NA
type2_data14[which(type2_data14$V2 == 1),6] <- "unique"
type2_data14[which(type2_data14$V2  > 1),6] <- "repeat"
type2_data14$four <- "type2"
type2_data15 <- read.table("S7leaf.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data15$type <- "data15"
type2_data15$select <- NA
type2_data15[which(type2_data15$V2 == 1),6] <- "unique"
type2_data15[which(type2_data15$V2  > 1),6] <- "repeat"
type2_data15$four <- "type2"
type2_data16 <- read.table("S8grain.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data16$type <- "data16"
type2_data16$select <- NA
type2_data16[which(type2_data16$V2 == 1),6] <- "unique"
type2_data16[which(type2_data16$V2  > 1),6] <- "repeat"
type2_data16$four <- "type2"
type2_data17 <- read.table("S8leaf.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data17$type <- "data17"
type2_data17$select <- NA
type2_data17[which(type2_data17$V2 == 1),6] <- "unique"
type2_data17[which(type2_data17$V2  > 1),6] <- "repeat"
type2_data17$four <- "type2"
type2_data18 <- read.table("S9grain.type2.both.allLine",header=F,stringsAsFactors = F)
type2_data18$type <- "data18"
type2_data18$select <- NA
type2_data18[which(type2_data18$V2 == 1),6] <- "unique"
type2_data18[which(type2_data18$V2  > 1),6] <- "repeat"
type2_data18$four <- "type2"
#读文件type3####
type3_data1 <- read.table("S1coleoptile.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data1$type <- "data1"
type3_data1$select <- NA
type3_data1[which(type3_data1$V2 == 1),6] <- "unique"
type3_data1[which(type3_data1$V2  > 1),6] <- "repeat"
type3_data1$four <- "type3"
type3_data2 <- read.table("S1root.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data2$type <- "data2"
type3_data2$select <- NA
type3_data2[which(type3_data2$V2 == 1),6] <- "unique"
type3_data2[which(type3_data2$V2  > 1),6] <- "repeat"
type3_data2$four <- "type3"
type3_data3 <- read.table("S2leaf.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data3$type <- "data3"
type3_data3$select <- NA
type3_data3[which(type3_data3$V2 == 1),6] <- "unique"
type3_data3[which(type3_data3$V2  > 1),6] <- "repeat"
type3_data3$four <- "type3"
type3_data4 <- read.table("S2root.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data4$type <- "data4"
type3_data4$select <- NA
type3_data4[which(type3_data4$V2 == 1),6] <- "unique"
type3_data4[which(type3_data4$V2  > 1),6] <- "repeat"
type3_data4$four <- "type3"
type3_data5 <- read.table("S3leaf.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data5$type <- "data5"
type3_data5$select <- NA
type3_data5[which(type3_data5$V2 == 1),6] <- "unique"
type3_data5[which(type3_data5$V2  > 1),6] <- "repeat"
type3_data5$four <- "type3"
type3_data6 <- read.table("S4Awn.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data6$type <- "data6"
type3_data6$select <- NA
type3_data6[which(type3_data6$V2 == 1),6] <- "unique"
type3_data6[which(type3_data6$V2  > 1),6] <- "repeat"
type3_data6$four <- "type3"
type3_data7 <- read.table("S4leaf.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data7$type <- "data7"
type3_data7$select <- NA
type3_data7[which(type3_data7$V2 == 1),6] <- "unique"
type3_data7[which(type3_data7$V2  > 1),6] <- "repeat"
type3_data7$four <- "type3"
type3_data8 <- read.table("S4spike.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data8$type <- "data8"
type3_data8$select <- NA
type3_data8[which(type3_data8$V2 == 1),6] <- "unique"
type3_data8[which(type3_data8$V2  > 1),6] <- "repeat"
type3_data8$four <- "type3"
type3_data9 <- read.table("S4stem.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data9$type <- "data9"
type3_data9$select <- NA
type3_data9[which(type3_data9$V2 == 1),6] <- "unique"
type3_data9[which(type3_data9$V2  > 1),6] <- "repeat"
type3_data9$four <- "type3"
type3_data10 <- read.table("S5Anthers.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data10$type <- "data10"
type3_data10$select <- NA
type3_data10[which(type3_data10$V2 == 1),6] <- "unique"
type3_data10[which(type3_data10$V2  > 1),6] <- "repeat"
type3_data10$four <- "type3"
type3_data11 <- read.table("S5leaf.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data11$type <- "data11"
type3_data11$select <- NA
type3_data11[which(type3_data11$V2 == 1),6] <- "unique"
type3_data11[which(type3_data11$V2  > 1),6] <- "repeat"
type3_data11$four <- "type3"
type3_data12 <- read.table("S6grain.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data12$type <- "data12"
type3_data12$select <- NA
type3_data12[which(type3_data12$V2 == 1),6] <- "unique"
type3_data12[which(type3_data12$V2  > 1),6] <- "repeat"
type3_data12$four <- "type3"
type3_data13 <- read.table("S6leaf.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data13$type <- "data13"
type3_data13$select <- NA
type3_data13[which(type3_data13$V2 == 1),6] <- "unique"
type3_data13[which(type3_data13$V2  > 1),6] <- "repeat"
type3_data13$four <- "type3"
type3_data14 <- read.table("S7grain.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data14$type <- "data14"
type3_data14$select <- NA
type3_data14[which(type3_data14$V2 == 1),6] <- "unique"
type3_data14[which(type3_data14$V2  > 1),6] <- "repeat"
type3_data14$four <- "type3"
type3_data15 <- read.table("S7leaf.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data15$type <- "data15"
type3_data15$select <- NA
type3_data15[which(type3_data15$V2 == 1),6] <- "unique"
type3_data15[which(type3_data15$V2  > 1),6] <- "repeat"
type3_data15$four <- "type3"
type3_data16 <- read.table("S8grain.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data16$type <- "data16"
type3_data16$select <- NA
type3_data16[which(type3_data16$V2 == 1),6] <- "unique"
type3_data16[which(type3_data16$V2  > 1),6] <- "repeat"
type3_data16$four <- "type3"
type3_data17 <- read.table("S8leaf.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data17$type <- "data17"
type3_data17$select <- NA
type3_data17[which(type3_data17$V2 == 1),6] <- "unique"
type3_data17[which(type3_data17$V2  > 1),6] <- "repeat"
type3_data17$four <- "type3"
type3_data18 <- read.table("S9grain.type3.both.allLine",header=F,stringsAsFactors = F)
type3_data18$type <- "data18"
type3_data18$select <- NA
type3_data18[which(type3_data18$V2 == 1),6] <- "unique"
type3_data18[which(type3_data18$V2  > 1),6] <- "repeat"
type3_data18$four <- "type3"
#读文件type4####
type4_data1 <- read.table("S1coleoptile.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data1$type <- "data1"
type4_data1$select <- NA
type4_data1[which(type4_data1$V2 == 1),6] <- "unique"
type4_data1[which(type4_data1$V2  > 1),6] <- "repeat"
type4_data1$four <- "type4"
type4_data2 <- read.table("S1root.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data2$type <- "data2"
type4_data2$select <- NA
type4_data2[which(type4_data2$V2 == 1),6] <- "unique"
type4_data2[which(type4_data2$V2  > 1),6] <- "repeat"
type4_data2$four <- "type4"
type4_data3 <- read.table("S2leaf.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data3$type <- "data3"
type4_data3$select <- NA
type4_data3[which(type4_data3$V2 == 1),6] <- "unique"
type4_data3[which(type4_data3$V2  > 1),6] <- "repeat"
type4_data3$four <- "type4"
type4_data4 <- read.table("S2root.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data4$type <- "data4"
type4_data4$select <- NA
type4_data4[which(type4_data4$V2 == 1),6] <- "unique"
type4_data4[which(type4_data4$V2  > 1),6] <- "repeat"
type4_data4$four <- "type4"
type4_data5 <- read.table("S3leaf.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data5$type <- "data5"
type4_data5$select <- NA
type4_data5[which(type4_data5$V2 == 1),6] <- "unique"
type4_data5[which(type4_data5$V2  > 1),6] <- "repeat"
type4_data5$four <- "type4"
type4_data6 <- read.table("S4Awn.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data6$type <- "data6"
type4_data6$select <- NA
type4_data6[which(type4_data6$V2 == 1),6] <- "unique"
type4_data6[which(type4_data6$V2  > 1),6] <- "repeat"
type4_data6$four <- "type4"
type4_data7 <- read.table("S4leaf.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data7$type <- "data7"
type4_data7$select <- NA
type4_data7[which(type4_data7$V2 == 1),6] <- "unique"
type4_data7[which(type4_data7$V2  > 1),6] <- "repeat"
type4_data7$four <- "type4"
type4_data8 <- read.table("S4spike.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data8$type <- "data8"
type4_data8$select <- NA
type4_data8[which(type4_data8$V2 == 1),6] <- "unique"
type4_data8[which(type4_data8$V2  > 1),6] <- "repeat"
type4_data8$four <- "type4"
type4_data9 <- read.table("S4stem.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data9$type <- "data9"
type4_data9$select <- NA
type4_data9[which(type4_data9$V2 == 1),6] <- "unique"
type4_data9[which(type4_data9$V2  > 1),6] <- "repeat"
type4_data9$four <- "type4"
type4_data10 <- read.table("S5Anthers.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data10$type <- "data10"
type4_data10$select <- NA
type4_data10[which(type4_data10$V2 == 1),6] <- "unique"
type4_data10[which(type4_data10$V2  > 1),6] <- "repeat"
type4_data10$four <- "type4"
type4_data11 <- read.table("S5leaf.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data11$type <- "data11"
type4_data11$select <- NA
type4_data11[which(type4_data11$V2 == 1),6] <- "unique"
type4_data11[which(type4_data11$V2  > 1),6] <- "repeat"
type4_data11$four <- "type4"
type4_data12 <- read.table("S6grain.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data12$type <- "data12"
type4_data12$select <- NA
type4_data12[which(type4_data12$V2 == 1),6] <- "unique"
type4_data12[which(type4_data12$V2  > 1),6] <- "repeat"
type4_data12$four <- "type4"
type4_data13 <- read.table("S6leaf.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data13$type <- "data13"
type4_data13$select <- NA
type4_data13[which(type4_data13$V2 == 1),6] <- "unique"
type4_data13[which(type4_data13$V2  > 1),6] <- "repeat"
type4_data13$four <- "type4"
type4_data14 <- read.table("S7grain.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data14$type <- "data14"
type4_data14$select <- NA
type4_data14[which(type4_data14$V2 == 1),6] <- "unique"
type4_data14[which(type4_data14$V2  > 1),6] <- "repeat"
type4_data14$four <- "type4"
type4_data15 <- read.table("S7leaf.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data15$type <- "data15"
type4_data15$select <- NA
type4_data15[which(type4_data15$V2 == 1),6] <- "unique"
type4_data15[which(type4_data15$V2  > 1),6] <- "repeat"
type4_data15$four <- "type4"
type4_data16 <- read.table("S8grain.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data16$type <- "data16"
type4_data16$select <- NA
type4_data16[which(type4_data16$V2 == 1),6] <- "unique"
type4_data16[which(type4_data16$V2  > 1),6] <- "repeat"
type4_data16$four <- "type4"
type4_data17 <- read.table("S8leaf.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data17$type <- "data17"
type4_data17$select <- NA
type4_data17[which(type4_data17$V2 == 1),6] <- "unique"
type4_data17[which(type4_data17$V2  > 1),6] <- "repeat"
type4_data17$four <- "type4"
type4_data18 <- read.table("S9grain.type4.both.allLine",header=F,stringsAsFactors = F)
type4_data18$type <- "data18"
type4_data18$select <- NA
type4_data18[which(type4_data18$V2 == 1),6] <- "unique"
type4_data18[which(type4_data18$V2  > 1),6] <- "repeat"
type4_data18$four <- "type4"
#画图#####
all <- rbind(type1_data1,type1_data2,type1_data3,type1_data4,type1_data5,type1_data6,type1_data7,type1_data8,type1_data9,type1_data10,type1_data11,type1_data12,type1_data13,type1_data14,type1_data15,type1_data16,type1_data17,type1_data18,type2_data1,type2_data2,type2_data3,type2_data4,type2_data5,type2_data6,type2_data7,type2_data8,type2_data9,type2_data10,type2_data11,type2_data12,type2_data13,type2_data14,type2_data15,type2_data16,type2_data17,type2_data18,type3_data1,type3_data2,type3_data3,type3_data4,type3_data5,type3_data6,type3_data7,type3_data8,type3_data9,type3_data10,type3_data11,type3_data12,type3_data13,type3_data14,type3_data15,type3_data16,type3_data17,type3_data18,type4_data1,type4_data2,type4_data3,type4_data4,type4_data5,type4_data6,type4_data7,type4_data8,type4_data9,type4_data10,type4_data11,type4_data12,type4_data13,type4_data14,type4_data15,type4_data16,type4_data17,type4_data18)
all$type <- factor(all$type, levels=c("data1","data2","data3","data4","data5","data6","data7","data8","data9","data10","data11","data12","data13","data14","data15","data16","data17","data18"))
#胚，根，叶，种子
all <- rbind(type1_data1,type1_data4,type1_data5,type1_data12,type2_data1,type2_data4,type2_data5,type2_data12,type3_data1,type3_data4,type3_data5,type3_data12,type4_data1,type4_data4,type4_data5,type4_data12)
all$type <- factor(all$type, levels=c("data1","data4","data5","data12"))
#3,7,13,17不同时期的叶子
all <- rbind(type1_data3,type1_data7,type1_data13,type1_data17,type2_data3,type2_data7,type2_data13,type2_data17,type3_data3,type3_data7,type3_data13,type3_data17,type4_data3,type4_data7,type4_data13,type4_data17)
all$type <- factor(all$type, levels=c("data3","data7","data13","data17"))
ggplot(data=all, aes(x = V3,group=type,color=type)) +
  facet_grid(four~select)+
  geom_density(binwidth = 0.01)+
  ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("")
ggplot(data=all, aes(x = V3, group = select, color = select)) +
  facet_grid(type~.,scales="free")+
  geom_density(binwidth = 0.01)+
  ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("")
t.test(data[which(data$type == "unique"),4],data[which(data$type == "repeat"),4])

#受选择基因的GO富集分析####
setwd("/Users/guoyafei/Desktop/network/GO")
data <- read.table("type1.repeat.GO",header=T,sep="\t",stringsAsFactors = F)
ggplot(data=data)+
  geom_bar(aes(y= Description, x=Count, fill=pvalue), stat='identity') +
  #coord_flip() +
  scale_fill_gradient(expression(pvalue),low="#FDAE6B", high = "#D94801") +
  ylab("Gene number") +
  xlab("GO term description") +
  theme(
    axis.text.x=element_text(angle=45,hjust=1,color="black"),
    axis.text.y=element_text(angle=45,color="black", size=rel(0.3)),
    axis.title.x = element_text(angle=45,color="black", size=rel(5)),
    axis.title.y = element_blank(),
    legend.text=element_text(color="black",size=rel(0.2)),
    legend.title = element_text(color="black",size=rel(0.7))
  )+
  #ylim(0,10)+
  theme_bw()

#受选择基因的表达网络连接度分析####
setwd("/Users/guoyafei/Desktop/network/type")

data <- read.table("S2leaf.type2.addTime.txt", header=T,stringsAsFactors = F)
data$SelectTime <- as.numeric(data$SelectTime)

data$type <- NA
data[which(data$SelectTime == 1), 8] <- "unique"
data[which(data$SelectTime > 1), 8] <- "repeat"
ggplot(data=data, aes(x = Eigenvector_Centrality,group=type,color=type)) +
  #facet_grid(four~select)+
  geom_density(binwidth = 0.01)+
  ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("")






#受选择基因的表达网络连接度分析2####
setwd("/Users/guoyafei/Desktop/network/S2leaf/")
data <- read.table("type2.S2leaf.uniqGene.network", header=T,stringsAsFactors = F)
data$SelectTime <- as.numeric(data$SelectTime)

data$type <- NA
data[which(data$SelectTime == 1), 18] <- "unique"
data[which(data$SelectTime > 1), 18] <- "repeat"
ggplot(data=data, aes(x = min_GWAS_P,color=max_eBPis)) +
  #facet_grid(type~.)+
  geom_density(binwidth = 0.01)+
  ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("")

#ggplot(data = sub, aes(x = min_GWAS_P, y=absmax_GWAS_BETA, color=max_BF)) +
ggplot(data = sub, aes(x = Eigenvector_Centrality, y=median_GWAS_BETA, color=max_BF)) +
    #facet_grid(type~.)+
  geom_point()+
  ggtitle("A boxplot with jitter") +
  theme_classic()

colnames(data)
[1] "gene"                   "max_BF"                 "absmax_BETA"            "min_SD_Beta_is"         "max_eBPis"             
[6] "median_GWAS_BETA"       "absmax_GWAS_BETA"       "median_GWAS_P"          "min_GWAS_P"             "Node"                  
[11] "Betweenness_Centrality" "Closeness_Centrality"   "Eigenvector_Centrality" "SelectTime"             "MeanR"                 
[16] "MeanLine"               "gene.1"                 "type" 

sub <- data[which(data$min_GWAS_P < 0.05),]

#受选择基因的表达网络连接度分析3####
#环境相关的基因在共表达网络上的富集情况
library(ggplot2)
setwd("/Users/guoyafei/Desktop/network/module/demo")
pdf("module04.pdf",height=3.5, width=4.5)
data <- read.table("module04_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=564)
ggplot(data,aes(V1))+
  geom_density(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module04")+
  geom_vline(xintercept = 564, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module06.pdf",height=3.5, width=4.5)
data <- read.table("module06_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=503)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module06")+
  geom_vline(xintercept = 503, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module07.pdf",height=3.5, width=4.5)
data <- read.table("module07_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=461)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module07")+
  geom_vline(xintercept = 461, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module08.pdf",height=3.5, width=4.5)
data <- read.table("module08_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=439)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module08")+
  geom_vline(xintercept = 439, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module10.pdf",height=3.5, width=4.5)
data <- read.table("module10_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=388)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module10")+
  geom_vline(xintercept = 388, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module13.pdf",height=3.5, width=4.5)
data <- read.table("module13_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=314)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module13")+
  geom_vline(xintercept = 314, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module14.pdf",height=3.5, width=4.5)
data <- read.table("module14_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=238)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module14")+
  geom_vline(xintercept = 238, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module21.pdf",height=3.5, width=4.5)
data <- read.table("module21_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=138)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module21")+
  geom_vline(xintercept = 138, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module25.pdf",height=3.5, width=4.5)
data <- read.table("module25_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=113)
ggplot(data,aes(V1))+
  geom_density(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module25")+
  geom_vline(xintercept = 113, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module26.pdf",height=3.5, width=4.5)
data <- read.table("module26_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=108)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module26")+
  geom_vline(xintercept = 108, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module27.pdf",height=3.5, width=4.5)
data <- read.table("module27_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=100)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module27")+
  geom_vline(xintercept = 100, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module28.pdf",height=3.5, width=4.5)
data <- read.table("module28_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=93)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module28")+
  geom_vline(xintercept = 93, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module32.pdf",height=3.5, width=4.5)
data <- read.table("module32_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=80)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module32")+
  geom_vline(xintercept = 80, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()


pdf("module45.pdf",height=3.5, width=4.5)
data <- read.table("module45_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=47)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module45")+
  geom_vline(xintercept = 47, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module46.pdf",height=3.5, width=4.5)
data <- read.table("module46_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=42)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module46")+
  geom_vline(xintercept = 42, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module50.pdf",height=3.5, width=4.5)
data <- read.table("module50_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=41)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module50")+
  geom_vline(xintercept = 41, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module51.pdf",height=3.5, width=4.5)
data <- read.table("module51_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=39)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module51")+
  geom_vline(xintercept = 39, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module54.pdf",height=3.5, width=4.5)
data <- read.table("module54_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=30)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module54")+
  geom_vline(xintercept = 30, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module57.pdf",height=3.5, width=4.5)
data <- read.table("module57_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=24)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module57")+
  geom_vline(xintercept = 24, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()

pdf("module58.pdf",height=3.5, width=4.5)
data <- read.table("module58_shuf.txt", header=F, stringsAsFactors = F)
t.test(data$V1,mu=26)
ggplot(data,aes(V1))+
  geom_histogram(bins=20,aes(x=V1, y=..density..),position="identity",alpha = 0.6)+
  theme_bw()+
  xlab("Number of Genes in Module58")+
  geom_vline(xintercept = 26, color = '#a50026', size = 1.5,linetype = "dashed") +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
dev.off()
#新的一部分########
data <- read.table("/Users/guoyafei/Desktop/network/module/module07.r03.network_result.txt3", header=F,stringsAsFactors = F)
data$V5 <- as.factor(data$V5)
ggplot(data=data, aes(y = V2, color=V5,group=V5)) +
  #facet_grid(type~.)+
  geom_boxplot(binwidth = 0.01)+
  #ggtitle("A boxplot with jitter") +
  theme_classic()+
  ylab("Betweenness Centrality")+
  scale_color_manual(values = c("#034E7B","#99000D")) +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

ggplot(data=data, aes(x=V7,y = V2)) +
  #facet_grid(type~.)+
  geom_point(binwidth = 0.01)+
  #ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("Edge number")+
  ylab("Betweenness Centrality")+
  scale_color_manual(values = c("#034E7B","#99000D")) +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

ggplot(data=data, aes(x=V7,y = V3)) +
  #facet_grid(type~.)+
  geom_point(binwidth = 0.01)+
  #ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("Edge number")+
  ylab("Closeness Centrality")+
  scale_color_manual(values = c("#034E7B","#99000D")) +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

ggplot(data=data, aes(x=V7,y = V4)) +
  #facet_grid(type~.)+
  geom_point(binwidth = 0.01)+
  #ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("Edge number")+
  ylab("Eigenvector Centrality")+
  scale_color_manual(values = c("#034E7B","#99000D")) +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

ggplot(data=data, aes(x=V2, y = V4)) +
  #facet_grid(type~.)+
  geom_point(binwidth = 0.01)+
  #ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("Betweenness Centrality")+
  ylab("Eigenvector Centrality")+
  scale_color_manual(values = c("#034E7B","#99000D")) +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

ggplot(data=data, aes(x=V3, y = V4)) +
  #facet_grid(type~.)+
  geom_point(binwidth = 0.01)+
  #ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("Closeness Centrality")+
  ylab("Eigenvector Centrality")+
  scale_color_manual(values = c("#034E7B","#99000D")) +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

ggplot(data=data, aes(x=V2, y = V3)) +
  #facet_grid(type~.)+
  geom_point(binwidth = 0.01)+
  #ggtitle("A boxplot with jitter") +
  theme_classic()+
  xlab("Betweenness Centrality")+
  ylab("Closeness Centrality")+
  scale_color_manual(values = c("#034E7B","#99000D")) +
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))


t.test(data[which(data$V5== "0"),4], data[which(data$V5== "1"),4])






