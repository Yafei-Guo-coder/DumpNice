setwd("/Users/guoyafei/Documents/lab/EGWAS/forYafei/")
#data <- read.csv("R1_count_statics.txt",header=T,stringsAsFactors = F,sep="\t")
#data <- read.csv("bam_size.txt",header=T,stringsAsFactors = F,sep="\t")
data <- read.csv("pileup_statics.txt",header=T,stringsAsFactors = F,sep="\t")
library(ggplot2)
#ggplot(data,aes(Reads.count))+geom_histogram(binwidth=0.01)

data$Taxon.name.5700 <- as.factor(data$Taxon.name.5700)
data$Flowcell.ID <- as.factor(data$Flowcell.ID)
data$Lane <- as.factor(data$Lane)
data$Library.ID <- as.factor(data$Library.ID)
data$Company <- as.factor(data$Company)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#直方图

facet_grid(Company~.)+
  scale_color_manual(values=cbPalette)

p+labs(x = "log(Reads count)", y = "Individual count", title = "Reads count at SNP sites in 5700 Individuals")+
  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),plot.title = element_text(hjust = 0.5,vjust=0.5,size=20),axis.title.x = element_text(size = 17),axis.title.y = element_text(size = 17))+ 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  geom_vline(aes(xintercept=2.4548), colour="#990000", linetype="dashed")+
  geom_vline(aes(xintercept=2.7559), colour="#E69F00", linetype="dashed")

#箱线图
pbox <- ggplot(data,aes(Company,Sequencing_error,color=Company))+
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+
  scale_color_manual(values=cbPalette)+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = "Error Frrequency Distribution of 124 Libraries")+
  theme(axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),plot.title = element_text(hjust = 0.5,vjust=0.5,size=20),axis.title.x = element_text(size = 17),axis.title.y = element_text(size = 17),legend.text = element_text(size = 13))+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))

#多因素方差分析
Result<-aov(Reads.count~Library.ID+Company+Library.ID:Company,data=data)
anova(Result)
interaction.plot(data$Company,data$Library.ID,data$log.bam_Mb.,type="b",legend =F,xlab="Company",ylab="log bamfile size(Mb)",main="The bamfile size(Mb) Distribution of 11,976 Individuals",cex.lab=1.5, cex.axis=1.5, cex.main=1.7, cex.sub=1.5)
#单因素方差分析
OneWay<-aov(Reads.count~Company,data=my)
anova(OneWay)
summary(OneWay)

library("gplots")
plotmeans(log.Reads.count~Company,data=data,p=0.95,use.t=TRUE)


