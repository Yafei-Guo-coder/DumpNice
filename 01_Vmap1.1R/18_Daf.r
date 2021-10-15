setwd("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/01_BasicStatistic/07_Daf")
A <- read.table("A_daf",header=F,stringsAsFactors = F)
B <- read.table("B_daf",header=F,stringsAsFactors = F)
D <- read.table("D_daf",header=F,stringsAsFactors = F)
library("ggplot2")
ggplot(all, aes(x=V1,fill=lineage)) +
  facet_grid(lineage ~ .)+
  geom_histogram(aes(y=..density..),binwidth=0.02, alpha=.7)+ 
  theme_classic()+
  theme(axis.text.x=element_text(size=200),
        axis.text.y=element_text(size=300,face="plain"),
        axis.title.y=element_text(size = 300,face="plain"))
