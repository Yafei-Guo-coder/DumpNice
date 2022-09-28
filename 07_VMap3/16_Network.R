setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network")
data <- read.table("shuf.degree.txt", header=F,stringsAsFactors = F)
library(ggplot2)
x<-data$V1
funss<-function(x) 1/sqrt(2*pi)*exp(-1/2*x^2)
ggplot(data, aes(x=V1))+geom_histogram(bins=50,aes(y=..density..))+stat_function(fun=funss,geom="line",colour="red")

