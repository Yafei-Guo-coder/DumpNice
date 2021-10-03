setwd("/Users/guoyafei/Downloads/Results/")
data <- read.csv("test.filter.prune.in.27.Q",header = F, sep = " ", stringsAsFactors = F)
dim(data)
colnames(data) <- paste("Q",1:27,sep="")
decathlon2 <- read.table("data", header = T, stringsAsFactors = F)
dim(decathlon2)
names <- colnames(decathlon2)[10:755]