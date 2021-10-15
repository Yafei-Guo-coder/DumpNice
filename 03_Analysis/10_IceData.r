#冰川数据
setwd("/Users/guoyafei/Downloads/stack/")
data <- read.table("LR04stack.txt",header=T,stringsAsFactors = F,sep="\t")
ggplot(data, aes(x=data$Time..ka., y=data$Benthic.d18O..per.mil.)) + 
  geom_line() +
  scale_x_log10()+scale_y_reverse() +
  geom_hline(aes(yintercept=4), colour="#990000", linetype="dashed")