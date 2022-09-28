#Rare Allele
setwd("/Users/guoyafei/Documents/02_VmapIII/10_RareAllele")
data <- read.table("sum2.txt", header=T, stringsAsFactors = F)

library(reshape2)
md <- melt(data, id="id")
ggplot(md, aes(x=id, y=value, color=id)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(variable~.)+
  #scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  xlab("")+ylab("")+theme_bw()+
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.position = "null",axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 10))
