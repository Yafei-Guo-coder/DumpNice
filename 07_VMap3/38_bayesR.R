#bayesR
data <- read.table("bayesR.txt", header=F,stringsAsFactors = F)
data2 <- read.table("chrom.txt", header=T,stringsAsFactors = F)

ggplot(sub, aes(length, Va,color=type)) +
  geom_point( size=3)+
  scale_color_manual(values = c("#999999","#0072B2","#E69F00"))+
  theme_bw()
