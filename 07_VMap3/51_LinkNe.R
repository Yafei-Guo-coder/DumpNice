#LinkNe
setwd("/Users/guoyafei/Desktop/LinkNe/")
library(ggplot2)
data <- read.table("ne.merge.25pop", header=F,stringsAsFactors = F)
sub <- data[which(data$V1 %in% c(400,200)),]
sub_decline <- sub[which(sub$V6 %in% c("pop1","pop2","pop3","pop4","pop6","pop8","pop9","pop10","pop11","pop13","pop14","pop16","pop17","pop19","pop20","pop22","pop23")),]

sub <- data[which(data$V1 %in% c(400,200)),]
sub_no_decline <- sub[which(sub$V6 %in% c("pop5","pop7","pop12","pop15","pop18","pop21","pop24","pop25")),]
         
pdf("test1.pdf",height=30,width=10)
ggplot(sub_no_decline, aes(V1, V3,type=V6)) +
  geom_point() +
  geom_line() +
  theme_classic()+
  ylim(150,1400)+
  labs(x = 'generation', y = 'Ne')
print(p)
dev.off()

a <- paste("pop",c(9,16,17,25,2,7,20,19),sep="")
sub <- data[which(data$V3 %in% a),]
b <- paste("pop",c(9,16,17,25,2,7,20,19,15,22,1,5,8,11,6,13,4,10.21),sep="")
sub <- data[which(data$V3 %in% b),]
