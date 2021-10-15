library(gggplot2)
library("ggpubr")

data <- read.table("/Users/guoyafei/Desktop/相关性.txt", header=T, stringsAsFactors = F)
sub <- data[!is.na(data$株高.郑老师),]
sub2 <- sub[!is.na(sub$气孔导度),]
ggscatter(data, x = "株高.郑老师", y = "气孔导度", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Plant height", ylab = "Stomatal conductance")
ggscatter(data, x = "ZX株高", y = "气孔导度", color = "orgCty",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Tiller", ylab = "Stomatal conductance")
ggplot(data=data, aes(x=ZX株高, y=气孔导度 ))+
  geom_point(size=3,aes(x=ZX株高, y=气孔导度,color=orgCty ))+
  stat_smooth(method="lm")+
  xlab("plant height")+
  ylab("Stomatal conductance")+
  theme_classic()

ggplot(data,aes(x = ZX株高,color=orgCty)) +
  geom_histogram(aes(x = ZX株高),stat="bin",binwidth=1, boundary = 0)+
  theme_classic()+
  xlab("plant height")

geom_text(aes(label=as.character(round(..density..,2))),stat="bin",binwidth=0.01,boundary = 0,vjust=-0.5)


