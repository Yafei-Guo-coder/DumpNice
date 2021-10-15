setwd("/Users/guoyafei/Documents/个人项目/傅老师/20210311/GWAS_sign_gene/")
data <- read.table("VariantComponent.txt",header=T,stringsAsFactors = F)
data$X1.428374375 <- as.factor(data$X1.428374375)
data$X3.106681694 <- as.factor(data$X3.106681694)
data$X3.107363563 <- as.factor(data$X3.107363563)
data$X16.83205051 <- as.factor(data$X16.83205051)
pdf("SNP_geno_trait.pdf")
ggplot(data, aes(x = X1.428374375, y = Weight))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="NGR5_1(1-428374375)") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))
ggplot(data, aes(x = X3.106681694, y = Weight))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="PIF_2(3-106681694)") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))
#ggplot(data, aes(x = X3.107363563, y = Weight))+ 
#  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="PIF_2(3-107363563)") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))
ggplot(data, aes(x = X16.83205051, y = Weight))+ 
  geom_boxplot(position=position_dodge(0.5),width=0.6) + geom_point() +theme_classic() + labs(x="GA2ox3-B1(16-83205051)") + theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))
dev.off()
