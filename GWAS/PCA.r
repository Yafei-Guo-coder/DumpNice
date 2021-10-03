setwd("/Users/guoyafei/Documents/个人项目/傅老师/GWAS/Gapit/")
data <- read.table("PCA.txt",header=T, stringsAsFactors = F)

group <- NA
group[data$orgCty == "CHN"] <- 1
group[data$orgCty == "AM"] <- 2
group[data$orgCty == "SCA"] <- 3
group[data$orgCty == "WA"] <- 4
group[data$orgCty == "IND"] <- 5
group[data$orgCty == "JPN"] <- 6

pairs(~PC1+PC2+PC3+PC4+PC5, col = c("red", "cornflowerblue", "purple", "yellow", "black","orange")[group], upper.panel = NULL,data=data,main="PC of 329 samples",legend.plot=TRUE)
