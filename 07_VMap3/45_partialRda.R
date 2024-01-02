#partial RDA
library(vegan)
#individual----
setwd("/Users/guoyafei/Desktop/冗余分析/样本")
phylum <- read.delim('/Users/guoyafei/Desktop/冗余分析/样本/all_20k.txt', row.names = 1, sep = '\t', header=F, stringsAsFactors = FALSE, check.names = FALSE)
colnames(phylum) <- c(1:dim(phylum)[2])
env <- read.delim('/Users/guoyafei/Desktop/冗余分析/样本/471_baypass_taxa.Info', row.names = 1, header=T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
geo <- data.frame(env[,c(1:2)])
env_all <- data.frame(env[,c(3,6:40)])
Q <- read.table("/Users/guoyafei/Desktop/冗余分析/样本/all_200k.5.Q", header=F, stringsAsFactors = F)
rownames(Q) <- rownames(phylum)
con <- cbind(env,Q)
#full model
spe.partial.rda <- rda(phylum ~ Latitude + Longitude + V1 + V2 + V3 + V4 + V5 + Elevation + prec1 + prec2 + prec3 + prec4 + prec5 + prec6 + prec7 + prec8 + soil1 + soil10+ soil11+ soil12+ soil13+ soil2 + soil3 + soil4 + soil5 + soil6 + soil7 + soil8 + soil9 + solar1+ solar2+ solar3+ temp1 + temp10+ temp11+ temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8 + temp9, data = con)
r_full <- as.numeric(RsquareAdj(spe.partial.rda)[1])
adjr_full <- as.numeric(RsquareAdj(spe.partial.rda)[2])
#p <- anova.cca(spe.partial.rda, step = 10)
#控制变量看环境
spe.partial.rda <- rda(phylum ~ Elevation + prec1 + prec2 + prec3 + prec4 + prec5 + prec6 + prec7 + prec8 + soil1 + soil10+ soil11+ soil12+ soil13+ soil2 + soil3 + soil4 + soil5 + soil6 + soil7 + soil8 + soil9 + solar1+ solar2+ solar3+ temp1 + temp10+ temp11+ temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8 + temp9 + Condition(Latitude + Longitude + V1 + V2 + V3 + V4 + V5), data = con)
r_clim <- as.numeric(RsquareAdj(spe.partial.rda)[1])
adjr_clim <- as.numeric(RsquareAdj(spe.partial.rda)[2])
#p <- anova.cca(spe.partial.rda, step = 10)
#控制变量看地理位置
spe.partial.rda <- rda(phylum ~ Latitude + Longitude + Condition(Elevation + prec1 + prec2 + prec3 + prec4 + prec5 + prec6 + prec7 + prec8 + soil1 + soil10+ soil11+ soil12+ soil13+ soil2 + soil3 + soil4 + soil5 + soil6 + soil7 + soil8 + soil9 + solar1+ solar2+ solar3+ temp1 + temp10+ temp11+ temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8 + temp9 + V1 + V2 + V3 + V4 + V5), data = con)
r_geo <- as.numeric(RsquareAdj(spe.partial.rda)[1])
adjr_geo <- as.numeric(RsquareAdj(spe.partial.rda)[2])
#p <- anova.cca(spe.partial.rda, step = 10)
#控制变量看群体结构
spe.partial.rda <- rda(phylum ~ V1 + V2 + V3 + V4 + V5 + Condition(Elevation + prec1 + prec2 + prec3 + prec4 + prec5 + prec6 + prec7 + prec8 + soil1 + soil10+ soil11+ soil12+ soil13+ soil2 + soil3 + soil4 + soil5 + soil6 + soil7 + soil8 + soil9 + solar1+ solar2+ solar3+ temp1 + temp10+ temp11+ temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8 + temp9 + Latitude + Longitude), data = con)
r_str <- as.numeric(RsquareAdj(spe.partial.rda)[1])
adjr_str <- as.numeric(RsquareAdj(spe.partial.rda)[2])
#p <- anova.cca(spe.partial.rda, step = 10)
r <- c(r_full, r_clim, r_geo, r_str)
adjr <- c(adjr_full, adjr_clim, adjr_geo, adjr_str)
out <- cbind(r,adjr)
colnames(out) <- c("r","adjr")
write.table(out, "random.r.txt", quote=F, row.names = F,sep="\t")

#cluster-------
library("data.table")
setwd("/Users/guoyafei/Desktop/冗余分析/群体")
con <- read.delim('env.10.Q.sum.txt', row.names = 1, sep = '\t', header=T, stringsAsFactors = FALSE, check.names = FALSE)
phylum_fst <- fread("all_fst.txt", header = T, stringsAsFactors = F)
rownames(phylum_fst) <- rownames(con)

#full model
spe.partial.rda <- rda(phylum ~ Latitude + Longitude + V1 + V2 + V3+ Elevation + prec1 + prec2 + prec3 + prec4 + prec5 + prec6 + prec7 + prec8 + soil1 + soil10+ soil11+ soil12+ soil13+ soil2 + soil3 + soil4 + soil5 + soil6 + soil7 + soil8 + soil9 + solar1+ solar2+ solar3+ temp1 + temp10+ temp11+ temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8 + temp9, data = con)
r_full <- as.numeric(RsquareAdj(spe.partial.rda)[1])
adjr_full <- as.numeric(RsquareAdj(spe.partial.rda)[2])
#p <- anova.cca(spe.partial.rda, step = 10)

#控制变量看环境
spe.partial.rda <- rda(phylum ~ Elevation + prec1 + prec2 + prec3 + prec4 + prec5 + prec6 + prec7 + prec8 + soil1 + soil10+ soil11+ soil12+ soil13+ soil2 + soil3 + soil4 + soil5 + soil6 + soil7 + soil8 + soil9 + solar1+ solar2+ solar3+ temp1 + temp10+ temp11+ temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8 + temp9 + Condition(Latitude + Longitude + V1 + V2 + V3 +V4 +V5+V6 + V7 + V8 + V9+V10), data = con)
r_clim <- as.numeric(RsquareAdj(spe.partial.rda)[1])
adjr_clim <- as.numeric(RsquareAdj(spe.partial.rda)[2])
#p <- anova.cca(spe.partial.rda, step = 10)

#控制变量看地理位置
spe.partial.rda <- rda(phylum ~ Latitude + Longitude + Condition(Elevation + prec1 + prec2 + prec3 + prec4 + prec5 + prec6 + prec7 + prec8 + soil1 + soil10+ soil11+ soil12+ soil13+ soil2 + soil3 + soil4 + soil5 + soil6 + soil7 + soil8 + soil9 + solar1+ solar2+ solar3+ temp1 + temp10+ temp11+ temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8 + temp9 + V1 + V2 + V3 +V4 +V5+V6 + V7 + V8 + V9+V10), data = con)
r_geo <- as.numeric(RsquareAdj(spe.partial.rda)[1])
adjr_geo <- as.numeric(RsquareAdj(spe.partial.rda)[2])
#p <- anova.cca(spe.partial.rda, step = 10)

#控制变量看群体结构
spe.partial.rda <- rda(phylum ~ Latitude + Longitude +V1 + V2 + V3 + V4 + V5 +V6 + V7 + V8 + V9+V10+ Condition(Elevation + prec1 + prec2 + prec3 + prec4 + prec5 + prec6 + prec7 + prec8 + soil1 + soil10+ soil11+ soil12+ soil13+ soil2 + soil3 + soil4 + soil5 + soil6 + soil7 + soil8 + soil9 + solar1+ solar2+ solar3+ temp1 + temp10+ temp11+ temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8 + temp9), data = con)
r_str <- as.numeric(RsquareAdj(spe.partial.rda)[1])
adjr_str <- as.numeric(RsquareAdj(spe.partial.rda)[2])
#p <- anova.cca(spe.partial.rda, step = 10)
r <- c(r_full, r_clim, r_geo, r_str)
adjr <- c(adjr_full, adjr_clim, adjr_geo, adjr_str)
out <- cbind(r,adjr)
colnames(out) <- c("r","adjr")
write.table(out, "random.r.txt", quote=F, row.names = F,sep="\t")

#画图(individual的结果)----
library(reshape2)
library(ggplot2)
setwd("/Users/guoyafei/Desktop/冗余分析/样本/V9")
data <- read.table("plot_V9.r.txt", header = T, stringsAsFactors = F)
cats <- melt(data,id="type")
sub <- cats[which(cats$type != "geog"),]
sub$type <- factor(sub$type, levels = c("full","clim","struct"))
ggplot(data=sub,mapping=aes(x=variable,y=value,fill=type))+
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.6)+
  theme_bw()





