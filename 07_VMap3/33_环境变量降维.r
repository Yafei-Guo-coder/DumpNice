
library(corrgram)
data <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)
sub <- data[,c(3,6:40)]
sub2 <- sub[,c(1,4,5,8,11,13,16,24,29,34,35,36)]
corrgram(sub, order = F,lower.panel = panel.shade,upper.panel = panel.cor,text.panel = panel.txt,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")

library(usdm)
vif(sub2)
v1 <- vifcor(sub2,th=0.7)

#环境变量的筛选
library(ggplot2)
setwd("/Users/guoyafei/Desktop/baypass")
data <- read.table("471_baypass_taxa.Info", header=T,row.names = 1,stringsAsFactors = F)
env <- data[,c(3:4,6:40)]
M <- aggregate(env, by=list(env$type),FUN=mean)
solar <- M[,25:27]
temp <- M[,c(2,28:38)]
prec <- M[,4:11]
soil <- M[,c(12:24)]

solar <- env[,24:26]
temp <- env[,c(1,27:37)]
prec <- env[,3:10]
soil <- env[,c(11:23)]

#PCA
library(FactoMineR)
library(factoextra)
library(psych)

#查看loading
PC_temp <- principal(temp, nfactors=3,rotate="none")
PC_solar <- principal(solar, nfactors=3,rotate="none")
PC_prec <- principal(prec, nfactors=3,rotate="none")
PC_soil <- principal(soil, nfactors=3,rotate="none")

#个体相关性：
temp1(0.99);temp6(0.92);temp10(0.90);temp11(0.90)
solar1(0.99);solar3(0.94)
prec1(0.97);prec5(0.84);prec6(0.83);prec3(0.81);prec2(0.80)
soil9(0.94);soil12(0.93)

#群体相关性：
PC1:temp11(0.98);temp1(0.96);temp6(0.96) 
solar3(0.99);solar1(0.94)
prec1(0.98);prec5(0.85);prec6(0.85);prec3(0.83);prec2(0.80)
soil9(0.94);soil12(0.94);soil13(0.93);soil11(-0.83)

sub <- solar
dt <- as.matrix(scale(sub))
rm1<- cor(dt)
rs1<- eigen(rm1)
val <- rs1$values
U <- as.matrix(rs1$vectors)
PC <- dt %*% U
PC <- as.data.frame(PC)
colnames(PC)[1:3] <- paste("PC",1:3,sep="")
write.table(PC,"/Users/guoyafei/Desktop/GF/等位基因频率变化/solar-PC.txt",quote=F,row.names = F)
#cor.test(temp$temp1,temp$PC1)




M <- aggregate(sub, by=list(sub$type),FUN=mean)
sub$type <- factor(sub$type, levels=rownames(M[order(M$PC1),]),ordered = TRUE)
#soil
sub$type <- factor(sub$type, levels=c("10", "19", "21", "6",  "4",  "15", "14", "23", "2",  "9",  "16", "18", "5",  "17", "8",  "1",  "25", "22", "11", "3",  "13", "20", "7","12"),ordered = TRUE)
ggplot(sub2, aes(x=sub2$type, y = sub2$PC1)) +
    geom_boxplot(color="blue",
                 fill="blue",
                 alpha=0.2) +
    ylab("soil$PC1") +
    xlab("Group")+
    theme_bw()




Proportion_of_Variance <- val/sum(val)
Cumulative_Proportion <- cumsum(Proportion_of_Variance)
#碎石图绘制;
par(mar=c(6,6,2,2))
plot(rs1$values,type="b",
     cex=2,
     cex.lab=2,
     cex.axis=2,
     lty=2,
     lwd=2,
     xlab = "PC",
     ylab="Eigenvalue (Principal Component Variance)")


df<-data.frame(PC, dt_all$ploidy, dt_all$Region)
library(ggplot2)

xlab<-paste0("PC1(",round(Proportion_of_Variance[1]*100,2),"%)")
ylab<-paste0("PC2(",round(Proportion_of_Variance[2]*100,2),"%)")

p1<-ggplot(data = PC,aes(x=PC1,y=PC2))+
  geom_point()+labs(x=xlab,y=ylab,color="")+
  theme_classic()

