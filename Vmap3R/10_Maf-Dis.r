library(ggplot2)
setwd("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Test_chr001")
data <- read.table("R2-Maf-Dis.txt",header=T,stringsAsFactors = F)
library(reshape)
data$R_bin <- NA
data$Maf_bin <- NA
data$Dis_bin <- NA

data[which(data$R.2 >=0 & data$R.2 < 0.5),4] <- "0<=r2<0.5"
data[which(data$R.2 >=0.5 & data$R.2 < 0.75),4] <- "0.5<=r2<0.75"
data[which(data$R.2 >=0.75 & data$R.2 < 1),4] <- "0.75<=r2<1"
data[which(data$R.2 == 1),4] <- "r2=1"

data[which(data$SNP_A_MAF >=0 & data$SNP_A_MAF < 0.05),5] <- "0-0.05"
data[which(data$SNP_A_MAF >=0.05 & data$SNP_A_MAF < 0.1),5] <- "0.05-0.1"
data[which(data$SNP_A_MAF >=0.1 & data$SNP_A_MAF < 0.15),5] <- "0.1-0.15"
data[which(data$SNP_A_MAF >=0.15 & data$SNP_A_MAF < 0.2),5] <- "0.15-0.2"
data[which(data$SNP_A_MAF >=0.2 & data$SNP_A_MAF < 0.25),5] <- "0.2-0.25"
data[which(data$SNP_A_MAF >=0.25 & data$SNP_A_MAF < 0.3),5] <- "0.25-0.3"
data[which(data$SNP_A_MAF >=0.3 & data$SNP_A_MAF < 0.35),5] <- "0.3-0.35"
data[which(data$SNP_A_MAF >=0.35 & data$SNP_A_MAF < 0.4),5] <- "0.35-0.4"
data[which(data$SNP_A_MAF >=0.4 & data$SNP_A_MAF < 0.45),5] <- "0.4-0.45"
data[which(data$SNP_A_MAF >=0.45 & data$SNP_A_MAF <= 0.5),5] <- "0.45-0.5"

data[which(data$Distance >=0 & data$Distance < 10),6] <- "0-10"
data[which(data$Distance >=10 & data$Distance < 20),6] <- "10-20"
data[which(data$Distance >=20 & data$Distance < 30),6] <- "20-30"
data[which(data$Distance >=30 & data$Distance < 40),6] <- "30-40"
data[which(data$Distance >=40 & data$Distance < 50),6] <- "40-50"
data[which(data$Distance >=50 & data$Distance < 60),6] <- "50-60"
data[which(data$Distance >=60 & data$Distance < 70),6] <- "60-70"
data[which(data$Distance >=70 & data$Distance < 80),6] <- "70-80"
data[which(data$Distance >=80 & data$Distance < 90),6] <- "80-90"
data[which(data$Distance >=90 & data$Distance <= 100),6] <- "90-100"

R <- data[,c(4,5)]
R$num <- 1
md <- melt(R,id=(c("Maf_bin","R_bin")))
md2 <- rbind(bin1[1:2000,],bin2[1:2000,],bin3[1:2000,],bin4[1:2000,],bin5[1:2000,],bin6[1:2000,],bin7,bin8,bin9)
newdata <- cast(md2,Maf_bin+R_bin~value)
colnames(newdata)[3] <- "num"


ggplot(data = newdata, mapping = aes(x = factor(Maf_bin), y = newdata$num, fill = newdata$R_bin)) + geom_bar(stat = 'identity', position = 'fill')+
  theme_classic() +
  labs(x = "Maf", y = "Propotion of SNPs(%)", title = "Allele frequency and R2 for the 10M regions with same snp numbers") + 
  theme(legend.text = element_text(size = 14),legend.title=element_blank(),plot.title = element_text(hjust = 0.5,size=20),axis.title = element_text(size = 20),axis.text.x = element_text(size = 15,vjust = 0.5, hjust = 0.5,angle = 45),axis.text.y = element_text(size = 15))


D <-  data[!is.na(data$Dis_bin),c(5,6)]
D$num <- 1
md <- melt(D,id=(c("Maf_bin","Dis_bin")))
newdata <- cast(md,Maf_bin+Dis_bin~value)
colnames(newdata)[3] <- "num"

ggplot(data = newdata, mapping = aes(x = factor(Maf_bin), y = newdata$num, fill = newdata$Dis_bin)) + geom_bar(stat = 'identity', position = 'fill')+
  theme_classic() +
  labs(x = "Maf", y = "Propotion of SNPs(%)", title = "Allele frequency and Distance of SNPs for the 10M regions") + 
  theme(legend.text = element_text(size = 14),legend.title=element_blank(),plot.title = element_text(hjust = 0.5,size=20),axis.title = element_text(size = 20),axis.text.x = element_text(size = 15,vjust = 0.5, hjust = 0.5,angle = 45),axis.text.y = element_text(size = 15))


