#snpAttr
#关联分析#############
setwd("/Users/guoyafei/Desktop/snpAttr")
data <- read.table("rand100k_max.txt", header=T, stringsAsFactors = F)
library(ggplot2)
# [1] "ID"                        "fst_pop01_pop02"           "fst_pop01_pop03"           "fst_pop01_pop04"           "fst_pop01_pop05"          
# [6] "fst_pop01_pop06"           "fst_pop02_pop03"           "fst_pop02_pop04"           "fst_pop02_pop05"           "fst_pop02_pop06"          
#[11] "fst_pop03_pop04"           "fst_pop03_pop05"           "fst_pop03_pop06"           "fst_pop04_pop05"           "fst_pop04_pop06"          
#[16] "fst_pop05_pop06"           "baypass_type1_BF"          "baypass_type1_Beta"        "baypass_type1_BetaSD"      "baypass_type1_eBPis"      
#[21] "envGWAS_type1_Beta_median" "envGWAS_type1_Beta_absmax" "envGWAS_type1_P_median"    "envGWAS_type1_P_min"       "baypass_type2_BF"         
#[26] "baypass_type2_Beta"        "baypass_type2_BetaSD"      "baypass_type2_eBPis"       "envGWAS_type2_Beta_median" "envGWAS_type2_Beta_absmax"
#[31] "envGWAS_type2_P_median"    "envGWAS_type2_P_min"       "baypass_type3_BF"          "baypass_type3_Beta"        "baypass_type3_BetaSD"     
#[36] "baypass_type3_eBPis"       "envGWAS_type3_Beta_median" "envGWAS_type3_Beta_absmax" "envGWAS_type3_P_median"    "envGWAS_type3_P_min"      
#[41] "baypass_type4_BF"          "baypass_type4_Beta"        "baypass_type4_BetaSD"      "baypass_type4_eBPis"       "envGWAS_type4_Beta_median"
#[46] "envGWAS_type4_Beta_absmax" "envGWAS_type4_P_median"    "envGWAS_type4_P_min"       "poppair01_3sstdev"         "poppair02_3sstdev"        
#[51] "poppair03_3sstdev"         "poppair04_3sstdev"         "poppair05_3sstdev"         "Average"                   "Median"                   
#[56] "Min"                       "Max"
sub <- data[which(data$baypass_type1_Beta != "NA" & data$envGWAS_type1_P_min < 0.01),]
sub$selec <- 0
sub[which(sub$poppair01_3sstdev == "1"), 58] <- "1"
sub[which(sub$poppair02_3sstdev == "1"), 58] <- "1"
sub[which(sub$poppair03_3sstdev == "1"), 58] <- "1"
sub[which(sub$poppair04_3sstdev == "1"), 58] <- "1"
sub[which(sub$poppair05_3sstdev == "1"), 58] <- "1"
sub1 <- sub[which(sub$poppair02_3sstdev == "1"),]
pdf("envgwas.pdf")
ggplot(data=sub, aes(x=abs(sub$envGWAS_type1_Beta_absmax), y=Max)) +
  geom_point(size=0.5, color="#bebada")+
  #geom_line()+
  theme_classic()+
  xlab("envGWAS BETA")+
  ylim(0,1)+
  ylab("Population Max Fst")+
  #geom_point(data=sub2,aes(x=abs(sub2$envGWAS_type1_Beta_absmax), y=sub2$Max), color="red",size=2)+
  geom_smooth(method=lm , color="red", fill="#bebada", se=TRUE) +
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
dev.off()

pdf("baypass.pdf")
ggplot(data=sub, aes(x=abs(sub$baypass_type1_Beta), y=Max)) +
  geom_point(size=0.5, color="#bebada")+
  #geom_line()+
  theme_classic()+
  xlab("Baypass BETA")+
  ylim(0,1)+
  ylab("Population Max Fst")+
  #geom_point(data=sub2,aes(x=abs(sub2$envGWAS_type1_Beta_absmax), y=sub2$Max), color="red",size=2)+
  geom_smooth(method=lm , color="red", fill="#bebada", se=TRUE) +
  theme(legend.position="none",plot.title = element_text(color="red", size=20, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))
dev.off()



ggplot(sub, aes(y = abs(envGWAS_type4_Beta_absmax), group = as.factor(poppair05_3sstdev), fill = as.factor(poppair05_3sstdev)))+
  geom_boxplot(notch = FALSE) +
  theme_classic() +  
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
  ) 



#画图baypass_Beta-----
sub <- data[which(data$baypass_type1_Beta != "NA" & data$envGWAS_type1_P_min < 0.05),]
a1 <- t.test(abs(sub[which(sub$poppair01_3sstdev == "1"), 18]), abs(sub[which(sub$poppair01_3sstdev == "0"), 18]))
b1 <- t.test(abs(sub[which(sub$poppair02_3sstdev == "1"), 18]), abs(sub[which(sub$poppair02_3sstdev == "0"), 18]))
c1 <- t.test(abs(sub[which(sub$poppair03_3sstdev == "1"), 18]), abs(sub[which(sub$poppair03_3sstdev == "0"), 18]))
d1 <- t.test(abs(sub[which(sub$poppair04_3sstdev == "1"), 18]), abs(sub[which(sub$poppair04_3sstdev == "0"), 18]))
e1 <- t.test(abs(sub[which(sub$poppair05_3sstdev == "1"), 18]), abs(sub[which(sub$poppair05_3sstdev == "0"), 18]))


sub <- data[which(data$baypass_type2_Beta != "NA" & data$envGWAS_type2_P_min < 0.05),]
a2 <- t.test(abs(sub[which(sub$poppair01_3sstdev == "1"), 26]), abs(sub[which(sub$poppair01_3sstdev == "0"), 26]))
b2 <- t.test(abs(sub[which(sub$poppair02_3sstdev == "1"), 26]), abs(sub[which(sub$poppair02_3sstdev == "0"), 26]))
c2 <- t.test(abs(sub[which(sub$poppair03_3sstdev == "1"), 26]), abs(sub[which(sub$poppair03_3sstdev == "0"), 26]))
d2 <- t.test(abs(sub[which(sub$poppair04_3sstdev == "1"), 26]), abs(sub[which(sub$poppair04_3sstdev == "0"), 26]))
e2 <- t.test(abs(sub[which(sub$poppair05_3sstdev == "1"), 26]), abs(sub[which(sub$poppair05_3sstdev == "0"), 26]))


sub <- data[which(data$baypass_type3_Beta != "NA" & data$envGWAS_type3_P_min < 0.05),]
a3 <- t.test(abs(sub[which(sub$poppair01_3sstdev == "1"), 34]), abs(sub[which(sub$poppair01_3sstdev == "0"), 34]))
b3 <- t.test(abs(sub[which(sub$poppair02_3sstdev == "1"), 34]), abs(sub[which(sub$poppair02_3sstdev == "0"), 34]))
c3 <- t.test(abs(sub[which(sub$poppair03_3sstdev == "1"), 34]), abs(sub[which(sub$poppair03_3sstdev == "0"), 34]))
d3 <- t.test(abs(sub[which(sub$poppair04_3sstdev == "1"), 34]), abs(sub[which(sub$poppair04_3sstdev == "0"), 34]))
e3 <- t.test(abs(sub[which(sub$poppair05_3sstdev == "1"), 34]), abs(sub[which(sub$poppair05_3sstdev == "0"), 34]))


sub <- data[which(data$baypass_type4_Beta != "NA" & data$envGWAS_type4_P_min < 0.05),]
a4 <- t.test(abs(sub[which(sub$poppair01_3sstdev == "1"), 42]), abs(sub[which(sub$poppair01_3sstdev == "0"), 42]))
b4 <- t.test(abs(sub[which(sub$poppair02_3sstdev == "1"), 42]), abs(sub[which(sub$poppair02_3sstdev == "0"), 42]))
c4 <- t.test(abs(sub[which(sub$poppair03_3sstdev == "1"), 42]), abs(sub[which(sub$poppair03_3sstdev == "0"), 42]))
d4 <- t.test(abs(sub[which(sub$poppair04_3sstdev == "1"), 42]), abs(sub[which(sub$poppair04_3sstdev == "0"), 42]))
e4 <- t.test(abs(sub[which(sub$poppair05_3sstdev == "1"), 42]), abs(sub[which(sub$poppair05_3sstdev == "0"), 42]))



p <- c(a1$p.value, b1$p.value,c1$p.value,d1$p.value,e1$p.value,a2$p.value, b2$p.value,c2$p.value,d2$p.value,e2$p.value,a3$p.value, b3$p.value,c3$p.value,d3$p.value,e3$p.value,a4$p.value, b4$p.value,c4$p.value,d4$p.value,e4$p.value)
select_mean <- c(a1$estimate[1], b1$estimate[1],c1$estimate[1],d1$estimate[1],e1$estimate[1],a2$estimate[1], b2$estimate[1],c2$estimate[1],d2$estimate[1],e2$estimate[1],a3$estimate[1], b3$estimate[1],c3$estimate[1],d3$estimate[1],e3$estimate[1],a4$estimate[1], b4$estimate[1],c4$estimate[1],d4$estimate[1],e4$estimate[1])
noselect_mean <- c(a2$estimate[2], b2$estimate[2],c2$estimate[2],d2$estimate[2],e2$estimate[2],a2$estimate[2], b2$estimate[2],c2$estimate[2],d2$estimate[2],e2$estimate[2],a3$estimate[2], b3$estimate[2],c3$estimate[2],d3$estimate[2],e3$estimate[2],a4$estimate[2], b4$estimate[2],c4$estimate[2],d4$estimate[2],e4$estimate[2])

all <- as.data.frame(cbind(p,select_mean,noselect_mean))
row.names(all) <- paste(rep(paste("type",seq(1:4),sep=""),each=5), rep(paste("route",seq(1:5),sep=""),4), sep="_")
all$type <- rep(paste("type",seq(1:4),sep=""),each=5)
all$route <- rep(paste("route",seq(1:5),sep=""),4)
all$logP <- -log10(all$p)
all$color <- 1

all[which(all$select_mean < all$noselect_mean),7] <- -1

all$size <- all$noselect_mean/all$select_mean

a <- all[,c(4,5,6)]
cats <- cast(a,route~type)
row.names(cats) <- cats[,1]
P <- cats[,c(2:5)]

b <- all[,c(4,5,8)]
size <- cast(b,route~type)
row.names(size) <- size[,1]
P <- size[,c(2:5)]

#100k
all[4,6] <- 20
all[14,6] <- 22
#500k
all[4,6] <- 70
all[14,6] <- 75
#all$logP
#[1]   3.2737316   3.8265883  26.8891561 207.1771299  62.2093594  29.9767688
#[7]   0.1843917  31.2065247  48.9082474  11.2716850   0.9271739   6.5335728
#[13]   7.5826685 299.0038437  58.3358235  14.2686604   2.8236534  67.0373461
#[19]  62.2731295  45.5857438

#p select_mean noselect_mean  type  route      logP color
#type2_route2 0.6540460  0.08297624    0.08383898 type2 route2 0.1843917    -1
#type3_route1 0.1182568  0.06654072    0.06846354 type3 route1 0.9271739    -1

pdf("baypass_500k.pdf", width=4,height=4)
ggplot( data=all, aes(x=type, y=route, size = size, color = logP,shape=as.factor(color))) +
  geom_point() +
  scale_size(range = c(2, 19), name="Population (M)") +
  theme_bw()+
  theme(legend.position="none")
dev.off()

pdf("baypass_500k_2.pdf", width=10,height=10)
ggplot( data=all, aes(x=type, y=route, size = size, color = logP,shape=as.factor(color))) +
  geom_point() +
  scale_size(range = c(2, 19), name="Population (M)") +
  theme_bw()
dev.off()
#画图envGWAS_Beta------
library(reshape)
library(ggplot2)
library(dplyr)
library(plotly)
library(viridis)
library(hrbrthemes)
sub <- data[which(data$baypass_type1_Beta != "NA" & data$envGWAS_type1_P_min < 0.5),]
a1 <- t.test(abs(sub[which(sub$poppair01_3sstdev == "1"), 22]), abs(sub[which(sub$poppair01_3sstdev == "0"), 22]))
b1 <- t.test(abs(sub[which(sub$poppair02_3sstdev == "1"), 22]), abs(sub[which(sub$poppair02_3sstdev == "0"), 22]))
c1 <- t.test(abs(sub[which(sub$poppair03_3sstdev == "1"), 22]), abs(sub[which(sub$poppair03_3sstdev == "0"), 22]))
d1 <- t.test(abs(sub[which(sub$poppair04_3sstdev == "1"), 22]), abs(sub[which(sub$poppair04_3sstdev == "0"), 22]))
e1 <- t.test(abs(sub[which(sub$poppair05_3sstdev == "1"), 22]), abs(sub[which(sub$poppair05_3sstdev == "0"), 22]))


sub <- data[which(data$baypass_type2_Beta != "NA" & data$envGWAS_type2_P_min < 0.5),]
a2 <- t.test(abs(sub[which(sub$poppair01_3sstdev == "1"), 30]), abs(sub[which(sub$poppair01_3sstdev == "0"), 30]))
b2 <- t.test(abs(sub[which(sub$poppair02_3sstdev == "1"), 30]), abs(sub[which(sub$poppair02_3sstdev == "0"), 30]))
c2 <- t.test(abs(sub[which(sub$poppair03_3sstdev == "1"), 30]), abs(sub[which(sub$poppair03_3sstdev == "0"), 30]))
d2 <- t.test(abs(sub[which(sub$poppair04_3sstdev == "1"), 30]), abs(sub[which(sub$poppair04_3sstdev == "0"), 30]))
e2 <- t.test(abs(sub[which(sub$poppair05_3sstdev == "1"), 30]), abs(sub[which(sub$poppair05_3sstdev == "0"), 30]))


sub <- data[which(data$baypass_type3_Beta != "NA" & data$envGWAS_type3_P_min < 0.5),]
a3 <- t.test(abs(sub[which(sub$poppair01_3sstdev == "1"), 38]), abs(sub[which(sub$poppair01_3sstdev == "0"), 38]))
b3 <- t.test(abs(sub[which(sub$poppair02_3sstdev == "1"), 38]), abs(sub[which(sub$poppair02_3sstdev == "0"), 38]))
c3 <- t.test(abs(sub[which(sub$poppair03_3sstdev == "1"), 38]), abs(sub[which(sub$poppair03_3sstdev == "0"), 38]))
d3 <- t.test(abs(sub[which(sub$poppair04_3sstdev == "1"), 38]), abs(sub[which(sub$poppair04_3sstdev == "0"), 38]))
e3 <- t.test(abs(sub[which(sub$poppair05_3sstdev == "1"), 38]), abs(sub[which(sub$poppair05_3sstdev == "0"), 38]))


sub <- data[which(data$baypass_type4_Beta != "NA" & data$envGWAS_type4_P_min < 0.5),]
a4 <- t.test(abs(sub[which(sub$poppair01_3sstdev == "1"), 46]), abs(sub[which(sub$poppair01_3sstdev == "0"), 46]))
b4 <- t.test(abs(sub[which(sub$poppair02_3sstdev == "1"), 46]), abs(sub[which(sub$poppair02_3sstdev == "0"), 46]))
c4 <- t.test(abs(sub[which(sub$poppair03_3sstdev == "1"), 46]), abs(sub[which(sub$poppair03_3sstdev == "0"), 46]))
d4 <- t.test(abs(sub[which(sub$poppair04_3sstdev == "1"), 46]), abs(sub[which(sub$poppair04_3sstdev == "0"), 46]))
e4 <- t.test(abs(sub[which(sub$poppair05_3sstdev == "1"), 46]), abs(sub[which(sub$poppair05_3sstdev == "0"), 46]))



p <- c(a1$p.value, b1$p.value,c1$p.value,d1$p.value,e1$p.value,a2$p.value, b2$p.value,c2$p.value,d2$p.value,e2$p.value,a3$p.value, b3$p.value,c3$p.value,d3$p.value,e3$p.value,a4$p.value, b4$p.value,c4$p.value,d4$p.value,e4$p.value)
select_mean <- c(a1$estimate[1], b1$estimate[1],c1$estimate[1],d1$estimate[1],e1$estimate[1],a2$estimate[1], b2$estimate[1],c2$estimate[1],d2$estimate[1],e2$estimate[1],a3$estimate[1], b3$estimate[1],c3$estimate[1],d3$estimate[1],e3$estimate[1],a4$estimate[1], b4$estimate[1],c4$estimate[1],d4$estimate[1],e4$estimate[1])
noselect_mean <- c(a2$estimate[2], b2$estimate[2],c2$estimate[2],d2$estimate[2],e2$estimate[2],a2$estimate[2], b2$estimate[2],c2$estimate[2],d2$estimate[2],e2$estimate[2],a3$estimate[2], b3$estimate[2],c3$estimate[2],d3$estimate[2],e3$estimate[2],a4$estimate[2], b4$estimate[2],c4$estimate[2],d4$estimate[2],e4$estimate[2])

all <- as.data.frame(cbind(p,select_mean,noselect_mean))
row.names(all) <- paste(rep(paste("type",seq(1:4),sep=""),each=5), rep(paste("route",seq(1:5),sep=""),4), sep="_")
all$type <- rep(paste("type",seq(1:4),sep=""),each=5)
all$route <- rep(paste("route",seq(1:5),sep=""),4)
all$logP <- -log10(all$p)
all$color <- 1

all[which(all$select_mean < all$noselect_mean),7] <- -1

all$size <- all$noselect_mean/all$select_mean


a <- all[,c(4,5,6)]
cats <- cast(a,route~type)
row.names(cats) <- cats[,1]
P <- cats[,c(2:5)]

b <- all[,c(4,5,8)]
size <- cast(b,route~type)
row.names(size) <- size[,1]
P <- size[,c(2:5)]
#100k
all[14,6] <- 15
#1000k
all[7,6] <- 70
all[15,6] <- 72
all[18,6] <- 75
all[10,6] <- 82
all[14,6] <- 86
#all$logP
#1.2090868  6.7190869 19.5329749 67.5561362  5.1916852  
#0.1408561 70.0000000  1.1023750  6.9580246 82.0000000  
#0.7817175  2.4052346 55.5332300 86.0000000 72.0000000  
#0.7129150  1.0026799 75.0000000 23.9448461 26.8344640
pdf("envGWAS_1000k_2.pdf", width=10,height=10)
ggplot( data=all, aes(x=type, y=route, size = size, color = logP,shape=as.factor(color))) +
  geom_point() +
  scale_size(range = c(2, 19), name="Population (M)") +
  theme_bw()
dev.off()

#画一下整体的分布
sub1 <- data[which(data$baypass_type1_Beta != "NA" & data$envGWAS_type1_P_min < 0.5),]
sub2 <- data[which(data$baypass_type2_Beta != "NA" & data$envGWAS_type1_P_min < 0.5),]
sub3 <- data[which(data$baypass_type3_Beta != "NA" & data$envGWAS_type1_P_min < 0.5),]
sub4 <- data[which(data$baypass_type4_Beta != "NA" & data$envGWAS_type1_P_min < 0.5),]
all <- rbind(sub1,sub2,sub3,sub4)
all$select <- 0
all[which(all$poppair01_3sstdev == 1 | all$poppair02_3sstdev == 1 |all$poppair03_3sstdev == 1 | all$poppair04_3sstdev == 1 |all$poppair05_3sstdev == 1 ), 58] <- 1

sub <- all[,c(18,26,34,42,57,58)]
sub$beta <- apply(sub[,1:4], 1, function(x) max(abs(x), na.rm = TRUE))
t.test(sub[which(sub$select == "1"), 7], sub[which(sub$select == "0"), 7])
pdf("baypass_all.pdf",height=5,width=6)
ggplot(sub, aes(y = beta, group = as.factor(select), fill = as.factor(select)))+
  geom_boxplot(notch = T) +
  theme_classic() +  
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
  ) 

sub <- all[,c(21,22,29,30,37,38,45,46,57,58)]
#sub$beta_median <- apply(sub[,c(1,3,5,7)], 1, function(x) max(abs(x), na.rm = TRUE))
sub$beta_max <- apply(sub[,c(2,4,6,8)], 1, function(x) max(abs(x), na.rm = TRUE))
t.test(sub[which(sub$select == "1"), 11], sub[which(sub$select == "0"), 11])
pdf("envgwas_all.pdf",height=5,width=6)
ggplot(sub, aes(y = beta_max, group = as.factor(select), fill = as.factor(select)))+
  geom_boxplot(notch = T) +
  theme_classic() +  
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
  ) 
dev.off()

##2018Science#############
setwd("/Users/guoyafei/Desktop/network/2018science")
data <- read.table("selection_Uniq.txt", header=F,stringsAsFactors = F)


data$logP <- -log(data$V3,10)
data$color <- 1
data[which(data$logP > 2),6] <- 0
data$ID <- paste(data$V1,data$V4,sep="_")
#Share
data$ID <- factor(data$ID,levels=c("module1_root","module15_root", "module30_root","module66_root","module8_spike","module15_spike","module21_spike", "module24_spike","module33_spike","module40_spike", "module48_spike","module3_leaf", "module26_leaf","module35_leaf","module6_grain", "module7_grain","module31_grain","module37_grain","module47_grain","module63_grain","module47_disease","module60_disease","module7_abiotic","module9_abiotic","module10_abiotic","module32_abiotic","module41_abiotic","module55_abiotic","module56_abiotic"))
#Unique
data$ID <- factor(data$ID,levels=c("module1_root","module8_root","module15_root","module1_spike","module7_spike","module21_spike","module24_spike", "module25_spike", "module1_leaf","module6_leaf","module21_leaf","module27_leaf", "module16_grain","module31_grain","module45_grain","module9_disease", "module10_disease", "module17_disease", "module46_disease", "module1_abiotic","module7_abiotic","module57_abiotic","module67_abiotic","module74_abiotic", "module100_abiotic"))
pdf("baypass_500k.pdf", width=4,height=4)
ggplot( data=data, aes(x=ID, y=V2, size=5,color = logP,shape=as.factor(color))) +
  geom_point() +
  scale_size(range = c(12,12), name="Population (M)") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  #theme(legend.position="none")
dev.off()

#不同类型基因的unique和shared情况####
setwd("/Users/guoyafei/Desktop/network/2018science/gene")
library(UpSetR)
type1 <- read.table("type1.3k.gene", header=F, stringsAsFactors = F)[,1]
type2 <- read.table("type2.3k.gene", header=F, stringsAsFactors = F)[,1]
type3 <- read.table("type3.3k.gene", header=F, stringsAsFactors = F)[,1]
type4 <- read.table("type4.3k.gene", header=F, stringsAsFactors = F)[,1]
type1 <- read.table("type1.pos.bed", header=F, stringsAsFactors = F)
type1_pos <- paste(type1$V1, type1$V3, sep="-")
type2 <- read.table("type2.pos.bed", header=F, stringsAsFactors = F)
type2_pos <- paste(type2$V1, type2$V3, sep="-")
type3 <- read.table("type3.pos.bed", header=F, stringsAsFactors = F)
type3_pos <- paste(type3$V1, type3$V3, sep="-")
type4 <- read.table("type4.pos.bed", header=F, stringsAsFactors = F)
type4_pos <- paste(type4$V1, type4$V3, sep="-")
#gene
a <- list(solar<-type1,temperature<-type2,precipitation<-type3,soil<-type4)
#pos
a <- list(solar<-type1_pos,temperature<-type2_pos,precipitation<-type3_pos,soil<-type4_pos)
names(a) <- c("solar","temperature","precipitation","soil")
tmp <- names(a)
p <- upset(fromList(a), keep.order = TRUE, sets=tmp, text.scale = c(2),point.size = 2.5, line.size = 1.5)

#分析不同的位点unique和shared类型的效应值有什么差异####

setwd("/Users/guoyafei/Desktop/network/2018science/gene")
dataM <- read.table("/Users/guoyafei/Desktop/network/2018science/gene/all.pos.timeM.out.txt", header=F, stringsAsFactors = F)
dataU <- read.table("/Users/guoyafei/Desktop/network/2018science/gene/all.pos.timeU.out.txt", header=F, stringsAsFactors = F)

M_type1 <- dataM[which(dataM$V17 != "NA"),c(18,22,57)]
U_type1 <- dataU[which(dataU$V17 != "NA"),c(18,22,57)]
M_type1$type <- "M"
U_type1$type <- "U"
colnames(M_type1) <- c("baypass_beta","envGWAS_beta","max","type")
colnames(U_type1) <- c("baypass_beta","envGWAS_beta","max","type")
all <- rbind(M_type1,U_type1)

M_type1 <- dataM[which(dataM$V25 != "NA"),c(26,30,57)]
U_type1 <- dataU[which(dataU$V25 != "NA"),c(26,30,57)]
M_type1$type <- "M"
U_type1$type <- "U"
colnames(M_type1) <- c("baypass_beta","envGWAS_beta","max","type")
colnames(U_type1) <- c("baypass_beta","envGWAS_beta","max","type")
all2 <- rbind(M_type1,U_type1)


M_type1 <- dataM[which(dataM$V33 != "NA"),c(34,38,57)]
U_type1 <- dataU[which(dataU$V33 != "NA"),c(34,38,57)]
M_type1$type <- "M"
U_type1$type <- "U"
colnames(M_type1) <- c("baypass_beta","envGWAS_beta","max","type")
colnames(U_type1) <- c("baypass_beta","envGWAS_beta","max","type")
all3 <- rbind(M_type1,U_type1)


M_type1 <- dataM[which(dataM$V41 != "NA"),c(42,46,57)]
U_type1 <- dataU[which(dataU$V41 != "NA"),c(42,46,57)]
M_type1$type <- "M"
U_type1$type <- "U"
colnames(M_type1) <- c("baypass_beta","envGWAS_beta","max","type")
colnames(U_type1) <- c("baypass_beta","envGWAS_beta","max","type")
all4 <- rbind(M_type1,U_type1)

all_f <- rbind(all,all2,all3,all4)
M <- all_f[sample(rownames(all_f[which(all_f$type == "M"),]),2000,replace = FALSE),]
U <- all_f[sample(rownames(all_f[which(all_f$type == "U"),]),2000,replace = FALSE),]
out <- rbind(M,U)

ggplot(out, aes(y =abs(envGWAS_beta), group=as.factor(type), fill=as.factor(type)))+
    geom_boxplot( notch = F) +
    theme_classic() + 
    scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
    theme(legend.title="")+
    theme(
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y=element_text(size=15),
      axis.text.x=element_blank(),
      axis.title.y=element_text(size = 15),
      axis.title.x=element_text(size = 15),
      legend.title=element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(color="Black", size=20, face="bold.italic")
    )
ggplot(out, aes(y =abs(max), group=as.factor(type), fill=as.factor(type)))+
  geom_boxplot( notch = F) +
  theme_classic() + 
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
  theme(legend.title="")+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
ggplot(all_f, aes(x=abs(envGWAS_beta),y =abs(max), color=as.factor(type),group=as.factor(type), fill=as.factor(type)))+
  geom_point(alpha=0.5) +
  theme_classic() + 
  scale_fill_manual(values = c("#ffffb3","#fccde5")) +
  scale_color_manual(values = c("#FDC086","#BEAED4")) +
  theme(legend.title="")+
  #geom_smooth(method = "lm", fill = "lightgray") +
  geom_density_2d()+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )


#分析不同的位点unique和shared类型在网络连接度上有什么差异####
setwd("/Users/guoyafei/Desktop/network/2018science/连接度")
library(ggplot2)
file <- c("correlation_grain.tpm.module16","correlation_grain.tpm.module31","correlation_grain.tpm.module37","correlation_grain.tpm.module45","correlation_grain.tpm.module47","correlation_grain.tpm.module7","correlation_leaf.tpm.module21","correlation_leaf.tpm.module26","correlation_leaf.tpm.module27","correlation_leaf.tpm.module35","correlation_root.tpm.module15","correlation_root.tpm.module30","correlation_root.tpm.module66","correlation_spike.tpm.module21","correlation_spike.tpm.module24","correlation_spike.tpm.module25","correlation_spike.tpm.module33","correlation_spike.tpm.module40","correlation_spike.tpm.module48","correlation_spike.tpm.module7","correlation_spike.tpm.module8")
file <- c("correlation_grain.tpm.module31","correlation_grain.tpm.module45","correlation_leaf.tpm.module21","correlation_leaf.tpm.module27","correlation_root.tpm.module66","correlation_spike.tpm.module21","correlation_spike.tpm.module33","correlation_spike.tpm.module40")
file <- c("correlation_grain.tpm.module16","correlation_grain.tpm.module37","correlation_grain.tpm.module45","correlation_grain.tpm.module47","correlation_grain.tpm.module7","correlation_spike.tpm.module25","correlation_spike.tpm.module33","correlation_spike.tpm.module7")
nw <- paste(file,".r03.out",sep="")
id <- paste(file,".ID",sep="")

type1_U <- read.table("type1.3k.time1.txt", header=F, stringsAsFactors = F)
type1_M <- read.table("type1.3k.timeM.txt", header=F, stringsAsFactors = F)
type2_U <- read.table("type2.3k.time1.txt", header=F, stringsAsFactors = F)
type2_M <- read.table("type2.3k.timeM.txt", header=F, stringsAsFactors = F)
type3_U <- read.table("type3.3k.time1.txt", header=F, stringsAsFactors = F)
type3_M <- read.table("type3.3k.timeM.txt", header=F, stringsAsFactors = F)
type4_U <- read.table("type4.3k.time1.txt", header=F, stringsAsFactors = F)
type4_M <- read.table("type4.3k.timeM.txt", header=F, stringsAsFactors = F)


pdf("spread_r05_介数连接度.pdf")
for(i in c(1:length(file))){
  network <- read.table(nw[i], header=F, check.names = F,stringsAsFactors = F)
  ID <- read.table(id[i], header=F, check.names = F,stringsAsFactors = F)
  ID$type <- "NO"
  ID[which(ID$V2 %in% type2_U[,1]),3] <- "U"
  ID[which(ID$V2 %in% type2_M[,1]),3] <- "M"
  #ID[which(ID$V2 %in% type4_B[,1]),3] <- "B"
  network$type <- "Z"
  network[which(network$V1 %in% ID[which(ID$type == "U"),1]),5] <- "U"
  network[which(network$V1 %in% ID[which(ID$type == "M"),1]),5] <- "M"
  #network[which(network$V1 %in% ID[which(ID$type == "B"),1]),5] <- "B"
  p <- ggplot(network, aes(y =V3, group=as.factor(type), fill=as.factor(type)))+
    geom_boxplot( notch = F) +
    theme_classic() + 
    scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
    theme(legend.title="")+
    ggtitle(file[i])+
    theme(
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y=element_text(size=15),
      axis.text.x=element_blank(),
      axis.title.y=element_text(size = 15),
      axis.title.x=element_text(size = 15),
      legend.title=element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(color="Black", size=20, face="bold.italic")
    )
  print(p)
}
dev.off()

t.test(network[which(network$type == "M"),3], network[which(network$type == "U"),3])

wilcox.test(network[which(network$type == "M"),3], network[which(network$type == "U"),3])


t.test(sub[which(sub$type == "M"),4], sub[which(sub$type == "U"),4])
wilcox.test(sub[which(sub$type == "M"),4], sub[which(sub$type == "U"),4])
table(sub$type)


#分析不同的位点unique和shared类型在网络中表达强度上有什么差异#####
expre <- read.table("/Users/guoyafei/Desktop/network/2018science/expression/grain.tpm.type2.expression", header=F,stringsAsFactors = F)
expre <- expre[which(expre$V2 > 0.05),]
shuf_expre <- read.table("/Users/guoyafei/Desktop/network/2018science/expression/grain.shuf1k.tpm.expression", header=F,stringsAsFactors = F)
shuf_expre <- shuf_expre[which(shuf_expre$V2 > 0.05),]
network$gene <- NA
rownames(ID) <- ID[,1]


all <- merge(network, ID, by.x="V1", by.y="V1")
all2 <- merge(all,expre,by.x="V2.y", by.y="V1")
colnames(all2) <- c("gene","V1","V2","V3","V4","V5","V6","V7","V8","V9")
shuf_expre$V4 <- NA
shuf_expre$V5 <- NA
shuf_expre$V6 <- NA
shuf_expre$V7 <- NA
shuf_expre$V8 <- "Z"
shuf_expre$V9 <- "Z"
shuf_expre$V10 <- NA
shuf_expre$V11 <- NA
shuf_expre2 <- shuf_expre[,c(1,4,5,6,7,8,9,10,2,3)]
colnames(shuf_expre2) <- c("gene","V1","V2","V3","V4","V5","V6","V7","V8","V9")
all3 <- rbind(all2,shuf_expre2)

sub <- all3[which(all3$V8 < 50),]
ggplot(sub, aes(y =V8, group=as.factor(V5), fill=as.factor(V5)))+
  geom_boxplot( notch = F) +
  theme_classic() + 
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
  theme(legend.title="")+
  #ggtitle(file[i])+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )

t.test(all3[which(all3$V5 == "U"),9], all3[which(all3$V5 == "M"),9])
t.test(sub[which(sub$V5 == "Z"),9], sub[which(sub$V5 == "M"),9])



wilcox.test(all3[which(all3$V5 == "M"),10], all3[which(all3$V5 == "U"),10])
table(expre$type)


setwd("/Users/guoyafei/Desktop/network/2018science/expression/")
library(ggplot2)

expre <- read.table("/Users/guoyafei/Desktop/network/2018science/expression/grain.tpm.type1.expression", header=F,stringsAsFactors = F)

type1_U <- read.table("../连接度/type1.3k.time1.txt", header=F, stringsAsFactors = F)
type1_M <- read.table("../连接度/type1.3k.timeM.txt", header=F, stringsAsFactors = F)
type2_U <- read.table("../连接度/type2.3k.time1.txt", header=F, stringsAsFactors = F)
type2_M <- read.table("../连接度/type2.3k.timeM.txt", header=F, stringsAsFactors = F)
type3_U <- read.table("../连接度/type3.3k.time1.txt", header=F, stringsAsFactors = F)
type3_M <- read.table("../连接度/type3.3k.timeM.txt", header=F, stringsAsFactors = F)
type4_U <- read.table("../连接度/type4.3k.time1.txt", header=F, stringsAsFactors = F)
type4_M <- read.table("../连接度/type4.3k.timeM.txt", header=F, stringsAsFactors = F)

expre$type <- NA
expre[which(expre$V1 %in% type1_U[,1]),4] <- "U"
expre[which(expre$V1 %in% type1_M[,1]),4] <- "M"

sub <- expre[which(expre$V2 < 25),]
ggplot(sub, aes(y =V2, group=as.factor(type), fill=as.factor(type)))+
    geom_boxplot( notch = F) +
    theme_classic() + 
    scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
    theme(legend.title="")+
    #ggtitle(file[i])+
    theme(
      axis.line.y = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y=element_text(size=15),
      axis.text.x=element_blank(),
      axis.title.y=element_text(size = 15),
      axis.title.x=element_text(size = 15),
      legend.title=element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(color="Black", size=20, face="bold.italic")
    )

t.test(expre[which(expre$type == "M"),2], expre[which(expre$type == "U"),2])
wilcox.test(expre[which(expre$type == "M"),2], expre[which(expre$type == "U"),2])
table(expre$type)



#分析不同的位点unique和shared类型在不同组织中表达量上有什么差异####
setwd("/Users/guoyafei/Desktop/network/2018science/expression")
library(ggplot2)
library(stringr)
data <- read.table("expression.summary", header=F,stringsAsFactors = F)
data <- read.table("expression.merge.summary", header=F,stringsAsFactors = F)
data <- read.table("expression.merge3.summary", header=F,stringsAsFactors = F)
data$type <- unlist(lapply(X = data$V4, FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))
data$tissue <- unlist(lapply(X = data$V4, FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))
data$tissue <- factor(data$tissue,levels =c("root","spike","leaf","grain","disease","abiotic"))

data$class <- NA
sub <- data[which(data$V1 == "U" | data$V1 == "M"), ] 
sub <- data[which(data$V1 == "A"), ] 
ggplot(data=data, mapping=aes(x = type, y = V2, fill=V1)) +
  scale_fill_manual(values = c("#A07D35","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00"))+ 
  geom_bar(stat="identity",position=position_dodge(0.85),width=0.75) +
  geom_errorbar(aes(ymax = V2+V3, ymin = V2-V3),
                position = position_dodge(0.85), width = 0.1)+
  #facet_grid(tissue~.)+
  theme_classic()
  
ggplot(data=data, mapping=aes(x = tissue, y = V2, fill=V1)) +
  geom_bar(stat="identity",position=position_dodge(0.85),width=0.75) +
  facet_grid(type~.)+
  theme_classic()

ggplot(data=data, mapping=aes(x = tissue, y = V2, fill=type)) +
  scale_fill_manual(values = c("grey","#A07D35","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00"))+ 
  geom_bar(stat="identity",position=position_dodge(0.85),width=0.75) +
  geom_errorbar(aes(ymax = V2+V3, ymin = V2-V3),
                position = position_dodge(0.85), width = 0.1)+
  #facet_grid(type~.)+
  theme_classic()

ggplot(data=sub, mapping=aes(x = V2, y = type, fill=tissue)) +
  geom_bar(stat="identity",position=position_dodge(0.85),width=0.75) +
  scale_fill_manual(values = c("#A07D35","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00"))+ 
  geom_errorbar(aes(ymax = V2+V3, ymin = V2-V3),
                position = position_dodge(0.85), width = 0.1)+
  #facet_grid(V1~.)+
  #geom_errorbar(aes(ymin=(V2-sqrt(V4)),ymax=(V2+sqrt(V4)),width=2,size=2))+
  theme_classic()



ggplot(data, aes(x = V1,y = V2, group=V1))+
  geom_bar()+ 
  #geom_errorbar(aes(ymin=(V2-V5),ymax=(V2+V5)),width=2,size=2)+
  #scale_fill_manual(values = c("red","blue"))+  
  labs(x = "Taxa Number",y = "Snp proportion") +
  facet_grid(type~tissue)+
  theme_bw() +
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  )+
  theme(legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 10))
#分析不同基因的表达广度####
type1_U <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/geneV2/type1.time1.gene",header=F,stringsAsFactors = F)
type1_M <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/geneV2/type1.timeM.gene",header=F,stringsAsFactors = F)
type2_U <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/geneV2/type2.time1.gene",header=F,stringsAsFactors = F)
type2_M <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/geneV2/type2.timeM.gene",header=F,stringsAsFactors = F)
type3_U <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/geneV2/type3.time1.gene",header=F,stringsAsFactors = F)
type3_M <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/geneV2/type3.timeM.gene",header=F,stringsAsFactors = F)
type4_U <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/geneV2/type4.time1.gene",header=F,stringsAsFactors = F)
type4_M <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/geneV2/type4.timeM.gene",header=F,stringsAsFactors = F)

all_U <- rbind(type1_U,type2_U,type3_U,type4_U)
all_U <- all_U[!duplicated(all_U),]
all_M <- rbind(type1_M,type2_M,type3_M,type4_M)
all_M <- all_M[!duplicated(all_M),]

type1_U <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/type1.3k.time1.txt", header=F, stringsAsFactors = F)
type1_2 <- read.table("type1.3k.time2.txt", header=F, stringsAsFactors = F)
type1_3 <- read.table("type1.3k.time3.txt", header=F, stringsAsFactors = F)
type1_4 <- read.table("type1.3k.time4.txt", header=F, stringsAsFactors = F)
type1_M <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/type1.3k.timeM.txt", header=F, stringsAsFactors = F)

type2_U <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/type2.3k.time1.txt", header=F, stringsAsFactors = F)
type2_2 <- read.table("type2.3k.time2.txt", header=F, stringsAsFactors = F)
type2_3 <- read.table("type2.3k.time3.txt", header=F, stringsAsFactors = F)
type2_4 <- read.table("type2.3k.time4.txt", header=F, stringsAsFactors = F)
type2_M <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/type2.3k.timeM.txt", header=F, stringsAsFactors = F)

type3_U <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/type3.3k.time1.txt", header=F, stringsAsFactors = F)
type3_2 <- read.table("type3.3k.time2.txt", header=F, stringsAsFactors = F)
type3_3 <- read.table("type3.3k.time3.txt", header=F, stringsAsFactors = F)
type3_4 <- read.table("type3.3k.time4.txt", header=F, stringsAsFactors = F)
type3_M <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/type3.3k.timeM.txt", header=F, stringsAsFactors = F)

type4_U <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/type4.3k.time1.txt", header=F, stringsAsFactors = F)
type4_2 <- read.table("type4.3k.time2.txt", header=F, stringsAsFactors = F)
type4_3 <- read.table("type4.3k.time3.txt", header=F, stringsAsFactors = F)
type4_4 <- read.table("type4.3k.time4.txt", header=F, stringsAsFactors = F)
type4_M <- read.table("/Users/guoyafei/Desktop/network/2018science/连接度/type4.3k.timeM.txt", header=F, stringsAsFactors = F)

all_U <- read.table("all.3k.timeU.txt", header=F, stringsAsFactors = F)
all_M <- read.table("all.3k.timeM.txt", header=F, stringsAsFactors = F)

select_U <- read.table("poppair.3k.timeU.gene", header=F, stringsAsFactors = F)
select_M <- read.table("poppair.3k.timeM.gene", header=F, stringsAsFactors = F)


setwd("/Users/guoyafei/Desktop/network/2018science/expression")
data1 <- read.table("type1.tau",header=F,stringsAsFactors = F)
data2 <- read.table("type2.tau",header=F,stringsAsFactors = F)
data3 <- read.table("type3.tau",header=F,stringsAsFactors = F)
data4 <- read.table("type4.tau",header=F,stringsAsFactors = F)
data <- rbind(data1,data2,data3,data4)
data <- data[!duplicated(data),]

shuf <- read.table("shuf10k.tau", header=F,stringsAsFactors = F)

rownames(data) <-data[,1]
#画图展示tau值
data_o <- data[order(data$V10), ]
list<-seq(1,2400,by=20)
sub <- as.matrix(data_o[list,6:9])
tau <- data_o[list,]
tau$V1 <- factor(tau$V1,levels = tau$V1)
ggplot(tau,aes(V10,V1))+geom_bar(stat="identity",fill="steelblue")+theme_classic()
ggplot(tau,aes(V1,V10))+ theme_bw()+ geom_histogram()
#画图展示U和M基因tau分布的差异
data$type <- "Z"
data[which(data$V1 %in% type1_U[,1]),11] <- "1"
data[which(data$V1 %in% type3_2[,1]),11] <- "2"
data[which(data$V1 %in% type3_3[,1]),11] <- "3"
data[which(data$V1 %in% type1_4[,1]),11] <- "4"

data$type <- "Z"
data[which(data$V1 %in% type1_U[,1]),11] <- "U"
data[which(data$V1 %in% type1_M[,1]),11] <- "M"

data$type <- "Z"
data[which(data$V1 %in% all_U[,1]),11] <- "U"
data[which(data$V1 %in% all_M[,1]),11] <- "M"

shuf$type <- "Z"
all <- rbind(data,shuf)

data[which(data$V1 %in% all_U[,1]),11] <- "U"
data[which(data$V1 %in% all_M[,1]),11] <- "M"

data[which(data$V1 %in% select_U[,1]),11] <- "U"
data[which(data$V1 %in% select_M[,1]),11] <- "M"

t.test(data[which(data$type=="Z"),10],data[which(data$type=="U"),10])

wilcox.test(data[which(data$type=="1"),10],data[which(data$type=="2"),10])

ggplot(all, aes(y =V10, group=as.factor(type), fill=as.factor(type)))+
  geom_boxplot( notch = F) +
  theme_classic() +   
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue")) +
  theme(legend.title="")+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
 


#分析不同基因的网络连接度有什么差别#######
setwd("/Users/guoyafei/Desktop/network/2018science/degree")
data <- read.table("gene.degree", header=F,stringsAsFactors = F)
#abiotic <- data[which(data$V4 == "grain"),]
#abiotic <- data
ggplot(data, aes(y =V5, group=as.factor(V2), fill=as.factor(V2)))+
  geom_boxplot( notch = F) +
  theme_classic() +   
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue")) +
  theme(legend.title="")+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
t.test(abiotic[which(abiotic$V2=="adapt"),5],abiotic[which(abiotic$V2=="shuf"),5])
wilcox.test(abiotic[which(abiotic$V2=="adapt"),5],abiotic[which(abiotic$V2=="shuf"),5])

data <- read.table("fourType_picMin_repeat_degree.txt", header=T, stringsAsFactors = F)
data <- read.table("test4", header=T, stringsAsFactors = F)
data <- read.table("tmp", header=T, stringsAsFactors = F)
data <- data[which(data$V9 > 0),]
data <- data[which(data$V4 < 0.05),]
b <- read.table("shuf.cons.txt", header=F,stringsAsFactors = F)
type1 <- data[which(data$Class == "adapt"),]
[1] "Gene"                "Type"                "Class"               "p"                   "q"                   "n_est"              
[7] "abiotic"             "disease"             "root"                "spike"               "leaf"                "grain"              
[13] "tau"                 "cons"                "Gene.1"              "fst_Average"         "fst_Median"          "fst_Min"            
[19] "fst_Max"             "baypass_BF"          "baypass_Beta"        "baypass_eBPis"       "envGWAS_Beta_median" "envGWAS_Beta_absmax"
[25] "envGWAS_P_median"    "envGWAS_P_min"       "fst_pop01_pop02"     "fst_pop02_pop03"     "fst_pop03_pop04"     "fst_pop03_pop05"    
[31] "fst_pop05_pop06"  

> colnames(data)
[1] "Gene"                "Type"                "Class"               "p"                   "q"                  
[6] "n_est"               "abiotic"             "disease"             "root"                "spike"              
[11] "leaf"                "grain"               "tau"                 "cons"                "Gene.1"             
[16] "fst_Average"         "fst_Median"          "fst_Min"             "fst_Max"             "baypass_BF"         
[21] "baypass_Beta"        "baypass_eBPis"       "envGWAS_Beta_median" "envGWAS_Beta_absmax" "envGWAS_P_median"   
[26] "envGWAS_P_min"  

summary(lm(data$fst_pop05_pop06 ~ data$cons))
summary(lm(abs(data$baypass_Beta) ~ abs(data$envGWAS_Beta_absmax)))
summary(lm(abs(data[,23])~data[,4]))
ggplot(data, aes(x =leaf, y=fst_pop01_pop02, color = as.factor(Type.1),group=as.factor(Type.1), fill=as.factor(Type.1)))+
#ggplot(data, aes(x = leaf, y=abs(envGWAS_Beta_absmax)))+
  geom_point() +
  facet_grid(Type.1~.)+
  theme_classic() +
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue")) +
  theme(legend.title="")+
  geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
summary(lm(data$p~data$cons))
ggplot(data, aes(x = p, y=cons))+
  geom_point() +
  #facet_grid(V9~.)+
  theme_classic() +
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue")) +
  theme(legend.title="")+
  geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
ggplot(data, aes(x = fst_pop02_pop03))+
  geom_histogram( notch = F) +
  #facet_grid(Type~.)+
  theme_classic() +   
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue")) +
  theme(legend.title="")+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )

#ggplot(data, aes(y =fst_pop01_pop02, group=as.factor(Type.1), fill=as.factor(Type.1)))+
ggplot(data, aes(y =cons, group=as.factor(Type), fill=as.factor(Type)))+
  geom_boxplot(notch = F) +
  theme_classic() +   
  #scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue")) +
  theme(legend.title="")+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )

#t.test([,4],type1[,9])
wilcox.test(type1[which(type1$V2=="adapt"),5],type1[which(type1$V2=="shuf"),5])
summary(lm(data[,13]~data[,14]))
t.test(data[which(data$V9 == "2"),10],data[which(data$V9 == "3"),10])
t.test(data[,7],data[,10])
t.test(data[which(data$Type.1 == "TimeU"),29],data[which(data$Type.1 == "TimeM"),29])

#meta分析####
setwd("/Users/guoyafei/Desktop/snpAttr/meta")

data <- read.table("type1.3k.gene.out", header=T, stringsAsFactors = F)

#计算tau, degree, fst_max, cons1与repeatability相关性热图，区分成4个模块
#clustering_distance_rows = 'euclidean', # 计算聚类间距的算法，可选'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
#clustering_method = 'complete', # 聚类方法, 可选'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'

#在使用 pheatmap 和手动执行层次聚类的方法中，结果可能不完全相同，这取决于多个因素，包括聚类方法、距离计算方法和聚类的细节处理。具体来说：
#聚类方法：pheatmap 默认使用 Ward 的方法 (ward.D2) 进行聚类，而我在手动聚类的示例中使用了“完全链接”（complete linkage）。这些方法在聚类合并步骤中使用的准则不同，因此可能导致不同的聚类结果。
#距离度量：在 pheatmap 的默认设置中，使用的是欧几里得距离（euclidean distance）。在手动示例中，我同样使用了欧几里得距离，所以在距离度量上应该是一致的。但是，用户可以在两种方法中选择不同的距离度量方式，这也可能导致结果差异。
#随机性：一些聚类方法可能涉及到随机初始状态的选择，尤其是在大规模数据集中或者使用一些特定的优化聚类算法时。但在标准层次聚类中，这通常不是一个因素。
#数值稳定性和计算精度：实现上的细微差异可能会导致结果的轻微变化，特别是在浮点运算中。
data <- read.table("all.adapt.txt", header=T, stringsAsFactors = F)
#cons1
sub <- t(scale(na.omit(data[,c(2,4)])))
pheatmap(sub,cluster_rows = T,cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "complete", border_color = "white",cutree_cols  = 4)
hc <- hclust(dist(t(sub)), method = "complete")
clusters <- cutree(hc, 4)
group1 <- data[which(clusters == 1),c(1,2,4)]
group1$type <- "group1"
group2 <- data[which(clusters == 2),c(1,2,4)]
group2$type <- "group2"
group3 <- data[which(clusters == 3),c(1,2,4)]
group3$type <- "group3"
group4 <- data[which(clusters == 4),c(1,2,4)]
group4$type <- "group4"
all <- rbind(group1,group2,group3,group4)
write.table(all, "cons1.4group.txt", row.names = F,quote=F, sep="\t")
#tau
sub <- t(scale(na.omit(data[,c(3,4)])))
pheatmap(sub,cluster_rows = T,cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "complete", border_color = "white",cutree_cols  = 4)
hc <- hclust(dist(t(sub)), method = "complete")
clusters <- cutree(hc, 4)
group1 <- data[which(clusters == 1),c(1,3,4)]
group1$type <- "group1"
group2 <- data[which(clusters == 2),c(1,3,4)]
group2$type <- "group2"
group3 <- data[which(clusters == 3),c(1,3,4)]
group3$type <- "group3"
group4 <- data[which(clusters == 4),c(1,3,4)]
group4$type <- "group4"
all <- rbind(group1,group2,group3,group4)
write.table(all, "tau.4group.txt", row.names = F,quote=F, sep="\t")
#root
sub <- t(scale(na.omit(data[,c(5,4)])))
pheatmap(sub,cluster_rows = T,cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "complete", border_color = "white",cutree_cols  = 4)
hc <- hclust(dist(t(sub)), method = "complete")
clusters <- cutree(hc, 4)
group1 <- data[which(clusters == 1),c(1,5,4)]
group1$type <- "group1"
group2 <- data[which(clusters == 2),c(1,5,4)]
group2$type <- "group2"
group3 <- data[which(clusters == 3),c(1,5,4)]
group3$type <- "group3"
group4 <- data[which(clusters == 4),c(1,5,4)]
group4$type <- "group4"
all <- rbind(group1,group2,group3,group4)
write.table(all, "rootdegree.4group.txt", row.names = F,quote=F, sep="\t")


data <- read.table("shuf.5k.gene.out", header=T, stringsAsFactors = F)
data1 <- read.table("type2.3k.time1.txt.out", header=T,strings=F)
data1$cat <- "uniq"
data2 <- read.table("type1.3k.timeM.txt.out", header=T,strings=F)
data2$cat <- "share"
data <- as.data.frame(rbind(data1,data2))
data[which(data$root_cc == 0),10] <- NA
sub <- na.omit(data[,c(2,4,8,5,21)])

sub2 <- na.omit(data[,c(2,4,8,5)])

sub$envGWAS_Beta_absmax <- abs(sub$envGWAS_Beta_absmax)
sub$baypass_Beta <- abs(sub$baypass_Beta)

data1 <- read.table("shuf.5k.gene.out", header=T, stringsAsFactors = F)
data1$type <- "shuf"

data2 <- read.table("type1.3k.gene.out", header=T, stringsAsFactors = F)
data2$type <- "adapt"

all <- as.data.frame(rbind(data1[,c(2,17)],data2[,c(2,34)]))

sub$envGWAS_Beta_absmax <- abs(sub$envGWAS_Beta_absmax)
sub$baypass_Beta <- abs(sub$baypass_Beta)

#colnames(data)
[1] "Gene"                "cons1"               "cons2"               "tau"                 "xp_p"                "poppair"            
[7] "popsingle"           "root"                "root_bc"             "root_cc"             "grain"               "grain_bc"           
[13] "grain_cc"            "abiotic"             "abiotic_bc"          "abiotic_cc"          "ID"                  "fst_Average"        
[19] "fst_Median"          "fst_Min"             "fst_Max"             "baypass_BF"          "baypass_Beta"        "baypass_eBPis"      
[25] "envGWAS_Beta_median" "envGWAS_Beta_absmax" "envGWAS_P_median"    "envGWAS_P_min"       "fst_pop01_pop02"     "fst_pop02_pop03"    
[31] "fst_pop03_pop04"     "fst_pop03_pop05"     "fst_pop05_pop06"

#colnames(data)
[1] "Gene"       "cons1"      "cons2"      "tau"        "xp_p"       "poppair"    "popsingle"  "root"       "root_bc"    "root_cc"    "grain"     
[12] "grain_bc"   "grain_cc"   "abiotic"    "abiotic_bc" "abiotic_cc"

library(corrgram)
corrgram(sub, order = F,lower.panel = panel.shade,upper.panel=panel.pie,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")
library(pheatmap)
data <- read.table("test", header=F,stringsAsFactors = F)
sub <- t(scale(na.omit(data[,c(2,4)])))

pheatmap(sub,cluster_rows = T,cluster_cols = T, border_color = "white",cutree_cols  = 4)

#相关性
library(GGally)
ggcorr(sub, method = c("everything","pearson"),  
       geom = "circle", max_size = 15,  # 使用圆形表示相关系数
       label = TRUE, label_size = 3,angle = 0,  # 设置圆形的角度
       palette = "RdYlBu",  # 设置调色板为红黄蓝
       hjust = 1, size = 4, color = "grey50",  # 设置相关系数标签的位置、大小和颜色
       layout.exp = 0.5) +  # 设置相关系数标签的名称为ρ
  #scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +  # 设置alpha值的映射关系，当系数满足条件时设置透明度为0.25，否则为0
  guides(alpha = FALSE) +  # 不显示alpha的图例
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),  # 设置图的边距
        legend.background = element_blank(),  # 不显示图例的背景
        legend.spacing.x = unit(0, "cm"))  # 设置图例水平间距为0cm

#3D散点图
library(plotly)
a <- as.data.frame(scale(sub))
m <- plot_ly(a, x = ~root, y = ~tau, z = ~date,
             #colors = c("#FF6DAE","#D4CA3A","#00BDFF"),
             marker = list(size = 2)) %>%
  add_markers()

add_lines(
  # plots one line per city since p knows city is a grouping variable
  add_lines(m, alpha = 0.2, name = "Texan Cities", hoverinfo = "none"),
  name = "Houston", data = filter(txhousing, city == "Houston")
)

layer_city <- function(plot, name) {
  plot %>% filter(city == name) %>% add_lines(name = name)
}

layer_iqr <- function(plot) {
  plot %>%
    group_by(date) %>% 
    summarise(
      q1 = quantile(median, 0.25, na.rm = TRUE),
      m = median(median, na.rm = TRUE),
      q3 = quantile(median, 0.75, na.rm = TRUE)
    ) %>%
    add_lines(y = ~m, name = "median", color = I("black")) %>%
    add_ribbons(ymin = ~q1, ymax = ~q3, name = "IQR", color = I("black"))
}

m %>%
  add_fun(layer_iqr)

allCities <- txhousing %>%
  group_by(city) %>%
  plot_ly(x = ~date, y = ~median) %>%
  add_lines(alpha = 0.2, name = "Texan Cities", hoverinfo = "none")

allCities %>%
  filter(city == "Houston") %>%
  add_lines(name = "Houston")


ggplot(data, aes(x=root, y = abiotic))+
  #ggplot(data, aes(x = leaf, y=abs(envGWAS_Beta_absmax)))+
  geom_point() +
  theme_classic() +
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue")) +
  theme(legend.title="")+
  geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
summary(lm(data$root ~ data$grain))

[1] "Gene"                "cons1"               "cons2"               "tau"                 "xp_p"               
[6] "poppair"             "popsingle"           "root"                "root_bc"             "root_cc"            
[11] "grain"               "grain_bc"            "grain_cc"            "abiotic"             "abiotic_bc"         
[16] "abiotic_cc"          "ID"                  "fst_Average"         "fst_Median"          "fst_Min"            
[21] "fst_Max"             "baypass_BF"          "baypass_Beta"        "baypass_eBPis"       "envGWAS_Beta_median"
[26] "envGWAS_Beta_absmax" "envGWAS_P_median"    "envGWAS_P_min"       "fst_pop01_pop02"     "fst_pop02_pop03"    
[31] "fst_pop03_pop04"     "fst_pop03_pop05"     "fst_pop05_pop06"     "type"  

all <- as.data.frame(rbind(data1[,c(9,17)],data2[,c(9,34)]))

#看选择重复性和选择次数的相关性
rep <- read.table("genome_selectTime_PicMin.txt", header=F, stringsAsFactors = F)
ggplot(all, aes(y = xp_p, group=as.factor(type),color =as.factor(type)))+
  geom_boxplot() +
  #facet_grid(Type~.)+
  theme_classic() +   
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue","#FF6DAE","#D4CA3A","#00BDFF")) +
  theme(legend.title="")+
  theme(
    panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )

#meta分析NEW#####
setwd("/Users/guoyafei/Desktop/snpAttr/meta")

data <- read.table("type1.3k.gene.out", header=T, stringsAsFactors = F)
data <- read.table("all.adapt.txt", header=T, stringsAsFactors = F)
#cons1
sub <- t(scale(na.omit(data[,c(2,4)])))
pheatmap(sub,cluster_rows = T,cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "complete", border_color = "white",cutree_cols  = 4)
hc <- hclust(dist(t(sub)), method = "complete")
clusters <- cutree(hc, 4)
group1 <- data[which(clusters == 1),c(1,2,4)]
group1$type <- "group1"
group2 <- data[which(clusters == 2),c(1,2,4)]
group2$type <- "group2"
group3 <- data[which(clusters == 3),c(1,2,4)]
group3$type <- "group3"
group4 <- data[which(clusters == 4),c(1,2,4)]
group4$type <- "group4"
all <- rbind(group1,group2,group3,group4)
write.table(all, "cons1.4group.txt", row.names = F,quote=F, sep="\t")
#tau
sub <- t(scale(na.omit(data[,c(3,4)])))
pheatmap(sub,cluster_rows = T,cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "complete", border_color = "white",cutree_cols  = 4)
hc <- hclust(dist(t(sub)), method = "complete")
clusters <- cutree(hc, 4)
group1 <- data[which(clusters == 1),c(1,3,4)]
group1$type <- "group1"
group2 <- data[which(clusters == 2),c(1,3,4)]
group2$type <- "group2"
group3 <- data[which(clusters == 3),c(1,3,4)]
group3$type <- "group3"
group4 <- data[which(clusters == 4),c(1,3,4)]
group4$type <- "group4"
all <- rbind(group1,group2,group3,group4)
write.table(all, "tau.4group.txt", row.names = F,quote=F, sep="\t")
#root
sub <- t(scale(na.omit(data[,c(5,4)])))
pheatmap(sub,cluster_rows = T,cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "complete", border_color = "white",cutree_cols  = 4)
hc <- hclust(dist(t(sub)), method = "complete")
clusters <- cutree(hc, 4)
group1 <- data[which(clusters == 1),c(1,5,4)]
group1$type <- "group1"
group2 <- data[which(clusters == 2),c(1,5,4)]
group2$type <- "group2"
group3 <- data[which(clusters == 3),c(1,5,4)]
group3$type <- "group3"
group4 <- data[which(clusters == 4),c(1,5,4)]
group4$type <- "group4"
all <- rbind(group1,group2,group3,group4)
write.table(all, "rootdegree.4group.txt", row.names = F,quote=F, sep="\t")


data <- read.table("shuf.5k.gene.out", header=T, stringsAsFactors = F)
data1 <- read.table("type2.3k.time1.txt.out", header=T,strings=F)
data1$cat <- "uniq"
data2 <- read.table("type1.3k.timeM.txt.out", header=T,strings=F)
data2$cat <- "share"
data <- as.data.frame(rbind(data1,data2))
data[which(data$root_cc == 0),10] <- NA
sub <- na.omit(data[,c(2,4,8,5,21)])

sub2 <- na.omit(data[,c(2,4,8,5)])

sub$envGWAS_Beta_absmax <- abs(sub$envGWAS_Beta_absmax)
sub$baypass_Beta <- abs(sub$baypass_Beta)

data1 <- read.table("shuf.5k.gene.out", header=T, stringsAsFactors = F)
data1$type <- "shuf"

data2 <- read.table("type1.3k.gene.out", header=T, stringsAsFactors = F)
data2$type <- "adapt"

all <- as.data.frame(rbind(data1[,c(2,17)],data2[,c(2,34)]))

sub$envGWAS_Beta_absmax <- abs(sub$envGWAS_Beta_absmax)
sub$baypass_Beta <- abs(sub$baypass_Beta)

#colnames(data)
[1] "Gene"                "cons1"               "cons2"               "tau"                 "xp_p"                "poppair"            
[7] "popsingle"           "root"                "root_bc"             "root_cc"             "grain"               "grain_bc"           
[13] "grain_cc"            "abiotic"             "abiotic_bc"          "abiotic_cc"          "ID"                  "fst_Average"        
[19] "fst_Median"          "fst_Min"             "fst_Max"             "baypass_BF"          "baypass_Beta"        "baypass_eBPis"      
[25] "envGWAS_Beta_median" "envGWAS_Beta_absmax" "envGWAS_P_median"    "envGWAS_P_min"       "fst_pop01_pop02"     "fst_pop02_pop03"    
[31] "fst_pop03_pop04"     "fst_pop03_pop05"     "fst_pop05_pop06"

#colnames(data)
[1] "Gene"       "cons1"      "cons2"      "tau"        "xp_p"       "poppair"    "popsingle"  "root"       "root_bc"    "root_cc"    "grain"     
[12] "grain_bc"   "grain_cc"   "abiotic"    "abiotic_bc" "abiotic_cc"

library(corrgram)
corrgram(sub, order = F,lower.panel = panel.shade,upper.panel=panel.pie,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")
library(pheatmap)
data <- read.table("test", header=F,stringsAsFactors = F)
sub <- t(scale(na.omit(data[,c(2,4)])))

pheatmap(sub,cluster_rows = T,cluster_cols = T, border_color = "white",cutree_cols  = 4)

#相关性
library(GGally)
ggcorr(sub, method = c("everything","pearson"),  
       geom = "circle", max_size = 15,  # 使用圆形表示相关系数
       label = TRUE, label_size = 3,angle = 0,  # 设置圆形的角度
       palette = "RdYlBu",  # 设置调色板为红黄蓝
       hjust = 1, size = 4, color = "grey50",  # 设置相关系数标签的位置、大小和颜色
       layout.exp = 0.5) +  # 设置相关系数标签的名称为ρ
  #scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +  # 设置alpha值的映射关系，当系数满足条件时设置透明度为0.25，否则为0
  guides(alpha = FALSE) +  # 不显示alpha的图例
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),  # 设置图的边距
        legend.background = element_blank(),  # 不显示图例的背景
        legend.spacing.x = unit(0, "cm"))  # 设置图例水平间距为0cm

#3D散点图
library(plotly)
a <- as.data.frame(scale(sub))
m <- plot_ly(a, x = ~root, y = ~tau, z = ~date,
             #colors = c("#FF6DAE","#D4CA3A","#00BDFF"),
             marker = list(size = 2)) %>%
  add_markers()

add_lines(
  # plots one line per city since p knows city is a grouping variable
  add_lines(m, alpha = 0.2, name = "Texan Cities", hoverinfo = "none"),
  name = "Houston", data = filter(txhousing, city == "Houston")
)

layer_city <- function(plot, name) {
  plot %>% filter(city == name) %>% add_lines(name = name)
}

layer_iqr <- function(plot) {
  plot %>%
    group_by(date) %>% 
    summarise(
      q1 = quantile(median, 0.25, na.rm = TRUE),
      m = median(median, na.rm = TRUE),
      q3 = quantile(median, 0.75, na.rm = TRUE)
    ) %>%
    add_lines(y = ~m, name = "median", color = I("black")) %>%
    add_ribbons(ymin = ~q1, ymax = ~q3, name = "IQR", color = I("black"))
}

m %>%
  add_fun(layer_iqr)

allCities <- txhousing %>%
  group_by(city) %>%
  plot_ly(x = ~date, y = ~median) %>%
  add_lines(alpha = 0.2, name = "Texan Cities", hoverinfo = "none")

allCities %>%
  filter(city == "Houston") %>%
  add_lines(name = "Houston")


ggplot(data, aes(x=root, y = abiotic))+
  #ggplot(data, aes(x = leaf, y=abs(envGWAS_Beta_absmax)))+
  geom_point() +
  theme_classic() +
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue")) +
  theme(legend.title="")+
  geom_smooth(se = TRUE, method = "gam", formula = y ~ x)+
  theme(
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )
summary(lm(data$root ~ data$grain))

[1] "Gene"                "cons1"               "cons2"               "tau"                 "xp_p"               
[6] "poppair"             "popsingle"           "root"                "root_bc"             "root_cc"            
[11] "grain"               "grain_bc"            "grain_cc"            "abiotic"             "abiotic_bc"         
[16] "abiotic_cc"          "ID"                  "fst_Average"         "fst_Median"          "fst_Min"            
[21] "fst_Max"             "baypass_BF"          "baypass_Beta"        "baypass_eBPis"       "envGWAS_Beta_median"
[26] "envGWAS_Beta_absmax" "envGWAS_P_median"    "envGWAS_P_min"       "fst_pop01_pop02"     "fst_pop02_pop03"    
[31] "fst_pop03_pop04"     "fst_pop03_pop05"     "fst_pop05_pop06"     "type"  

all <- as.data.frame(rbind(data1[,c(9,17)],data2[,c(9,34)]))

#看选择重复性和选择次数的相关性
rep <- read.table("genome_selectTime_PicMin.txt", header=F, stringsAsFactors = F)
ggplot(all, aes(y = xp_p, group=as.factor(type),color =as.factor(type)))+
  geom_boxplot() +
  #facet_grid(Type~.)+
  theme_classic() +   
  scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5","blue","#FF6DAE","#D4CA3A","#00BDFF")) +
  theme(legend.title="")+
  theme(
    panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_text(size=15),
    axis.text.x=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color="Black", size=20, face="bold.italic")
  )



