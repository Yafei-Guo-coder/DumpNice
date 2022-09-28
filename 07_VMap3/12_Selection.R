#Fst
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/Fst")
domemmer <- c("domemmer.sub1_sub2.windowed.weir.fst", "domemmer.sub1_sub3.windowed.weir.fst", "domemmer.sub1_sub4.windowed.weir.fst", "domemmer.sub1_sub5.windowed.weir.fst", "domemmer.sub2_sub3.windowed.weir.fst", "domemmer.sub2_sub4.windowed.weir.fst", "domemmer.sub2_sub5.windowed.weir.fst", "domemmer.sub3_sub4.windowed.weir.fst", "domemmer.sub3_sub5.windowed.weir.fst", "domemmer.sub4_sub5.windowed.weir.fst")
freethresh <- c("freethresh.sub1_sub2.windowed.weir.fst", "freethresh.sub1_sub3.windowed.weir.fst", "freethresh.sub1_sub4.windowed.weir.fst", "freethresh.sub2_sub3.windowed.weir.fst", "freethresh.sub2_sub4.windowed.weir.fst", "freethresh.sub3_sub4.windowed.weir.fst")
wildemmer <- c("wildemmer.sub1_sub2.windowed.weir.fst", "wildemmer.sub1_sub3.windowed.weir.fst", "wildemmer.sub1_sub4.windowed.weir.fst", "wildemmer.sub1_sub5.windowed.weir.fst", "wildemmer.sub1_sub6.windowed.weir.fst", "wildemmer.sub1_sub7.windowed.weir.fst", "wildemmer.sub2_sub3.windowed.weir.fst", "wildemmer.sub2_sub4.windowed.weir.fst", "wildemmer.sub2_sub5.windowed.weir.fst", "wildemmer.sub2_sub6.windowed.weir.fst", "wildemmer.sub2_sub7.windowed.weir.fst", "wildemmer.sub3_sub4.windowed.weir.fst", "wildemmer.sub3_sub5.windowed.weir.fst", "wildemmer.sub3_sub6.windowed.weir.fst", "wildemmer.sub3_sub7.windowed.weir.fst", "wildemmer.sub4_sub5.windowed.weir.fst", "wildemmer.sub4_sub6.windowed.weir.fst", "wildemmer.sub4_sub7.windowed.weir.fst", "wildemmer.sub5_sub6.windowed.weir.fst", "wildemmer.sub5_sub7.windowed.weir.fst", "wildemmer.sub6_sub7.windowed.weir.fst")
x=vector()
for (i in wildemmer) {
  data <- read.table(i, header = T,stringsAsFactors = F)
  x <- append(x,mean(data$MEAN_FST))
}

#PI
library(ggplot2)
library(RColorBrewer)
color <- brewer.pal(12, "Paired")[c(1,2,8,7)]
color <- brewer.pal(12, "Paired")[c(2,8)]
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/Pi")
data1 <- read.table("sub.cultivar-wildemmer.A.sites.pi", header = T, stringsAsFactors = F)
sub1 <- data1[which(data1$cultivar.wildemmer  < 2 ),,drop=F]
data2 <- read.table("sub.landrace-wildemmer.A.sites.pi", header = T, stringsAsFactors = F)
sub2 <- data2[which(data2$landrace.wildemmer  < 2 ),,drop=F]
data3 <- read.table("sub.domemmer-wildemmer.A.sites.pi", header = T, stringsAsFactors = F)
sub3 <- data3[which(data3$domemmer.wildemmer  < 2 ),,drop=F]
data4 <- read.table("sub.freethresh-wildemmer.A.sites.pi", header = T, stringsAsFactors = F)
sub4 <- data4[which(data4$freethresh.wildemmer  < 2 ),,drop=F]
data5 <- read.table("sub.cultivar-strangulata.D.sites.pi", header = T, stringsAsFactors = F)
sub5 <- data5[which(data5$cultivar.strangulata  < 2 ),,drop=F]
data6 <- read.table("sub.landrace-strangulata.D.sites.pi", header = T, stringsAsFactors = F)
sub6 <- data6[which(data6$landrace.strangulata  < 2),,drop=F]
sub1$type <- "cultivar/wildemmer(A)"
sub2$type <- "landrace/wildemmer(A)"
sub3$type <- "domemmer/wildemmer(A)"
sub4$type <- "freethresh/wildemmer(A)"
sub5$type <- "cultivar/strangulata(D)"
sub6$type <- "landrace/strangulata(D)"
colnames(sub1) <- c("group","type")
colnames(sub2) <- c("group","type")
colnames(sub3) <- c("group","type")
colnames(sub4) <- c("group","type")
colnames(sub5) <- c("group","type")
colnames(sub6) <- c("group","type")
sub_all <- rbind(sub1,sub2,sub3,sub4)
sub_all <- rbind(sub5,sub6)
sub_all <- rbind(sub1,sub2)
ggplot(sub_all, aes(x = group, group=type, color=type))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
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
  ) +
  theme(legend.text = element_text(size=12),legend.title=element_blank(), legend.position = c(.7, .8))

#Tajima'D
library(ggplot2)
library(RColorBrewer)
color <- brewer.pal(12, "Paired")[c(1,2,8,7)]
color <- brewer.pal(12, "Paired")[c(2,8)]
color <- brewer.pal(8, "Paired")[c(3,1,7)]
color <- brewer.pal(12, "Paired")[c(2,3,1,8,7)]
color <- brewer.pal(12, "Paired")[c(2,8,3)]
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/TajimaD")
data1 <- read.table("landrace_D_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub1 <- data1[which(data1$TajimaD != "NaN" ),,drop=F]
data2 <- read.table("cultivar_D_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub2 <- data2[which(data2$TajimaD != "NaN" ),,drop=F]
data3 <- read.table("wildemmer_B_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub3 <- data3[which(data3$TajimaD != "NaN" ),,drop=F]
data4 <- read.table("domemmer_B_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub4 <- data4[which(data4$TajimaD != "NaN" ),,drop=F]
data5 <- read.table("freethresh_B_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub5 <- data5[which(data5$TajimaD != "NaN" ),,drop=F]
data6 <- read.table("strangulata_D_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub6 <- data6[which(data6$TajimaD != "NaN" ),,drop=F]
sub1$type <- "landrace(D)"
sub2$type <- "cultivar(D)"
sub3$type <- "wildemmer(B)"
sub4$type <- "domemmer(B)"
sub5$type <- "freethresh(B)"
sub6$type <- "strangulata(D)"
sub_all <- rbind(sub1,sub2)
sub_all <- rbind(sub1,sub2,sub6)
sub_all <- rbind(sub3,sub4,sub5)
sub_all <- rbind(sub1,sub2,sub3,sub4,sub5)
ggplot(sub_all, aes(x = TajimaD, group=type, color=type))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
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
  ) +
  theme(legend.text = element_text(size=12),legend.title=element_blank(), legend.position = c(.8, .8))

#计算H12和H2/H1
#路径：/data4/home/yafei/plink_VCF/maf001r209/group/lineage
python ../gene.index.py example_012_pos example_gene.pos example_gene.index 
python ../geno2h.py example_012_geno example_gene.index example.h12.stats example.h12.hap 0  #仅基因
python ../geno2h.py example_012_geno example_gene.index example.h12.stats example.h12.hap 50 #基因及上下游50 SNPs
python ../hap.similarity.py example_012_geno example_gene.index example.h12.sim 100 # similarity between the first two haplotypes

#画图
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/H12/")
filename <- vector()
for(j in c("cultivar","landrace")){
  for(i in c(0,10,20,30,40,50,60,70,80,90,100)){
    A <- paste("hap",i,"/",j,"/chr",c(1,2,7,8,13,14,19,20,25,26,31,32,37,38),"_",j,".A.h12.",i,".stats",sep="")
    B <- paste("hap",i,"/",j,"/chr",c(3,4,9,10,15,16,21,22,27,28,33,34,39,40),"_",j,".B.h12.",i,".stats",sep="")
    D <- paste("hap",i,"/",j,"/chr",c(5,6,11,12,17,18,23,24,29,30,35,36,41,42),"_",j,".D.h12.",i,".stats",sep="")
    filename <- append(filename,c(A,B,D))
  }
}
for(j in c("wildemmer","domemmer","freethresh")){
  for(i in c(0,10,20,30,40,50,60,70,80,90,100)){
    A <- paste("hap",i,"/",j,"/chr",c(1,2,7,8,13,14,19,20,25,26,31,32,37,38),"_",j,".A.h12.",i,".stats",sep="")
    B <- paste("hap",i,"/",j,"/chr",c(3,4,9,10,15,16,21,22,27,28,33,34,39,40),"_",j,".B.h12.",i,".stats",sep="")
    filename <- append(filename,c(A,B))
  }
}
for(j in c("strangulata")){
  for(i in c(0,10,20,30,40,50,60,70,80,90,100)){
    D <- paste("hap",i,"/",j,"/chr",c(5,6,11,12,17,18,23,24,29,30,35,36,41,42),"_",j,".D.h12.",i,".stats",sep="")
    filename <- append(filename,c(D))
  }
}
all <- data.frame(Gene="", H12_1="", H2_1="", genome="", hap="",type="",stringsAsFactors=FALSE)
for(num in c(1:462)){
    out1 <- strsplit(gsub("cultivar.","/", filename[num]) , "/")[[1]][4]
    out2 <- strsplit(out1 , ".h12.")
    data <- read.table(filename[num], header=F,stringsAsFactors = F,sep=" ")
    colnames(data) <- c("Gene", "H12_1", "H2_1")
    data$genome <- out2[[1]][1]
    data$hap <- out2[[1]][2]
    data$type <- "cultivar"
    all <- rbind(all,data)
}
for(num in c(463:924)){
    out1 <- strsplit(gsub("landrace.","/", filename[num]) , "/")[[1]][4]
    out2 <- strsplit(out1 , ".h12.")
    data <- read.table(filename[num], header=F,stringsAsFactors = F,sep=" ")
    colnames(data) <- c("Gene", "H12_1", "H2_1")
    data$genome <- out2[[1]][1]
    data$hap <- out2[[1]][2]
    data$type <- "landrace"
    all <- rbind(all,data)
}
for(num in c(925:1232)){
  out1 <- strsplit(gsub("wildemmer.","/", filename[num]) , "/")[[1]][4]
  out2 <- strsplit(out1 , ".h12.")
  data <- read.table(filename[num], header=F,stringsAsFactors = F,sep=" ")
  colnames(data) <- c("Gene", "H12_1", "H2_1")
  data$genome <- out2[[1]][1]
  data$hap <- out2[[1]][2]
  data$type <- "wildemmer"
  all <- rbind(all,data)
}
for(num in c(1233:1540)){
  out1 <- strsplit(gsub("domemmer.","/", filename[num]) , "/")[[1]][4]
  out2 <- strsplit(out1 , ".h12.")
  data <- read.table(filename[num], header=F,stringsAsFactors = F,sep=" ")
  colnames(data) <- c("Gene", "H12_1", "H2_1")
  data$genome <- out2[[1]][1]
  data$hap <- out2[[1]][2]
  data$type <- "domemmer"
  all <- rbind(all,data)
}
for(num in c(1541:1848)){
  out1 <- strsplit(gsub("freethresh.","/", filename[num]) , "/")[[1]][4]
  out2 <- strsplit(out1 , ".h12.")
  data <- read.table(filename[num], header=F,stringsAsFactors = F,sep=" ")
  colnames(data) <- c("Gene", "H12_1", "H2_1")
  data$genome <- out2[[1]][1]
  data$hap <- out2[[1]][2]
  data$type <- "freethresh"
  all <- rbind(all,data)
}
for(num in c(1849:2002)){
  out1 <- strsplit(gsub("strangulata.","/", filename[num]) , "/")[[1]][4]
  out2 <- strsplit(out1 , ".h12.")
  data <- read.table(filename[num], header=F,stringsAsFactors = F,sep=" ")
  colnames(data) <- c("Gene", "H12_1", "H2_1")
  data$genome <- out2[[1]][1]
  data$hap <- out2[[1]][2]
  data$type <- "strangulata"
  all <- rbind(all,data)
}
write.table(all,"H12.txt",quote = F,sep="\t",row.names = F)

library(ggplot2)
library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(12, "Paired")[c(1,2,8,7)]
color <- brewer.pal(12, "Paired")[c(2,8)]
color <- brewer.pal(8, "Paired")[c(3,1,7)]
color <- brewer.pal(12, "Paired")
color <- brewer.pal(8, "Set3")[c(1,2,3,4,5,7)]
all2 <- read.table("H12.txt",header=T,stringsAsFactors = F)
all2$hap <- factor(all2$hap,levels = c("0.stats","10.stats","20.stats","30.stats","40.stats","50.stats","60.stats","70.stats","80.stats","90.stats","100.stats"))
hap <- c("0.stats","20.stats","40.stats","60.stats","80.stats","100.stats")
all <- all2[which(all2$hap %in% hap ),]
data <- as.data.frame(all2[which(all$type=="landrace"),])
p1 <- ggplot(data, aes(x = as.numeric(H12_1), group=hap, color=hap))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
  theme_classic() + 
  ggtitle("Landrace")+
  xlab("H12")+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  ) +
  theme(title = element_text(size=15),legend.text = element_text(size=15),plot.title = element_text(hjust = 0.5),legend.title=element_blank(), legend.position = c(.8, .6))
data <- as.data.frame(all[which(all$type=="cultivar"),])
p2 <- ggplot(data, aes(x = as.numeric(H12_1), group=hap, color=hap))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
  theme_classic() + 
  ggtitle("Cultivar")+
  xlab("H12")+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  ) +
  theme(title = element_text(size=15),legend.text = element_text(size=15),plot.title = element_text(hjust = 0.5),legend.title=element_blank(), legend.position = c(.8, .6))
data <- as.data.frame(all[which(all$type=="wildemmer"),])
p3 <- ggplot(data, aes(x = as.numeric(H12_1), group=hap, color=hap))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
  theme_classic() + 
  ggtitle("Wild Emmer")+
  xlab("H12")+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  ) +
  theme(title = element_text(size=15),legend.text = element_text(size=15),plot.title = element_text(hjust = 0.5),legend.title=element_blank(), legend.position = c(.8, .6))
data <- as.data.frame(all[which(all$type=="domemmer"),])
p4 <- ggplot(data, aes(x = as.numeric(H12_1), group=hap, color=hap))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
  theme_classic() + 
  ggtitle("Domesticated Emmer")+
  xlab("H12")+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  ) +
  theme(title = element_text(size=15),legend.text = element_text(size=15),plot.title = element_text(hjust = 0.5),legend.title=element_blank(), legend.position = c(.8, .6))
data <- as.data.frame(all[which(all$type=="freethresh"),])
p5 <- ggplot(data, aes(x = as.numeric(H12_1), group=hap, color=hap))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
  theme_classic() + 
  ggtitle("Free-threshing Tetraploid")+
  xlab("H12")+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  ) +
  theme(title = element_text(size=15),legend.text = element_text(size=15),plot.title = element_text(hjust = 0.5),legend.title=element_blank(), legend.position = c(.8, .6))
data <- as.data.frame(all[which(all$type=="strangulata"),])
p6 <- ggplot(data, aes(x = as.numeric(H12_1), group=hap, color=hap))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
  theme_classic() + 
  ggtitle("strangulata")+
  xlab("H12")+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  ) +
  theme(title = element_text(size=15),legend.text = element_text(size=15),plot.title = element_text(hjust = 0.5),legend.title=element_blank(), legend.position = c(.8, .6))

library(gridExtra)
pdf("plot_H12.pdf",height=12,width=18)
grid.arrange(p6,p3,p4,p5,p1,p2,nrow=2)
dev.off()

#H2/H1
library(ggplot2)
all2 <- read.table("H12.txt", header=T, stringsAsFactors = F)
all2$hap <- factor(all2$hap,levels = c("0","10","20","30","40","50","60","70","80","90","100"))
hap <- c("0","10","20","30","40","50","60","70","80","90")
all <- all2[which(all2$hap %in% hap ),]
data <- as.data.frame(all[which(all$type=="cultivar"),])
data$soft <- "NA"
data[which(data$H2_1 <= 0.059), 7] <- 0
data[which(data$H2_1 > 0.059), 7] <- 1
p1 <- ggplot(data = data, aes(x=hap), group=hap, fill=hap) + geom_bar(aes(fill=soft),position="fill")+theme_classic()+ggtitle("Cultivar")+
  theme(legend.position="none",title = element_text(size=20),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),  
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18))
data <- as.data.frame(all[which(all$type=="landrace"),])
data$soft <- "NA"
data[which(data$H2_1 <= 0.059), 7] <- 0
data[which(data$H2_1 > 0.059), 7] <- 1
p2 <- ggplot(data = data, aes(x=hap), group=hap, fill=hap) + geom_bar(aes(fill=soft),position="fill")+theme_classic()+ggtitle("Landrace")+
  theme(legend.position="none",title = element_text(size=20),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18))

data <- as.data.frame(all[which(all$type=="wildemmer"),])
data$soft <- "NA"
data[which(data$H2_1 <= 0.059), 7] <- 0
data[which(data$H2_1 > 0.059), 7] <- 1
p3 <- ggplot(data = data, aes(x=hap), group=hap, fill=hap) + geom_bar(aes(fill=soft),position="fill")+theme_classic()+ggtitle("Wild Emmer")+
  theme(legend.position="none",title = element_text(size=20),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18))
data <- as.data.frame(all[which(all$type=="domemmer"),])
data$soft <- "NA"
data[which(data$H2_1 <= 0.059), 7] <- 0
data[which(data$H2_1 > 0.059), 7] <- 1
p4 <- ggplot(data = data, aes(x=hap), group=hap, fill=hap) + geom_bar(aes(fill=soft),position="fill")+theme_classic()+ggtitle("Domesticated Emmer")+
  theme(legend.position="none",title = element_text(size=20),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18))

data <- as.data.frame(all[which(all$type=="freethresh"),])
data$soft <- "NA"
data[which(data$H2_1 <= 0.059), 7] <- 0
data[which(data$H2_1 > 0.059), 7] <- 1
p5 <- ggplot(data = data, aes(x=hap), group=hap, fill=hap) + geom_bar(aes(fill=soft),position="fill")+theme_classic()+ggtitle("Free-threshing Tetraploid")+
  theme(legend.position="none",title = element_text(size=20),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18))
data <- as.data.frame(all[which(all$type=="strangulata"),])
data$soft <- "NA"
data[which(data$H2_1 <= 0.059), 7] <- 0
data[which(data$H2_1 > 0.059), 7] <- 1
p6 <- ggplot(data = data, aes(x=hap), group=hap, fill=hap) + geom_bar(aes(fill=soft),position="fill")+theme_classic()+ggtitle("Strangulata")+
  theme(legend.position="none",title = element_text(size=20),plot.title = element_text(hjust = 0.5),legend.title=element_blank(),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 18),
        axis.title.x=element_text(size = 18))

library(gridExtra)
pdf("plot_H2_1.pdf",height=11,width=22)
grid.arrange(p6,p3,p4,p5,p2,p1,nrow=2)
dev.off()

#iHH
python ../geno2ihh.py example_012_geno example_gene.index example.gene.hap.ehh 0.05 20

library(ggplot2)
library(RColorBrewer)
color <- brewer.pal(12, "Paired")[c(1,2,8,7)]
color <- brewer.pal(12, "Paired")[c(2,8)]
color <- brewer.pal(8, "Paired")[c(3,1,7)]
color <- brewer.pal(12, "Paired")[c(2,3,1,8,7)]
color <- brewer.pal(12, "Paired")[c(2,8,3)]
#setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/TajimaD")
data1 <- read.table("landrace_D_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub1 <- data1[which(data1$TajimaD != "NaN" ),,drop=F]
data2 <- read.table("cultivar_D_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub2 <- data2[which(data2$TajimaD != "NaN" ),,drop=F]
data3 <- read.table("wildemmer_B_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub3 <- data3[which(data3$TajimaD != "NaN" ),,drop=F]
data4 <- read.table("domemmer_B_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub4 <- data4[which(data4$TajimaD != "NaN" ),,drop=F]
data5 <- read.table("freethresh_B_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub5 <- data5[which(data5$TajimaD != "NaN" ),,drop=F]
data6 <- read.table("strangulata_D_TajimaD.Tajima.D", header=T, stringsAsFactors = F)
sub6 <- data6[which(data6$TajimaD != "NaN" ),,drop=F]
sub1$type <- "landrace(D)"
sub2$type <- "cultivar(D)"
sub3$type <- "wildemmer(B)"
sub4$type <- "domemmer(B)"
sub5$type <- "freethresh(B)"
sub6$type <- "strangulata(D)"
sub_all <- rbind(sub1,sub2)
sub_all <- rbind(sub1,sub2,sub6)
sub_all <- rbind(sub3,sub4,sub5)
sub_all <- rbind(sub1,sub2,sub3,sub4,sub5)
ggplot(sub_all, aes(x = TajimaD, group=type, color=type))+
  #geom_boxplot(fill = '#f8766d', notch = TRUE)+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
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
  ) +
  theme(legend.text = element_text(size=12),legend.title=element_blank(), legend.position = c(.8, .8))




#模拟数据计算pi, tajima's D, 和h12
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/H12/simulation")
data1 <- read.table("landrace.ms.stat", header=F,stringsAsFactors = F)
sub1 <- data1[,c(2,6)]
sub1$type <- "landrace"
data2 <- read.table("wildemmer.ms.stat", header=F,stringsAsFactors = F)
sub2 <- data2[,c(2,6)]
sub2$type <- "wildemmer"
data3 <- read.table("domemmer.ms.stat", header=F,stringsAsFactors = F)
sub3 <- data3[,c(2,6)]
sub3$type <- "domemmer"
data4 <- read.table("freethresh.ms.stat", header=F,stringsAsFactors = F)
sub4 <- data4[,c(2,6)]
sub4$type <- "freethresh"
all <- rbind(sub1,sub2,sub3,sub4)
colnames(all) <- c("Pi","D","type")
all$type <- factor(all$type,levels=c("wildemmer","domemmer","freethresh","landrace"))
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(7, "Set2")
ggplot(all, aes(x = as.numeric(Pi), group=type, color=type))+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
  theme_classic() + 
  ggtitle("Pi Simulation")+
  xlab("Pi")+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  ) +
  theme(title = element_text(size=15),legend.text = element_text(size=15),plot.title = element_text(hjust = 0.5),legend.title=element_blank(), legend.position = c(.8, .8))

data1 <- read.table("landrace.sim.20k.h12", header=F,stringsAsFactors = F)
data1$type <- "landrace"
data2 <- read.table("wildemmer.sim.20k.h12", header=F,stringsAsFactors = F)
data2$type <- "wildemmer"
data3 <- read.table("domemmer.sim.20k.h12", header=F,stringsAsFactors = F)
data3$type <- "domemmer"
data4 <- read.table("freethresh.sim.20k.h12", header=F,stringsAsFactors = F)
data4$type <- "freethresh"
all <- rbind(data1,data2,data3,data4)
colnames(all) <- c("ID","H12","H2/H1","type")
all$type <- factor(all$type,levels=c("wildemmer","domemmer","freethresh","landrace"))
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(7, "Set2")

ggplot(all, aes(x = as.numeric(H12), group=type, color=type))+
  geom_density(size=0.7,geom="area", fill="grey",alpha=0.1) +
  scale_color_manual(values = color)+ 
  theme_classic() + 
  ggtitle("H12 Simulation")+
  xlab("H12")+
  theme(
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=17),
    axis.text.y=element_text(size=20),
    axis.title.y=element_text(size = 20),
    axis.title.x=element_text(size = 20),
  ) +
  theme(title = element_text(size=15),legend.text = element_text(size=15),plot.title = element_text(hjust = 0.5),legend.title=element_blank(), legend.position = c(.3, .85))

#LASSI:画曼哈顿图
#服务器：204:/data2/yafei/004_Vmap3/LASSI/out
#本地：/Users/guoyafei/Documents/02_VmapIII/06_Selection/LASSI/landrace
library(ggplot2)
library(gridExtra)
library(tidyverse)
centro <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv/V4/manhuttan/centro.txt", header=F, stringsAsFactors=F, row.names=1)

setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/LASSI/landrace")
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/LASSI/cultivar")
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/LASSI/freethresh")
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/LASSI/wildemmer")
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/LASSI/domemmer")
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/LASSI/strangulata")

#landrace
names <- paste("landrace_chr",rep(c(1:7),times=3),rep(c("A","B","D"),each=7),".txt",sep="")
title <- paste(rep(c(1:7),times=3),rep(c("A","B","D"),each=7),sep="")
#cultivar
names <- paste("cultivar_chr",rep(c(1:7),times=3),rep(c("A","B","D"),each=7),".txt",sep="")
title <- paste(rep(c(1:7),times=3),rep(c("A","B","D"),each=7),sep="")
#wildemmer
names <- paste("wildemmer_chr",rep(c(1:7),times=2),rep(c("A","B"),each=7),".txt",sep="")
title <- paste(rep(c(1:7),times=2),rep(c("A","B"),each=7),sep="")
#domemmer
names <- paste("domemmer_chr",rep(c(1:7),times=2),rep(c("A","B"),each=7),".txt",sep="")
title <- paste(rep(c(1:7),times=2),rep(c("A","B"),each=7),sep="")
#freethresh
names <- paste("freethresh_chr",rep(c(1:7),times=2),rep(c("A","B"),each=7),".txt",sep="")
title <- paste(rep(c(1:7),times=2),rep(c("A","B"),each=7),sep="")
#strangulata
names <- paste("strangulata_chr",c(1:7),rep("D",each=7),".txt",sep="")
title <- paste(rep(c(1:7),times=1),rep("D",each=7),sep="")

p<-list()
for (i in c(1:7)){
  cen <- centro[,1,drop=F]
  data <- read.table(names[i],header=F,stringsAsFactors=F)
  gwasResults <- data[sample(c(1:nrow(data)), 10000),]
  thresh <- gwasResults[order(gwasResults$V2),][nrow(gwasResults)*0.95,2]
  colnames(gwasResults) <- c("CHR", "P","ID1","ID2","ID3","ID4","BP")
  don <- gwasResults %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) 
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  len <- don %>% group_by(CHR) %>% summarize(center=max(BPcum))
  asx <- c(0,len$center[1:6])
  cen$V2 <- cen$V2+asx
  p[[i]] <- ggplot(don, aes(x=BPcum, y=log(P))) +
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.1) +
    #geom_line(aes(color=as.factor(CHR)))+
    scale_color_manual(values = rep(c("grey","skyblue","grey","skyblue"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
    # Custom the theme:
    #geom_vline(xintercept = 33952048, color = 'skyblue', size = 0.5) + 
    #geom_vline(xintercept = 33956269, color = 'skyblue', size = 0.5) + 
    geom_hline(yintercept = log(thresh), color = 'red', size = 0.5) + 
    geom_point(aes(cen[i,1],0),color="#BF812D")+
    #geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    scale_y_continuous(limits = c(-2,10))+
    theme_classic() +
    ylab("log(T)")+
    ggtitle(title[i])+
    xlab("")+
    theme( 
      plot.title = element_text(color="black", size=15, face="bold"),
      legend.position="none",
      #panel.border = element_blank(),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=11),
      axis.text.y=element_text(size=11),
      axis.title.y=element_text(size = 14),
      axis.title.x=element_text(size = 14),
    )
}
#21chr
pdf("lassi.pdf",height = 20, width = 20)
grid.arrange(p[[1]],p[[8]],p[[15]],p[[2]],p[[9]],p[[16]],p[[3]],p[[10]],p[[17]],p[[4]],p[[11]],p[[18]],p[[5]],p[[12]],p[[19]],p[[6]],p[[13]],p[[20]],p[[7]],p[[14]],p[[21]],nrow=7)
dev.off()
#14chr
pdf("lassi.pdf",height = 20, width = 14)
grid.arrange(p[[1]],p[[8]],p[[2]],p[[9]],p[[3]],p[[10]],p[[4]],p[[11]],p[[5]],p[[12]],p[[6]],p[[13]],p[[7]],p[[14]],nrow=7)
dev.off()
#7chr
pdf("lassi.pdf",height = 20, width = 7)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],nrow=7)
dev.off()

#选择信号鉴定
#测试文件目录
#PI: 204 /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/PI/domemmer-wildemmer.A.sites.pi  domemmer-wildemmer.B.sites.pi
#1	40648	NA
#1	54326	NA
#1	78468	1.49259
#1	380609	NA

#TajimaD: 204 /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/TajimaD/domemmer_A_TajimaD.Tajima.D  domemmer_B_TajimaD.Tajima.D  wildemmer_A_TajimaD.Tajima.D  wildemmer_B_TajimaD.Tajima.D
#CHROM	BIN_START	N_SNPS	TajimaD
#1	0	1	-0.581517
#1	300000	0	nan
#1	400000	5	-0.810035

#iHH: 204 /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/iHH/domemmer.A.gene.hap.ehh  domemmer.B.gene.hap.ehh  wildemmer.A.gene.hap.ehh  wildemmer.B.gene.hap.ehh
#TraesCS3A02G000300 17 4.705882352941177 11.764705882352942 0000000000000000000000000000
#TraesCS3A02G000400 26 17.04615384615385 14.769230769230766 0000000000000
#TraesCS3A02G000500 23 17.075098814229246 40.23715415019763 0000000000000000

#H12: /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/H12/chr1_domemmer.B.h12_40.stats   chr1_domemmer.A.h12_40.stat
#TraesCS2B02G317900 0.07766272189349113 0.4130434782608687
#TraesCS2B02G318000 0.132396449704142 0.25877192982456093
#TraesCS2B02G318100 0.1242603550295858 0.40686274509803827

#SF2: /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/SF2/domemmer.chr01.out wildemmer.chr01.out
#location	LR	alpha
#54326.000000	0.361143	5.024377e+01
#154354.728508	18.037109	3.058327e-06
#254383.457015	18.535522	3.412980e-06

#LASSI: /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/LASSI/wildemmer_chr1.txt window_centers_wildemmer_out_wildemmer_chr1.txt
#1 5.96190000482 2.0 0.0005 0.00393534518886 3 692052.5
#1 0.0 0.0 0.0005 0.00393534518886 3 719188.0
#1 3.24984873273 7.0 0.0005 0.00393534518886 3 793844.0

#692052.5
#719188.0
#793844.0

#写代码，提取受选择位点（区域，基因）
#PI: /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/PI/domemmer-wildemmer.A.005.pi.bed（在两种状态中不一样的位点）
#TajimaD: /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/TajimaD/domemmer_A_005.Tajima.D.bed（各自状态中受选择的区间）
#H12: /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/H12/domemmer.A.40.005.h12.gene（各自状态中不一样的基因）
#SF2: /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/SF2/domemmer.A.005.updown50k.bed（各自状态中不一样的区间）

#找基因
#Gene: /data1/home/yafei/004_Vmap3/selscan/wildemmer-domemmer/Gene
#gff: /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3

bedtools intersect -b ../PI/domemmer-wildemmer.A.005.pi.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > domemmer-wildemmer.A.005.pi.gene
bedtools intersect -b ../PI/domemmer-wildemmer.B.005.pi.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > domemmer-wildemmer.B.005.pi.gene

bedtools intersect -b ../TajimaD/domemmer_A_005.Tajima.D.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > domemmer_A_005.Tajima.D.gene
bedtools intersect -b ../TajimaD/domemmer_B_005.Tajima.D.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > domemmer_B_005.Tajima.D.gene
bedtools intersect -b ../TajimaD/wildemmer_A_005.Tajima.D.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > wildemmer_A_005.Tajima.D.gene
bedtools intersect -b ../TajimaD/wildemmer_B_005.Tajima.D.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > wildemmer_B_005.Tajima.D.gene

cp ../H12/domemmer.A.40.005.h12.gene domemmer.A.40.005.h12.gene
cp ../H12/domemmer.B.40.005.h12.gene domemmer.B.40.005.h12.gene
cp ../H12/wildemmer.A.40.005.h12.gene wildemmer.A.40.005.h12.gene
cp ../H12/wildemmer.B.40.005.h12.gene wildemmer.B.40.005.h12.gene

bedtools intersect -b ../SF2/domemmer.A.005.updown50k.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > domemmer.A.005.updown50k.gene
bedtools intersect -b ../SF2/domemmer.B.005.updown50k.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > domemmer.B.005.updown50k.gene
bedtools intersect -b ../SF2/wildemmer.A.005.updown50k.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > wildemmer.A.005.updown50k.gene
bedtools intersect -b ../SF2/wildemmer.B.005.updown50k.bed -a /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.gff3 -wa  | awk '{print $9}' | awk -F";" '{print $1}' | awk -F"=" '{print $2}' | sort | uniq > wildemmer.B.005.updown50k.gene

#positive control:
Btr1-A (Chr3A(chr13): 65,869,056−65,869,644) 
Btr1-B (Chr3B(chr15): 88,971,298–88,971,838)
Btr1-D (Chr3D(chr17): 55,779,918–55,780,576)
TaTB1-4A: TraesCS4A02G271300()
TaTB1-4B: TraesCS4B02G042700()
TaTB1-4D: TraesCS4D02G040100()

Rht-A1: TraesCS4A02G271000()
Rht-B1: TraesCS4B02G043100()	
Rht-D1: TraesCS4D02G040400()

#画图：UpSetR
library(UpSetR)
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/Gene")
domeA_Name <- c("domemmer.A.005.50k.SF2.gene", "domemmer.A.40.005.h12.gene", "domemmer_A_005.Tajima.D.gene")
domeB_Name <- c("domemmer.B.005.50k.SF2.gene", "domemmer.B.40.005.h12.gene", "domemmer_B_005.Tajima.D.gene")
wildA_Name <- c("wildemmer.A.005.50k.SF2.gene", "wildemmer.A.40.005.h12.gene", "wildemmer_A_005.Tajima.D.gene")
wildB_Name <- c("wildemmer.B.005.50k.SF2.gene", "wildemmer.B.40.005.h12.gene", "wildemmer_B_005.Tajima.D.gene")

p <- list()
for(i in c(1:3)){
  data <- read.table(wildB_Name[i],header=F,stringsAsFactors = F)
  p[[i]] <- data
}

listInput <- list(SF2 = p[[1]][,1], H12 = p[[2]][,1], TajimaD = p[[3]][,1])
upset(fromList(listInput),  order.by = "freq",text.scale = c(2),point.size = 3.5, line.size = 2)

#共表达网络
#204: /data1/home/yafei/004_Vmap3/Network/WGCNA/Leaf4/ModuleGo/module01.txt
