#GWAS#####
setwd("/Users/guoyafei/Desktop/coTian/线虫GWAS/V3")
library(qqman)
library(tidyverse)
gwasResults <- read.table("file.V2.txt", header=T, stringsAsFactors = F)
gwasResults2 <- read.table("output.fastGWA", header=T, stringsAsFactors = F)
gwasResults <- gwasResults2[!is.na(gwasResults2$P),]
colnames(gwasResults)[3] <- "BP"
colnames(gwasResults) <- c("CHR", "BP","SNP", "P","data")
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
pdf("test2.pdf",height = 4.5,width = 9)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
      # Show all points
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 7)) +
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) ) +     
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size = 15),
        axis.title.x=element_text(size = 15),
      )+
      scale_y_continuous(limits = c(0,7))+
    #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
    geom_hline(yintercept = -log10(0.05), colour="orange",linetype=2, size=0.6)+
    geom_hline(yintercept = -log10(0.01), colour="red",linetype=2, size=0.6)
dev.off()

png("V3.png",height=350, width=1500, res = 300)
p <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=0.06) +
  scale_color_manual(values = rep(c("#E6767B", "#BA5481"), 7)) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     
  theme_classic() +
  ylab("LOD score")+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=8),
    axis.text.y=element_text(size=8),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size = 8),
  )+
  scale_y_continuous(limits = c(0,7))+
  #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
  geom_hline(yintercept = 4, colour="#000000",linetype=2, size=0.3)
  #geom_hline(yintercept = , colour="#000000",linetype=2, size=0.3)
print(p)
dev.off()
#qqplot
qq_dat <- data.frame(chr=gwasResults[order(gwasResults$P,decreasing=FALSE),1],
                     BP=gwasResults[order(gwasResults$P,decreasing=FALSE),2],
                     SNP=gwasResults[order(gwasResults$P,decreasing=FALSE),3],
                     obs=-log10(sort(gwasResults$P,decreasing=FALSE)),
                     exp=-log10(ppoints(length(gwasResults$P))))

a <- qq_dat[which(qq_dat$exp > 4 & qq_dat$obs > qq_dat$exp),]
write.table(a,"qq.point.txt", sep="\t",quote=F,row.names = F)
pd_qq <- ggplot(data=qq_dat,aes(exp,obs))+
  geom_point(alpha=0.7,color="#7F7F7FFF")+
  geom_abline(color="#D62728FF")+
  xlab("Expected -log10(P-value)")+
  ylab("Observed -log10(P-value)")+
  scale_x_continuous(limits = c(0,7))+
  scale_y_continuous(limits = c(0,7))+
  theme(
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    #axis.line = element_line(size=0.8,color="black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())

png("height_qq.png")
pd_qq
dev.off()

##具体的一些基因的snp的情况#####
setwd("/Users/guoyafei/Desktop/coTian/基因单倍型分析")
#画单倍型图V1######
#GWAS
#name <- read.table("/Users/guoyafei/Desktop/coTian/life.sort.txt", header=F, stringsAsFactors = F)
#QTL
name <- read.table("/Users/guoyafei/Desktop/coTian/QTL/QTL/QTL1.sort.txt", header=F, stringsAsFactors = F)
name <- read.table("/Users/guoyafei/Desktop/coTian/QTL/QTL/QTL2.sort.txt", header=F, stringsAsFactors = F)
#V1
#file <- c("WBGene00000536","WBGene00001500","WBGene00001501","WBGene00001520","WBGene00001752","WBGene00001758","WBGene00001792","WBGene00002162","WBGene00004804","WBGene00012553","WBGene00012895","WBGene00022170")
#file <- c("WBGene00001792","WBGene00004804","WBGene00012895","WBGene00022170")
#file <- c("WBGene00022170")
#V2
#setwd("/Users/guoyafei/Desktop/coTian/newThree")
#file <- c("WBGene00010822","WBGene00010825","WBGene00012961")

#QTL
library(qqman)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
#setwd("/Users/guoyafei/Desktop/coTian/QTL/QTL/box")

#WBGene00001793####
setwd("/Users/guoyafei/Desktop/coTian/filter3/WBGene00001793")
data <- read.table("gsy-1.QTL1.plot.txt", header=T,stringsAsFactors = F,check.names = F)
data <- read.table("gsy-1.QTL2.plot.txt", header=T,stringsAsFactors = F,check.names = F)

data <- read.table("gsy-1.exon.QTL1.plot.txt", header=T,stringsAsFactors = F,check.names = F)
data <- read.table("gsy-1.exon.QTL2.plot.txt", header=T,stringsAsFactors = F,check.names = F)

data <- read.table("gsy-1.intron.QTL1.plot.txt", header=T,stringsAsFactors = F,check.names = F)
data <- read.table("gsy-1.intron.QTL2.plot.txt", header=T,stringsAsFactors = F,check.names = F)
#file <- c("WBGene00000536", "WBGene00012895", "WBGene00021800", "WBGene00022170")
file <- c("gsy-1.QTL1.plot", "gsy-1.QTL2.plot")

outfile <- paste(file,".pdf",sep="")
#for(i in c(2:length(file))){
for(i in c(2:3)){  
  #filename <- paste(file[i],".QTL1.plot.txt",sep="")
  filename <- paste(file[i],".txt",sep="")
  data <- read.table(filename, header=T, check.names = F,stringsAsFactors = F)
  row.names(data) <- data$`CHROM-POS`
  data_sort_r <- as.data.frame(t(data[name$V1,c(2:(dim(data)[2]-1)),drop=F]), stringsAsFactors=FALSE)
  data_sort_n <- as.data.frame(lapply(data_sort_r, as.numeric))
  row.names(data_sort_n) <- row.names(data_sort_r)
  #cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
  cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C")
  #pdf(outfile[i])
  #pdf(outfile[i], height=4,width=6)
  pdf(outfile[i], height=30,width=6)
  #pdf(outfile[i], height=70,width=6)
  pheatmap(data_sort_n, show_rownames=T, show_colnames=T, color = cols,legend_breaks = c(0:2), legend_labels = c("0/0", "0/1", "1/1"),  cluster_col = T, cluster_row = FALSE)

  #pheatmap(data_sort_n, show_rownames=T, show_colnames=T, color = cols,legend_breaks = c(-1:2), legend_labels = c("./.", "0/0", "0/1", "1/1"),  cluster_col = F, cluster_row = FALSE)
  dev.off()
}
hap1 <- c("N2A","A30","A12","A32","A7","A3","A35")
hap2 <- c("A29","A20","A17","A8","A15","A16","A25","A34","A28","A13","A14","A11","A6","A1","A18","A24","A4","A5","A31","A23","A9","A19")
hap3 <- c("A21","A2","A27","A36","A10","A22","KR314A","A33","A26")

sub <- as.data.frame(data[-c(1:2),c(1,19),drop=F])
sub$hap <- NA
sub[which(sub$`CHROM-POS` %in% hap1),3] <- "hap1"
sub[which(sub$`CHROM-POS` %in% hap2),3] <- "hap2"
sub[which(sub$`CHROM-POS` %in% hap3),3] <- "hap3"
sub$`line-40` <- as.numeric(sub$`line-40`)
ggplot(sub, aes(y = `line-40`, group=as.factor(sub$hap), fill=as.factor(sub$hap)))+
  geom_boxplot(notch = F) +
  theme_classic() + 
  #scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
  scale_fill_manual(values = c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C"))+
  ylab("Lifespan")+
  theme(legend.title="")+
  #ggtitle(colnames(out)[m])+
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
#otherGene####
setwd("/Users/guoyafei/Desktop/coTian/filter3/otherGene/")
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C")
file <- c("WBGene00001793","WBGene00000870", "WBGene00003144", "WBGene00009245", "WBGene00013261", "WBGene00013591")
fileInput <- c(paste(c(paste(file,"exon",sep=".")),"QTL1",sep="."))
fileOutput <- paste(fileInput, ".heatmap.pdf", sep="")
for(i in c(1:length(fileInput))){ 
  filename <- paste(fileInput[i],".plot.txt",sep="")
  data <- read.table(filename, header=T, check.names = F, stringsAsFactors = F)
  row.names(data) <- data$`CHROM-POS`
  data_sort_r <- as.data.frame(t(data[name$V1,c(2:(dim(data)[2]-1)),drop=F]), stringsAsFactors=FALSE)
  data_sort_n <- as.data.frame(lapply(data_sort_r, as.numeric))
  row.names(data_sort_n) <- row.names(data_sort_r)
  cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C")
  pdf(fileOutput[i], height=3,width=5)
  pheatmap(data_sort_n, show_rownames=T, show_colnames=T, color = cols,legend_breaks = c(0:2), legend_labels = c("0/0", "0/1", "1/1"),  cluster_col = T, cluster_row = FALSE)
  dev.off()
}
#exon
#WBGene00001793
hap1 <- c("N2A","A30","A12","A32","A7","A3","A35")
hap2 <- c("A29","A20","A17","A8","A15","A16","A25","A34","A28","A13","A14","A11","A6","A1","A18","A24","A4","A5","A31","A23","A9","A19")
hap3 <- c("A21","A2","A27","A36","A10","A22","KR314A","A33","A26")
WBGene00001793 <- list(hap1,hap2,hap3)
#WBGene00000870
hap1 <- c("N2A","A30", "A12", "A32","A7","A3","A35")
hap2 <- c("A17","A8","A20", "A29", "A13", "A19","A6","A9","A11","A1","A14", "A15", "A16", "A24", "A25","A34","A4","A5","A23", "A28", "A31", "A18")
hap3 <- c("A21","A2","A27", "A36", "A10", "A22", "KR314A", "A33", "A26")
WBGene00000870 <- list(hap1,hap2,hap3)
#WBGene00003144
hap1 <- c("N2A", "A30", "A12", "A32","A7","A3","A35")
hap2 <- c("A17","A8","A20","A29", "A13", "A19","A6","A9","A11","A1","A14", "A15", "A16", "A18", "A24", "A25", "A34","A4","A5","A23", "A28", "A31")
hap3 <- c("A21","A2","A27","A36", "A10", "A22", "KR314A", "A33","A26")
WBGene00003144 <- list(hap1,hap2,hap3)
#WBGene00009245
hap1 <- c("N2A","A30","A12","A7","A3","A35")
hap2 <- c("A17","A8","A20","A29","A13","A19","A6","A9","A11","A16","A24","A25","A34","A4","A5","A23","A28","A31","A14","A18","A32")
hap3 <- c("A1","A2","A15","A21","A27","A10","A22","KR314A","A33","A26","A36")
WBGene00009245 <- list(hap1,hap2,hap3)
#WBGene00013261
hap1 <- c("N2A","A30","A12","A32","A7","A3","A35","A31")
hap2 <- c("A29","A8","A5","A13","A24","A9","A15","A19","A20","A17","A11","A25","A23","A14","A18","A16","A34","A6","A4","A1")
hap3 <- c("A21","A2","A36","A10","A22","A28","KR314A","A33","A26","A27")
WBGene00013261 <- list(hap1,hap2,hap3)
#WBGene00013591
hap1 <- c("N2A","A30","A12","A32","A7","A3","A35","A18")
hap2 <- c("A17","A8","A20","A29","A13","A6","A9","A11","A16","A24","A25","A4","A5","A23","A28","A31","A19","A34","A14")
hap3 <- c("A21","A2","A27","A36","A10","A22","KR314A","A33","A26","A1","A15")
WBGene00013591 <- list(hap1,hap2,hap3)
gene <- list(WBGene00001793, WBGene00000870, WBGene00003144, WBGene00009245,WBGene00013261, WBGene00013591,WBGene00001793,WBGene00000870, WBGene00003144, WBGene00009245,WBGene00013261, WBGene00013591)

file <- c("WBGene00001793","WBGene00000870", "WBGene00003144", "WBGene00009245", "WBGene00013261", "WBGene00013591")
fileInput <- c(c(paste(file,"exon.QTL1",sep="."), paste(file,"exon.QTL2",sep=".")))
fileOutput <- paste(fileInput, ".pdf",sep="")

for(i in c(1:length(fileInput))){ 
  filename <- paste(fileInput[i],".plot.txt",sep="")
  data <- read.table(filename, header=T, check.names = F, stringsAsFactors = F)
  sub <- as.data.frame(data[-c(1:2),c(1,dim(data)[2]),drop=F])
  sub$hap <- NA
  sub[which(sub$`CHROM-POS` %in% gene[[i]][[1]]),3] <- "hap1"
  sub[which(sub$`CHROM-POS` %in% gene[[i]][[2]]),3] <- "hap2"
  sub[which(sub$`CHROM-POS` %in% gene[[i]][[3]]),3] <- "hap3"
  sub$`line-40` <- as.numeric(sub$`line-40`)
  p <- ggplot(sub, aes(y = `line-40`, group=as.factor(sub$hap), fill=as.factor(sub$hap)))+
    geom_boxplot(notch = F) +
    theme_classic() + 
    #scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
    scale_fill_manual(values = c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C"))+
    ylab("Lifespan")+
    theme(legend.title="")+
    #ggtitle(colnames(out)[m])+
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
  
  pdf(fileOutput[i], height=4,width=5)
  print(p)
  dev.off()
}


#统计每一个位点和表型的相关性以及画箱线图进行展示##########
library(ggplot2)
#GWAS
#name <- read.table("/Users/guoyafei/Desktop/coTian/life.sort.txt", header=F, stringsAsFactors = F)
#QTL
name <- read.table("/Users/guoyafei/Desktop/coTian/QTL/QTL/QTL1.sort.txt", header=F, stringsAsFactors = F)
name <- read.table("/Users/guoyafei/Desktop/coTian/QTL/QTL/QTL2.sort.txt", header=F, stringsAsFactors = F)
#V1
#file <- c("WBGene00000536","WBGene00001500","WBGene00001501","WBGene00001520","WBGene00001752","WBGene00001758","WBGene00001792","WBGene00002162","WBGene00004804","WBGene00012553","WBGene00012895","WBGene00022170")
#V2
#setwd("/Users/guoyafei/Desktop/coTian/newThree")
#file <- c("WBGene00010822","WBGene00010825","WBGene00012961")
#QTL
setwd("/Users/guoyafei/Desktop/coTian/QTL/QTL/box")
file <- c("WBGene00000536", "WBGene00012895", "WBGene00021800", "WBGene00022170")

#outfile <- paste(file,".F2.txt",sep="")
outfileBox <- paste(file,".LOD2.box.pdf",sep="")
for(i in c(1:length(file))){
  mR <- vector()
  adR <- vector()
  P <- vector()
  filename <- paste(file[i],".QTL2.plot.txt",sep="")
  data <- read.table(filename, header=T, check.names = F,stringsAsFactors = F)
  row.names(data) <- data$`CHROM-POS`
  data_sort_r <- as.data.frame(t(data[name$V1,c(2:(dim(data)[2]))]), stringsAsFactors=FALSE)
  data_sort_n <- as.data.frame(lapply(data_sort_r, as.numeric))
  rownames(data_sort_n) <- rownames(data_sort_r)
  for(j in c(1:(dim(data_sort_n)[1]-1))){
    a <- summary(lm(as.numeric(data_sort_n[j,])~as.numeric(data_sort_n[dim(data_sort_n)[1],])))
    mR <- append(mR,as.numeric(a$r.squared))
    adR <- append(adR,as.numeric(a$adj.r.squared))
    P <- append(P, a$coefficients[2,4])
  }
  #all <- as.data.frame(cbind(row.names(data_sort_r)[1:(dim(data_sort_r)[1]-1)],mR,adR,P),stringsAsFactors = F)
  #colnames(all) <- c("SNP","R-squared", "Adjusted-R-squared", "F-statistic-p-value")
  #write.table(all,outfile[i], quote=F, row.names = F,sep="\t")
  all <- as.data.frame(cbind(row.names(data_sort_r)[1:(dim(data_sort_r)[1]-1)],mR,P),stringsAsFactors = F)
  colnames(all) <- c("SNP","R-squared", "F-statistic-p-value")
  #write.table(all,outfile[i], quote=F, row.names = F,sep="\t")
  #sub <- all[which(all$`F-statistic-p-value` < 0.05),]
  sub <- all
  if (dim(sub)[1] == 0) { 
    next
  }
  out <- as.data.frame(t(data_sort_n[c(sub$SNP,"line-40"),]))
  pdf(outfileBox[i])
  for(m in c(1:(dim(out)[2]-1))){
    figure <- ggplot(out, aes(y = out$`line-40`, group=as.factor(out[,m]), fill=as.factor(out[,m])))+
      geom_boxplot( notch = F) +
      theme_classic() + 
      #scale_fill_manual(values = c("#FDC086","#BEAED4","#ffffb3","#fccde5")) +
      scale_fill_manual(values = c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C"))+
      ylab("Lifespan")+
      theme(legend.title="")+
     ggtitle(colnames(out)[m])+
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
    print(figure)
  }
  
  dev.off()
}

#QTL#########
library(qqman)
library(tidyverse)
#QTL1
setwd("/Users/guoyafei/Desktop/coTian/QTL/LOD1/")

setwd("/Users/guoyafei/Desktop/coTian/filter3")
gwasResults2 <- read.table("hk.LOD1.out.txt", header=T, stringsAsFactors = F)
gwasResults2 <- read.table("out1.hk.LOD1_filter2.txt", header=T, stringsAsFactors = F)
gwasResults2 <- read.table("em.LOD1.out.txt", header=T, stringsAsFactors = F)
gwasResults2 <- read.table("/Users/guoyafei/Desktop/coTian/filter2/hk.LOD1_out.txt", header=T, stringsAsFactors = F)
#gwasResults2 <- read.table("QTL1.fastGWA", header=T, stringsAsFactors = F)
gwasResults <- gwasResults2[!is.na(gwasResults2$LOD),]
colnames(gwasResults) <- c("SNP","CHR-o", "BP","P")

gwasResults <- gwasResults2[!is.na(gwasResults2$lod),]
colnames(gwasResults) <- c("CHR-o", "BP","P","SNP")
#colnames(gwasResults) <- c("CHR-o","SNP", "BP","A1","A2","N","AF1","BETA","SE","P")
gwasResults$CHR <- NA
gwasResults[which(gwasResults$`CHR-o` == "1"),5] <- "I"
gwasResults[which(gwasResults$`CHR-o` == "2"),5] <- "II"
gwasResults[which(gwasResults$`CHR-o` == "3"),5] <- "III"
gwasResults[which(gwasResults$`CHR-o` == "4"),5] <- "IV"
gwasResults[which(gwasResults$`CHR-o` == "5"),5] <- "V"
gwasResults[which(gwasResults$`CHR-o` == "6"),5] <- "X"
gwasResults[which(gwasResults$`CHR-o` == "7"),5] <- "MtDNA"
gwasResults$CHR <- factor(gwasResults$CHR,levels = c("I","II","III","IV","V","X","MtDNA"))
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
png("/Users/guoyafei/Desktop/coTian/filter3/hk_lod1.png",height=350, width=1500,res = 300)
#png("QTL1.png",height=350, width=1500,res = 300)
p <- ggplot(don, aes(x=BPcum, y=P)) +
  # Show all points
  geom_line( aes(color=as.factor(CHR)), alpha=1, size=0.6) +
  scale_color_manual(values = rep(c("#10557B", "#69BAB3"), 7)) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     
  theme_classic() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=8),
    axis.text.y=element_text(size=8),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size = 8),
  )+
  ylab("LOD score")+
  scale_y_continuous()+
  #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
  #geom_hline(yintercept = 5.298, colour="#000000",linetype=2, size=0.3)+
  #geom_hline(yintercept = 6.048, colour="#000000",linetype=2, size=0.3)
  #geom_hline(yintercept = 5.12, colour="#000000",linetype=2, size=0.3)+
  #geom_hline(yintercept = 5.81, colour="#000000",linetype=2, size=0.3)
  geom_hline(yintercept = 2.68, colour="#000000",linetype=2, size=0.3)+
  geom_hline(yintercept = 3.5, colour="#000000",linetype=2, size=0.3)

  #geom_hline(yintercept = -log10(0.01), colour="red",linetype=2, size=0.6)
print(p)
dev.off()

#QTL2#####
setwd("/Users/guoyafei/Desktop/coTian/QTL/LOD2/")
gwasResults2 <- read.table("ehk.LOD2.out.txt", header=T, stringsAsFactors = F)

gwasResults2 <- read.table("hk.LOD2.out.txt", header=T, stringsAsFactors = F)
gwasResults2 <- read.table("out1.hk.LOD2_filter2.txt", header=T, stringsAsFactors = F)

gwasResults2 <- read.table("/Users/guoyafei/Desktop/coTian/filter2/hk.LOD2_out.txt", header=T, stringsAsFactors = F)
gwasResults2 <- read.table("/Users/guoyafei/Desktop/coTian/filter2/out1.hk.LOD2_filter2.txt", header=T, stringsAsFactors = F)

#gwasResults2 <- read.table("QTL1.fastGWA", header=T, stringsAsFactors = F)
gwasResults <- gwasResults2[!is.na(gwasResults2$LOD),]
colnames(gwasResults) <- c("SNP","CHR-o", "BP","P")

gwasResults <- gwasResults2[!is.na(gwasResults2$lod),]
colnames(gwasResults) <- c("CHR-o", "BP","P","SNP")

gwasResults <- gwasResults2[!is.na(gwasResults2$lod),]
colnames(gwasResults) <- c("CHR-o", "BP","P","SNP")
#colnames(gwasResults) <- c("CHR-o","SNP", "BP","A1","A2","N","AF1","BETA","SE","P")
gwasResults$CHR <- NA
gwasResults[which(gwasResults$`CHR-o` == "1"),5] <- "I"
gwasResults[which(gwasResults$`CHR-o` == "2"),5] <- "II"
gwasResults[which(gwasResults$`CHR-o` == "3"),5] <- "III"
gwasResults[which(gwasResults$`CHR-o` == "4"),5] <- "IV"
gwasResults[which(gwasResults$`CHR-o` == "5"),5] <- "V"
gwasResults[which(gwasResults$`CHR-o` == "6"),5] <- "X"
gwasResults[which(gwasResults$`CHR-o` == "7"),5] <- "MtDNA"
gwasResults$CHR <- factor(gwasResults$CHR,levels = c("I","II","III","IV","V","X","MtDNA"))

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

png("hk_lod2_filter.png",height=350, width=1500,res = 300)
p <- ggplot(don, aes(x=BPcum, y=P)) +
  # Show all points
  geom_line( aes(color=as.factor(CHR)), alpha=1, size=0.6) +
  scale_color_manual(values = rep(c("#E6767B", "#BA5481"), 7)) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     
  theme_classic() +
  ylab("LOD score")+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=8),
    axis.text.y=element_text(size=8),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size = 8),
  )+
  #scale_y_continuous(limits = c(0,7))+
  #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
  #geom_hline(yintercept = 5.005, colour="#000000",linetype=2, size=0.3)+
  #geom_hline(yintercept = 5.713, colour="#000000",linetype=2, size=0.3)
  geom_hline(yintercept = 2.63, colour="#000000",linetype=2, size=0.3)+
  geom_hline(yintercept = 3.47, colour="#000000",linetype=2, size=0.3)
print(p)
dev.off()

#qqplot####
qq_dat <- data.frame(obs=-log10(sort(gwasResults$P,decreasing=FALSE)),
                     exp=-log10(ppoints(length(gwasResults$P))))

#a <- qq_dat[which(qq_dat$exp > 4 & qq_dat$obs > qq_dat$exp),]
#write.table(a,"qq.point.txt", sep="\t",quote=F,row.names = F)
png("QTL1-1.png",height = 3,width = 3)
ggplot(data=qq_dat,aes(exp,obs))+
  geom_point(alpha=0.7,color="#7F7F7FFF")+
  geom_abline(color="#D62728FF")+
  xlab("Expected -log10(P-value)")+
  ylab("Observed -log10(P-value)")+
  scale_x_continuous(limits = c(0,7))+
  scale_y_continuous(limits = c(0,7))+
  theme(
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    #axis.line = element_line(size=0.8,color="black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())
dev.off()
