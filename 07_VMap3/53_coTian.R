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

png("V3.png",height=350, width=1500,res = 300)
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
setwd("/Users/guoyafei/Desktop/coTian/QTL/QTL/box")
file <- c("WBGene00000536", "WBGene00012895", "WBGene00021800", "WBGene00022170")

outfile <- paste(file,".LOD1.pdf",sep="")
#for(i in c(2:length(file))){
for(i in c(2:3)){  
  filename <- paste(file[i],".QTL1.plot.txt",sep="")
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
  pheatmap(data_sort_n, show_rownames=T, show_colnames=T, color = cols,legend_breaks = c(0:2), legend_labels = c("0/0", "0/1", "1/1"),  cluster_col = F, cluster_row = FALSE)
  
  #pheatmap(data_sort_n, show_rownames=T, show_colnames=T, color = cols,legend_breaks = c(-1:2), legend_labels = c("./.", "0/0", "0/1", "1/1"),  cluster_col = F, cluster_row = FALSE)
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
gwasResults2 <- read.table("em.LOD1.out.txt", header=T, stringsAsFactors = F)
#gwasResults2 <- read.table("QTL1.fastGWA", header=T, stringsAsFactors = F)
gwasResults <- gwasResults2[!is.na(gwasResults2$P),]
colnames(gwasResults) <- c("SNP","CHR-o", "BP","P")
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
png("em.png",height=350, width=1500,res = 300)
#png("QTL1.png",height=350, width=1500,res = 300)
p <- ggplot(don, aes(x=BPcum, y=P)) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=0.06) +
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
  scale_y_continuous()+
  #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
  geom_hline(yintercept = 5.298, colour="#000000",linetype=2, size=0.3)+
  geom_hline(yintercept = 6.048, colour="#000000",linetype=2, size=0.3)

  #geom_hline(yintercept = -log10(0.01), colour="red",linetype=2, size=0.6)
print(p)
dev.off()

#QTL2#####
setwd("/Users/guoyafei/Desktop/coTian/QTL/LOD2/")
gwasResults2 <- read.table("ehk.LOD2.out.txt", header=T, stringsAsFactors = F)
#gwasResults2 <- read.table("QTL1.fastGWA", header=T, stringsAsFactors = F)
gwasResults <- gwasResults2[!is.na(gwasResults2$LOD),]
colnames(gwasResults) <- c("SNP","CHR-o", "BP","P")
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

png("ehk.png",height=350, width=1500,res = 300)
p <- ggplot(don, aes(x=BPcum, y=P)) +
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
  geom_hline(yintercept = 5.005, colour="#000000",linetype=2, size=0.3)+
  geom_hline(yintercept = 5.713, colour="#000000",linetype=2, size=0.3)
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
