#################################################################适应性位点在全基因组上的物理分布################################################################
library(CMplot)
end <- read.table("chrend.txt",header=F,stringsAsFactors = F)
colnames(end) <- c("snp","chr","pos")
setwd("/Users/guoyafei/Desktop/3type/location")
mydata <- read.table("type1.21chr.bed", header=F, stringsAsFactors = F)
mydata <- read.table("type2.21chr.bed", header=F, stringsAsFactors = F)
mydata <- read.table("type3.21chr.bed", header=F, stringsAsFactors = F)
colnames(mydata) <- c("chr","start","pos","snp")
data <- mydata[,c(4,1,3)]
all <- rbind(data,end)
head(data)
# snp         chr       pos
# snp1_1    1        2041
CMplot(all,plot.type="d",bin.size=1e2,col=c("darkgreen","yellow", "red"),file="pdf",memo="snp_density",dpi=300)
#all position
mydata1 <- read.table("type1.21chr.bed", header=F, stringsAsFactors = F)
mydata2 <- read.table("type2.21chr.bed", header=F, stringsAsFactors = F)
mydata3 <- read.table("type3.21chr.bed", header=F, stringsAsFactors = F)
colnames(mydata1) <- c("chr","start","pos","snp")
colnames(mydata2) <- c("chr","start","pos","snp")
colnames(mydata3) <- c("chr","start","pos","snp")
my2 <- mydata2[,c(4,1,3)]
my2$type <- "type2"
end$type <- "type2"
type2 <- rbind(my2,end)
my3 <- mydata3[,c(4,1,3)]
my3$type <- "type3"
end$type <- "type3"
type3 <- rbind(my3,end)
all <- rbind(type1,type2,type3)
#但是我不会把他们画在同一条染色体上


####################################################################曼哈顿图(lfmm & baypass)并标注基因###########################################################
library(qqman)
library(tidyverse)
library(ggrepel)
setwd("/Users/guoyafei/Desktop/3type/曼哈顿")

input <- c("baypass.solar1_21chr.txt","lfmm.solar1_21chr.txt")
output <- c("baypass.solar1_21chr.pdf","lfmm.solar1_21chr.pdf")

point <- read.table("select.30k.21chr.bed", header=F,stringsAsFactors = F)
point$SNP <- paste(point$V1,point$V2,sep="-")
snpsOfInterest <- point$SNP
anno <- point[,c(1,2,7)]
anno$P <- 30
colnames(anno) <- c("CHR","BP","SNP","P")

for (i in 1){
  gwasResults2 <- read.table(input[i], header=T, stringsAsFactors = F)
  gwasResults3 <- gwasResults2[sample(nrow(gwasResults2), 25000), ]
  gwasResults <- rbind(gwasResults3,anno)
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
    mutate( BPcum=BP+tot) %>%
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) 

  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  
  data <- subset(don, is_highlight=="yes")
  all <- merge(data,point,by="SNP")
  p <- ggplot(don, aes(x=BPcum, y=P)) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
    #scale_color_manual(values = rep(c("grey","skyblue","#E69F00"), 7)) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +   
    #Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    theme_bw() +
    geom_vline(xintercept=all$BPcum) +
    scale_color_manual(values = rep(c("grey","#CC79A7","#A07D35"), 7)) +
    #scale_color_manual(values = c("grey","skyblue","#E69F00")) +
    #geom_point(data=all, aes(x=BPcum, y=P,color=as.factor(V6), size=2)) +
    new_scale_color() +
    geom_label_repel(data=all, aes(label=V5,color=as.factor(V6)), size=2,max.overlaps=20) +
    scale_color_manual(values = c("#009E73","#0072B2","#D55E00")) +
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
    xlab("CHR")+
    ylab("Bayes Factor")
  #scale_y_continuous(limits = c(0,7))+
  #geom_point(data=point,aes(x=BPcum,y=-log10(P)),color="red")
  pdf(output[i],height = 2.5,width = 15)
  print(p)
  dev.off()
}

  
  
  
  
  
  
