#SDS
#66:/data1/home/yafei/004_Vmap3/SDS-master
library(qqman)
library(tidyverse)
require(gridExtra)
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/SDS")
chr1 <- read.table("test.txt", header=T, stringsAsFactors = F)

#manhattan
gwasResults2 <- chr1
thresh <- gwasResults2[order(gwasResults2$P),][nrow(gwasResults)*0.05,4]
#gwasResults2$CHR = 1
colnames(gwasResults2) <- c("CHR","BP","SNP","P")
  gwasResults <- gwasResults2
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
p <- ggplot(don, aes(x=BPcum, y=-log(P))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    geom_hline(yintercept = -log(thresh), color = 'red', size = 0.5) +   
    scale_y_continuous(limits = c(0, 100))+
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
    )

pdf("Landrace.pdf",height = 9,width = 9)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=3)
dev.off()
pdf("Cultivar.pdf",height = 9,width = 9)
grid.arrange(p[[4]],p[[5]],p[[6]],nrow=3)
dev.off()


#qq

qq_dat <- data.frame(obs=-log10(sort(data$P,decreasing=FALSE)),
                         exp=-log10( ppoints(length(data$P))))
ggplot(data=qq_dat,aes(exp,obs))+
      geom_point(alpha=0.7,color="#7F7F7FFF")+
      geom_abline(color="#D62728FF")+
      xlab("Expected -log10(P-value)")+
      ylab("Observed -log10(P-value)")+
      scale_x_continuous(limits = c(0,7))+
      scale_y_continuous(limits = c(0,7))+
      #ggtitle(name)+
      theme(
        plot.title = element_text(color="red", size=20, face="bold.italic"),
        axis.title = element_text(size=12,face="bold"),
        axis.text = element_text(face="bold",size=8,color = "black"),
        #axis.line = element_line(size=0.8,color="black"),
        axis.ticks= element_line(size=0.8,colour = "black"),
        panel.grid =element_blank(),
        panel.border = element_rect(fill=NA,size = 0.8),
        panel.background = element_blank())
pdf(name1,width = 6,height = 3)
grid.arrange(p[[1]],p[[2]],nrow=1)
dev.off()

