library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Desktop/lassip/pdf")
path <- "/Users/guoyafei/Desktop/lassip/plot"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F)})

#threshold
#library(gdata)
#thresh <- read.xls("thresh.xlsx",sheet=1,row.name=1,na.strings=c("NA","#DIV/0!"))
for (i in seq(1:length(data))){
  filename <- strsplit(names(data)[i], "_")[[1]][1]
  p <- list()
  gwasResults <- data[[i]]
  colnames(gwasResults)[1:2] <- c("CHR", "BP")
  don <- gwasResults %>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) 
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  p[[i]] <- ggplot(don, aes(x=BPcum, y=pop1_T)) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "#56B4E9", "#D55E00"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +   
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=8),
      axis.text.y=element_text(size=10),
      axis.title.y=element_text(size = 10),
      axis.title.x=element_text(size = 10),
    )
  #scale_y_continuous(limits = c(0,7))+
  #geom_hline(yintercept = -log10(0.0699), colour="red",linetype=2, size=1)
  pdf(paste(filename,"_T.pdf",sep=""),height = 3,width = 6.5)
  print(p[[i]])
  dev.off()
}

all <- data[[11]]
ggplot(all, aes(pop1_g2g1, pop1_g123)) +
  geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
  #geom_point(data = sub_gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
  geom_vline(xintercept=0.059,color='red',linetype = "dotted")+
  #geom_hline(yintercept=0.059,color='red',linetype = "dotted")+
  #geom_text_repel(data = sub_gene1,aes(V2, V3, label = name),max.overlaps = 100) +
  xlab("H2/H1")+
  ylab("H12")+
  theme_bw()+
  theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))


