library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Desktop/RAiSD/pdf")
path <- "/Users/guoyafei/Desktop/RAiSD/plot"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors = F)})

#threshold
#library(gdata)
#thresh <- read.xls("thresh.xlsx",sheet=1,row.name=1,na.strings=c("NA","#DIV/0!"))
for (i in seq(1:length(data))){
  filename <- strsplit(names(data)[i], "_")[[1]][1]
  p <- list()
  gwasResults <- data[[i]]
  colnames(gwasResults) <- c("CHR", "BP", "P")
  don <- gwasResults %>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) 
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  p[[i]] <- ggplot(don, aes(x=BPcum, y=P)) +
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
  pdf(paste(filename,".pdf",sep=""),height = 3,width = 6.5)
  print(p[[i]])
  dev.off()
}

#temp selection regions
library(UpSetR)

setwd("/Users/guoyafei/Desktop/RAiSD/")
#data <- read.table("temp.all.plot.txt", header=F, stringsAsFactors = F)

EU <- read.table("temp.EU.gene.plot.txt", header=F,stringsAsFactors = F)
EU <- EU[,1]
#EU <- data[,2]
#EU <- EU[which(EU!="0")]
WA <- read.table("temp.WA.gene.plot.txt", header=F,stringsAsFactors = F)
WA <- WA[,1]
#WA <- data[,2]
#WA <- WA[which(WA!="0")]
CA <- read.table("temp.CA.gene.plot.txt", header=F,stringsAsFactors = F)
CA <- CA[,1]
#CA <- data[,3]
#CA <- CA[which(CA!="0")]
SA <- read.table("temp.SA.gene.plot.txt", header=F,stringsAsFactors = F)
SA <- SA[,1]
#SA <- data[,4]
#SA <- SA[which(SA!="0")]
EA <- read.table("temp.EA.gene.plot.txt", header=F,stringsAsFactors = F)
EA <- EA[,1]
#EA <- data[,5]
#EA <- EA[which(EA!="0")]
AM <- read.table("temp.AM.gene.plot.txt", header=F,stringsAsFactors = F)
AM <- AM[,1]
#AM <- data[,6]
#AM <- EA[which(EA!="0")]


all <- list(
  EU <- EU,
  WA <- WA,
  CA <- CA,
  SA <- SA,
  EA <- EA
)

names(all) <- c("EU","WA","CA","SA","EA")
upset(fromList(all), order.by = "freq")

upset(fromList(all), keep.order = TRUE,text.scale = c(1),point.size = 1.5, line.size = 1)



