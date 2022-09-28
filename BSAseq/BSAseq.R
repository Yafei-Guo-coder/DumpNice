#run script following qtl-seq
#install.packages("qqman")
#BiocManager::install("qqman")
library(qqman)
library(tidyverse)
require(gridExtra)
setwd("/Users/hangyuancheng/Desktop/BSA")
snpindex <- read.table("snp_index.txt", header=T, stringsAsFactors = F)
slidewindow <- read.table("sliding_window.txt", header=T, stringsAsFactors = F)
#manhattan
colnames(snpindex) <- c("CHR","POSI","VAR","DEPTH1","DEPTH2","p99","p95","SNP_index1","SNP_index2","delta_SNP_index")
colnames(slidewindow) <- c("CHR","POSI","mean_p99","mean_p95","mean_index1","mean_index2","mean_delta_index")
win <- slidewindow %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(POSI)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset 
  left_join(slidewindow, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, POSI) %>%
  mutate( POSIcum=POSI+tot) 

don <- snpindex %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(POSI)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset 
  left_join(snpindex, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, POSI) %>%
  mutate( POSIcum=POSI+tot) 
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(POSIcum) + min(POSIcum) )/ 2)

#draw snp-index plot
snpindex1 <- ggplot(don, aes(x=POSIcum, y=SNP_index1)) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("darkgreen","grey", "skyblue", "orange","#F4A2A3"),1)) +   #repeat colors number here should match chromosome number
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  coord_cartesian(ylim = c(0,1))+            # reset y range
  geom_line(data = win, aes(x = POSIcum, y = mean_index1) ,colour = "black",size = 1)+
  xlab("Chromosome")+
  ylab("(SNP + InDel) - index")+
  theme_bw() +
  ggtitle("HR Manhattan Plot")+
  theme(
    plot.title = element_text(hjust = 0.5), #title to middle
    legend.position="none",
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
  )
snpindex1

snpindex2 <- ggplot(don, aes(x=POSIcum, y=SNP_index2)) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("darkgreen","grey", "skyblue", "orange","#F4A2A3"),1)) +   #repeat colors number here should match chromosome number
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  coord_cartesian(ylim = c(0,1))+            # reset y range
  geom_line(data = win, aes(x = POSIcum, y = mean_index2) ,colour = "black",size = 1)+
  xlab("Chromosome")+
  ylab("(SNP + InDel) - index")+
  theme_bw() +
  ggtitle("NHR Manhattan Plot")+
  theme(
    plot.title = element_text(hjust = 0.5), #title to middle
    legend.position="none",
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
  )
snpindex2

#draw delta snp-index plot
delta <- ggplot(don, aes(x=POSIcum, y=delta_SNP_index)) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("darkgreen","grey", "skyblue", "orange","#F4A2A3"),1)) +   #repeat colors number here should match chromosome number
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  coord_cartesian(ylim = c(-1,1))+            # reset y range
  geom_line(data = win, aes(x = POSIcum, y = mean_delta_index) ,colour = "black",size = 1)+
  geom_line(data = win, aes(x = POSIcum, y = mean_p95) ,colour = "blue",size = 1)+
  geom_line(data = win, aes(x = POSIcum, y = mean_p99) ,colour = "red",size = 1)+
  geom_line(data = win, aes(x = POSIcum, y = -mean_p95) ,colour = "blue",size = 1)+
  geom_line(data = win, aes(x = POSIcum, y = -mean_p99) ,colour = "red",size = 1)+
  xlab("Chromosome")+
  ylab("Δ(SNP + InDel) - index")+
  theme_bw() +
  ggtitle("Manhattan Plot")+
  theme(
    plot.title = element_text(hjust = 0.5), #title to middle
    legend.position="none",
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
  )
delta

  