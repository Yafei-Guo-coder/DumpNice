#sweed.r
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/sweed")

#sweed----
file <- c("SweeD.cultivar.A","SweeD.cultivar.B","SweeD.cultivar.D","SweeD.domemmer.A","SweeD.domemmer.B","SweeD.freethresh.A","SweeD.freethresh.B","SweeD.landrace.A","SweeD.landrace.B","SweeD.landrace.D","SweeD.strangulata.D","SweeD.wildemmer.A","SweeD.wildemmer.B")
input <- paste(file,".test",sep="")
output <- paste(file,".pdf",sep="")

for (i in c(1:length(input))) {
#p <- list()
data <- read.table(input[i],header=F,stringsAsFactors = F)
colnames(data)[c(1,2,3)] <- c("CHR", "BP","P")
#data <- data[which(data$piR != "Inf"),]
#a <- sample(c(1:dim(data)[1]),10000)
#b <- data[a[order(a)],]
thresh1 <- as.numeric(quantile(data$P,0.99))
#colnames(b)[c(1,2,7)] <- c("CHR", "BP","P")
gwasResults <- data
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
p <- ggplot(don, aes(x=BPcum, y=P)) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  # Add highlighted points
  #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
  # Custom the theme:
  theme_bw() +
  #ylab("log(piR)")+
  ylab("CLR")+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.title.y=element_text(size = 12),
    axis.title.x=element_text(size = 12),
  )+
  #scale_y_continuous(limits = c(0,1))+
  geom_hline(yintercept=thresh1,color='#E69F00',linetype = "dashed")
#geom_hline(yintercept=thresh2,color='#E69F00',linetype = "dashed")
pdf(output[i],width = 7.3,height = 2)
#grid.arrange(p[[1]],p[[2]],nrow=2)
print(p)
dev.off()
}

###############
#画图：UpSetR
library(UpSetR)
#sweed
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/sweed/top5")
A_Name <- c("SweeD.wildemmer.A.gene.uniq","SweeD.domemmer.A.gene.uniq","SweeD.freethresh.A.gene.uniq","SweeD.landrace.A.gene.uniq","SweeD.cultivar.A.gene.uniq")
B_Name <- c("SweeD.wildemmer.B.gene.uniq","SweeD.domemmer.B.gene.uniq","SweeD.freethresh.B.gene.uniq","SweeD.landrace.B.gene.uniq","SweeD.cultivar.B.gene.uniq")
D_Name <- c("SweeD.strangulata.D.gene.uniq","SweeD.landrace.D.gene.uniq","SweeD.cultivar.D.gene.uniq")
p <- list()
for(i in c(1:length(D_Name))){
  data <- read.table(D_Name[i],header=F,stringsAsFactors = F)
  p[[i]] <- data
}

listInput <- list(domemmer = p[[2]][,1], freethresh = p[[3]][,1],landrace = p[[4]][,1], cultivar=p[[5]][,1])
listInput <- list(strangulata = p[[1]][,1], landrace = p[[2]][,1], cultivar=p[[3]][,1])
upset(fromList(listInput), order.by = "freq",text.scale = c(2),point.size = 3.5, line.size = 2,keep.order=T)
 
#piR-fst
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/fst-piR/top5")
A_Name <- c("domemmer-wildemmer.A.gene.uniq","freethresh-domemmer.A.gene.uniq","landrace-freethresh.A.gene.uniq","cultivar-landrace.A.gene.uniq")
B_Name <- c("domemmer-wildemmer.B.gene.uniq","freethresh-domemmer.B.gene.uniq","landrace-freethresh.B.gene.uniq","cultivar-landrace.B.gene.uniq")
D_Name <- c("landrace-strangulata.D.gene.uniq","cultivar-landrace.D.gene.uniq")
p <- list()
for(i in c(1:length(D_Name))){
  data <- read.table(D_Name[i],header=F,stringsAsFactors = F)
  p[[i]] <- data
}

#listInput <- list(domemmer = p[[1]][,1], freethresh = p[[2]][,1],landrace = p[[3]][,1], cultivar=p[[4]][,1])
listInput <- list(landrace = p[[1]][,1], cultivar=p[[2]][,1])
upset(fromList(listInput), order.by = "freq",text.scale = c(2),point.size = 3.5, line.size = 2,keep.order=T)

#inter
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/inter")
name <- c("domemmer.A.gene","domemmer.B.gene","freethresh.A.gene","freethresh.B.gene","landrace.A.gene","landrace.B.gene","landrace.D.gene","cultivar.A.gene","cultivar.B.gene","cultivar.D.gene")
p <- list()
for(i in c(1:length(name))){
  data <- read.table(name[i],header=F,stringsAsFactors = F)
  p[[i]] <- data
}

#listInput <- list(domemmer = p[[1]][,1], freethresh = p[[2]][,1],landrace = p[[3]][,1], cultivar=p[[4]][,1])
listInput <- list(domemmerA =p[[1]][,1], domemmerB = p[[2]][,1], freethreshA = p[[3]][,1],freethreshB = p[[4]][,1],landraceA = p[[5]][,1], landraceB = p[[6]][,1], landracD = p[[7]][,1], cultivarA=p[[8]][,1],cultivarB=p[[9]][,1],cultivarD=p[[10]][,1])
upset(fromList(listInput), order.by = "freq",text.scale = c(2),point.size = 3.5, line.size = 2,keep.order=T)



