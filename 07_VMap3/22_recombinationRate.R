#recombination rate
library(ggplot2)
setwd("/Users/guoyafei/公共资源/01_数据资源/1.0_recombination_rate")
iwgsc1 <- read.table("iwgsc_refseqv1.0_mapping_data.txt", header=T,stringsAsFactors = F)
ggplot(iwgsc1, aes(x=physicalPosition, y=geneticPosition,group=chromosome)) +
  geom_point(size=0.5)+
  facet_wrap(vars(chromosome), nrow = 7,scales = "free")+
  geom_line()+
  theme_bw()

iwgsc2 <- read.table("iwgsc_refseqv1.0_recombination_rate.txt", header=T,stringsAsFactors = F)
ggplot(iwgsc2, aes(x=intervalStart, y=nbOfSnps,group=chromosome)) +
  geom_point(size=0.3)+
  facet_wrap(vars(chromosome), nrow = 7,scales = "free")+
  geom_line()+
  theme_bw()

yafei <- read.table("RecombinationRate_1M_yafei.txt", header=T,stringsAsFactors = F)
ggplot(yafei, aes(x=intervalStart, y=recombinationRate,group=chromosome)) +
  geom_point(size=0.3)+
  facet_wrap(vars(chromosome), nrow = 7,scales = "free")+
  geom_line()+
  theme_bw()

lipeng <- read.table("recombination_rate_geva_lipeng.txt", header=T,stringsAsFactors = F)
ggplot(lipeng, aes(x=Position.bp., y=Map.cM.,group=Chromosome)) +
  geom_point(size=0.3)+
  facet_wrap(vars(Chromosome), nrow = 7,scales = "free")+
  geom_line()+
  theme_bw()

xue <- read.table("slidewindow_recomrate_updown_20M_xpclr_xuebo.txt", header=F,stringsAsFactors = F)
colnames(xue) <- c("chromosome","intervalStart","intervalEnd","nbOfSnps","recombinationRate")
ggplot(xue, aes(x=intervalStart, y=recombinationRate,group=chromosome)) +
  geom_point(size=0.3)+
  facet_wrap(vars(chromosome), nrow = 7,scales = "free")+
  geom_line()+
  theme_bw()

xue <- read.table("RecombinationRate_1M_xuebo.txt", header=F,stringsAsFactors = F)
colnames(xue) <- c("chromosome","intervalStart","intervalEnd","nbOfSnps","recombinationRate")
ggplot(xue, aes(x=intervalStart, y=nbOfSnps,group=chromosome)) +
  geom_point(size=0.3)+
  facet_wrap(vars(chromosome), nrow = 7,scales = "free")+
  geom_line()+
  theme_bw()

ggplot(D, mapping = aes(x = V4)) +
  geom_density( alpha = 0.5, color = "#FF7F00",fill="#FF7F00") +
  theme_bw()
p <- ggplot(don, aes(x=BPcum, y=P)) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  theme_bw() +
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
  geom_hline(yintercept=thresh1,color='#E69F00',linetype = "dashed")
ggplot(all, aes(x=V1, y=V2, color=type)) + 
  geom_bar(stat="identity",width=0.5,position='dodge') + 
  xlab("")+
  ylab("")+
  theme_bw()
  
