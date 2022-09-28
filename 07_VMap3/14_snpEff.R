setwd("/Users/guoyafei/Documents/02_VmapIII/09_snpEff/")
data <- read.table("57.genes.txt",header=T,stringsAsFactors = F)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(8, "Set2")[1:4]

md <- melt(data, id=c("GeneId","TranscriptId"))
md$Type <- "NA"
md[which(md$variable == "variants_impact_HIGH"),5] <- "High"
md[which(md$variable == "variants_impact_LOW"),5] <- "Low"
md[which(md$variable == "variants_impact_MODERATE"),5] <- "Moderate"
md[which(md$variable == "variants_impact_MODIFIER"),5] <- "Modifier"
ggplot(md, aes(x=Type, y=value)) + 
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +xlab("")+ylab("")+
  theme(plot.title = element_text(color="red", size=15, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  

#landrace和cultivar的变异分布
setwd("/Users/guoyafei/Documents/02_VmapIII/09_snpEff")
cultivar <- read.table("cultiavr.out2.txt",header=F,stringsAsFactors = F)
cultivar$type <- "cultivar"
a <- cultivar[which(cultivar$V2 == "LOW"),][sample(1:189673, 3000, replace = FALSE),]
b <- cultivar[which(cultivar$V2 == "HIGH"),][sample(1:6854, 3000, replace = FALSE),]
c <- cultivar[which(cultivar$V2 == "MODERATE"),][sample(1:14333, 3000, replace = FALSE),]
d <- cultivar[which(cultivar$V2 == "MODIFIER"),][sample(1:201019, 3000, replace = FALSE),]
all1 <- rbind(a,b,c,d)

landrace <- read.table("landrace.out2.txt",header=F,stringsAsFactors = F)
landrace$type <- "landrace"
a <- landrace[which(landrace$V2 == "LOW"),][sample(1:189673, 3000, replace = FALSE),]
b <- landrace[which(landrace$V2 == "HIGH"),][sample(1:6854, 3000, replace = FALSE),]
c <- landrace[which(landrace$V2 == "MODERATE"),][sample(1:14333, 3000, replace = FALSE),]
d <- landrace[which(landrace$V2 == "MODIFIER"),][sample(1:201019, 3000, replace = FALSE),]
all2 <- rbind(a,b,c,d)
all <- rbind(all1,all2)

all3 <- rbind(landrace,cultivar)
a <- all3[which(all3$V2 == "LOW"),][sample(1:379346, 3000, replace = FALSE),]
b <- all3[which(all3$V2 == "HIGH"),][sample(1:13708, 3000, replace = FALSE),]
c <- all3[which(all3$V2 == "MODERATE"),][sample(1:28666, 3000, replace = FALSE),]
d <- all3[which(all3$V2 == "MODIFIER"),][sample(1:402038, 3000, replace = FALSE),]
all4 <- rbind(a,b,c,d)
all4$type <- "all"

all5 <- rbind(all,all4)

library(ggplot2)
ggplot(all5, aes(x=V2, y=V3, color=type)) + 
  geom_boxplot() + 
  scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  #scale_color_brewer(brewer.pal(8, "Set2")[3:4])+
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +xlab("")+ylab("")+
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=10),legend.title=element_blank(),axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 10))

#57个seedsize相关基因的变异分布

landrace <- read.table("57.gene.landrace.out2.txt", header=F, stringsAsFactors = F)
cultivar <- read.table("57.gene.cultivar.out2.txt", header=F, stringsAsFactors = F)

landrace$type <- "landrace"
cultivar$type <- "cultivar"

all <- rbind(landrace,cultivar)
ggplot(all, aes(x=V1, y=V3, color=type)) + 
  geom_boxplot() + 
  facet_grid(V2~.)+
  scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
 xlab("")+ylab("")+theme_bw()+
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=15),legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))

#基因功能富集分析
data <- read.table("gene_Go.txt",header=T,sep="\t")
data <- data[1:50,1:9]

ggplot(data=data)+
    geom_bar(aes(x= Term, y=queryitem, fill=(-log(pvalue))), stat='identity') +
    coord_flip() +
    #facet_grid(.~Name,scales="free") +
    #facet_wrap(~Name,ncol = 1)+
    #A
    scale_fill_gradient(expression(-log(pvalue)),low="#FDAE6B", high = "#D94801") +
    #B
    #scale_fill_gradient(expression(p.adjust),low="#9ECAE1", high = "#2171B5") +
    #D
    #scale_fill_gradient(expression(p.adjust),low="#A1D99B", high = "#238B45") +
    ylab("Gene number") +
    xlab("GO term description") +
    #expand_limits(y=c(0,8))+
    #ylim(0,100)+
    theme(
      axis.text.x=element_text(color="black",size=rel(0.3)),
      axis.text.y=element_text(color="black", size=rel(0.3)),
      axis.title.x = element_text(color="black", size=rel(5)),
      axis.title.y = element_blank(),
      legend.text=element_text(color="black",size=rel(0.2)),
      legend.title = element_text(color="black",size=rel(0.7))
      #legend.position=c(0,1),legend.justification=c(-1,0)
      #legend.position="top",
    )+
    #  scale_x_discrete(limits= c("Domesticated_einkorn", "Domesticated_emmer","Durum","EA","EU","Indian_dwarf","Khorasan_wheat","Landrace","Macha","Persian_wheat","Polish_wheat","Rivet_wheat","SCA","Spelt","Tibetan_semi_wild","Urartu","Vavilovii","WA","Wild_Einkorn","Wild_emmer","Xinjiang_wheat"))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))











