#working directory
#yafei@204:/data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225/PPD
#PPD-D1
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/10_Gene")
data<- read.table("all.bin.depth.txt",header=T,row.names = 1,stringsAsFactors = F)
data <- read.table("ppd_331taxa.txt",header=T,stringsAsFactors = F)
sub <- data
rownames(data) <- round(data[,1])
sub <- data[,2:332]
names <- colnames(sub)
for(i in c(1:331)){
    sub[which(sub[,i] >0 & sub[,i] <10),i] <- 10
    sub[which(sub[,i] >10 & sub[,i] <30),i] <- 25
    sub[which(sub[,i] >30),i] <- 50
 }
annotation <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/select_taxa4.txt", header=T,  stringsAsFactors = F, sep="\t")
anno <- annotation[,c(2,9),drop=FALSE]

AB_anno <- anno[which(anno$VMap3 %in% colnames(sub)),]
rownames(AB_anno) <- AB_anno$VMap3
AB_anno <- AB_anno[,2,drop=F]
AB_anno <- AB_anno[!is.na(AB_anno$heatmap),1,drop=F]
labels_row =  rownames(sub)
labels_row[1:150] <- ""
labels_row[361:511] <- ""
labels_row[152:359] <- ""
ann_color = list(
  #Growing_Habit = c(Facultative = "yellow", Spring="orange", Winter="blue"),
  heatmap = c(Strangulata = "#8C510A", EU="#66C2A5", WA= "#FC8D62",CA="#8DA0CB", EA ="#E78AC3",SA="#A6D854",Tibet="#FFD92F", Cultivar="#999999"))
sub <- sub[, rownames(AB_anno)]
labels_row <- rep("",times=32)
labels_row[28] <- "16bp"
labels_row[29] <- "5bp"
labels_row[32] <- "2k"
pheatmap(sub, show_rownames=T, labels_row=labels_row,show_colnames=F, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), annotation = AB_anno,cluster_col = F,annotation_names_col = F,annotation_colors = ann_color, cluster_row = FALSE) 

#Haplotype Map: PPD-D1
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/10_Gene")
data <- read.table("ppd_331taxa.txt", header = T, stringsAsFactors = F)
taxa <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/10_Gene/cluster_70.txt",header=T,row.names = 1, stringsAsFactors = F)
taxa <- taxa[!is.na(taxa$VMap3),]
rownames(taxa) <- taxa[,1]
pdf("ppd.pdf")
for( j in c(1:dim(sub)[1])){
    loc <- taxa[,5:6]
    loc$type <- as.numeric(data[j,rownames(taxa)])
    loc$value <- 1
    wheat_reshape <- cast(loc,Latitude_70+Longitude_70~type) 
    wheat_reshape2 <- as.data.frame(wheat_reshape)
    name <- names(table(loc$type))
    if((j-1)%%4 !=0){
      mapPies(wheat_reshape2,xlim=c(-10,130),ylim=c(10,70),nameX="Longitude_70",nameY="Latitude_70",nameZs=name,symbolSize=1.5,
              zColours=brewer.pal(8, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9],main="tit")
    }else{
      par(mfrow=c(2,2),oma=c(1,1,1,1), mar=c(0,1,0,1), cex=1)
      mapPies(wheat_reshape2,xlim=c(-10,130),ylim=c(10,70),nameX="Longitude_70",nameY="Latitude_70",nameZs=name,symbolSize=1.5,
              zColours=brewer.pal(12, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9], main="tit")
    }
    }
dev.off()
#fst
library(ggplot2)
data <- read.table("topP.txt1",header=T,stringsAsFactors=F)
col3 <- colorRampPalette(c("white","blue"))
data$pop1 <- factor(stbp$V5, levels = c("EU","WA","CA","EA","SA"))
data$pop2 <- factor(stbp$V6, levels = c("EU","WA","CA","EA","SA"))

pdf("gene.pdf",height = 6,width = 7)
ggplot(data = data) + 
  geom_tile(aes(pop1, pop2, fill = gene)) +
  scale_fill_gradientn(colours = col3(100), limits = c(-0.1, 0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("")+
  ylab("")
dev.off()

pdf("16bp.pdf",height = 6,width = 7)
ggplot(data = data) + 
  geom_tile(aes(pop1, pop2, fill = sixbp)) +
  scale_fill_gradientn(colours = col3(100), limits = c(-0.1, 0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("")+
  ylab("")
dev.off()

pdf("5bp.pdf",height = 6,width = 7)
ggplot(data = data) + 
  geom_tile(aes(pop1, pop2, fill = fivebp)) +
  scale_fill_gradientn(colours = col3(100), limits = c(-0.1, 0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("")+
  ylab("")
dev.off()

pdf("2k.pdf",height = 6,width = 7)
ggplot(data = data) + 
  geom_tile(aes(pop1, pop2, fill = twobp)) +
  scale_fill_gradientn(colours = col3(100), limits = c(-0.1, 0.5))+
  theme_classic()+
  
  theme(axis.text.x = element_text(size = 25), axis.title.x = element_text(size = 25),axis.text.y = element_text(size = 25),axis.title.y = element_text(size = 25))+
  xlab("")+
  ylab("")
dev.off()

#大麦中检测
#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
setwd("/Users/guoyafei/Documents/01_Migration/07_ManuScript/NP/reviewer")

#A
annotation_col <- read.table("Qingke-anno.txt",header=T,stringsAsFactors = F,sep="\t")
rownames(annotation_col) = c(1:172)
anno2 <- annotation_col[which(annotation_col$classfication != "NA"),] 
seq <- anno2[,2]
data <- read.table("/Users/guoyafei/Documents/01_Migration/07_ManuScript/NP/reviewer/Qingke.PPD-H1.sample.10k.txt", header=F, stringsAsFactors = F)
colnames(data) <- c(1:172)
#plot haplotype heatmap
anno <- annotation_col[,3,drop=FALSE]

cols <- c("#225EA8","#DEEBF7","#FEB24C","#BD0026")

ann_color = list(
  Region_sub = c(cultivar = "#8C510A", landrace = "#DFC27D", qingkecultivar = "#F6E8C3", qingkelandrace="#66C2A5", Tibetanweedybarley= "#FC8D62",wild="#8DA0CB"))

pdf("A.pdf")
  data <- data[,seq]
  colnames(data) <- c(1:171)
  pheatmap(data, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = anno, annotation_colors = ann_color,annotation_names_col = F)
dev.off()

##reviewer
library(ggplot2)
setwd("/Users/guoyafei/Documents/01_Migration/07_ManuScript/NP/reviewer")
data <- read.table("gene_up100k.result.txt", header=T, stringsAsFactors = F)
sub1 <- data[-1,c(1,2,61)]
sub2 <- as.data.frame(t(data[c(1,3,62,16,20,17,51),1:2]))
data2 <- as.data.frame(t(data[c(1,3,62,16,20,17,51),3:60]))
colnames(data2) <- c("Altitude","pos_33953684","pos_34000377","pos_33963940","pos_33966039","pos_33964177","pos_34000269")
r2<-as.data.frame(t(data[c(1,3,62,16,20,17,51),61]))

library(gridExtra)
a <- ggplot(data2,aes(x=Altitude,y=pos_33953684))+geom_point()+stat_smooth(method=lm)+geom_text(x=1000, y=0.9, label="r2=0.778",size=3.5) +ggtitle("33953684")+ylab("Allele Frequency")+theme_classic()+theme(plot.title = element_text(colour = "red",size = 10,face = "bold"))

b <- ggplot(data2,aes(x=Altitude,y=pos_34000377))+geom_point()+stat_smooth(method=lm)+geom_text(x=3000, y=0.015, label="r2=0.758",size=3.5)  +ggtitle("34000377")+ylab("Allele Frequency")+theme_classic()+theme(plot.title = element_text(colour = "red",size = 10,face = "bold"))
c <- ggplot(data2,aes(x=Altitude,y=pos_33963940))+geom_point()+stat_smooth(method=lm)+geom_text(x=3000, y=0.018, label="r2=0.681",size=3.5)  +ggtitle("33963940")+ylab("Allele Frequency")+theme_classic()+theme(plot.title = element_text(colour = "red",size = 10,face = "bold"))
d <- ggplot(data2,aes(x=Altitude,y=pos_33966039))+geom_point()+stat_smooth(method=lm)+geom_text(x=3000, y=0.012, label="r2=0.619",size=3.5)  +ggtitle("33966039")+ylab("Allele Frequency")+theme_classic()+theme(plot.title = element_text(colour = "red",size = 10,face = "bold"))
e <- ggplot(data2,aes(x=Altitude,y=pos_33964177))+geom_point()+stat_smooth(method=lm)+geom_text(x=3000, y=0.025, label="r2=0.609",size=3.5)  +ggtitle("33964177")+ylab("Allele Frequency")+theme_classic()+theme(plot.title = element_text(colour = "red",size = 10,face = "bold"))
f <- ggplot(data2,aes(x=Altitude,y=pos_34000269))+geom_point()+stat_smooth(method=lm)+geom_text(x=3000, y=0.015, label="r2=0.593",size=3.5)  +ggtitle("34000269")+ylab("Allele Frequency")+theme_classic()+theme(plot.title = element_text(colour = "red",size = 10,face = "bold"))

grid.arrange(a,b,c,d,e,f,nrow=2)


###snpEff
library(ggplot2)
library(ggseqlogo)
library(cowplot)
substrRight <- function(x){
  num = nchar(x)-2
  gsub('[0-9.]', '', substr(x, 6, nchar(x)))
}
setwd("/Users/guoyafei/Documents/01_Migration/07_ManuScript/NP/reviewer/snpEff")
pdf("Glu-1A.pdf", width = 13, height = 8)

  fasta_file = read.table("ppd.logo.seq",header=F,stringsAsFactors = F)
  fasta = fasta_file[,2]
  snpEff <- read.table("ppd.snpEff",header=F,stringsAsFactors = F,fill=TRUE,sep="\t")
  p1 <- ggseqlogo(fasta,method="prob")+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position="none")+
    #labs(title = pdf_tit[i,1])+
    theme(plot.title = element_text(hjust = 0.5,size = 15))
  Ref <- as.data.frame(snpEff$V3)
  colnames(Ref) <- "letter"
  Alt <- as.data.frame(snpEff$V4)
  colnames(Alt) <- "letter"
  if(dim(snpEff)[2] >6){
    Old <- as.data.frame(substring(snpEff$V7,3,5)) 
    colnames(Old) <- "letter"
    Pos <- as.data.frame(gsub('[a-zA-Z.*]', '', snpEff$V7))
    colnames(Pos) <- "letter"
    New <- as.data.frame(substrRight(snpEff$V7))
    colnames(New) <- "letter"
    all <- rbind(Ref,Alt,Old,Pos,New)
    num <- dim(Ref)[1]
    aln <- data.frame(
      letter = all,
      Type=rep(c("Ref","Alt","Old","Pos","New"), each=num),
      x=rep(1:num,5)
    )
    p2 <- ggplot(aln, aes(x, Type)) +
      geom_text(aes(label=letter,)) +
      scale_y_discrete(limits=c("New","Pos","Old","Alt","Ref"))+
      scale_x_continuous(breaks=1:10, expand = c(0.07, 0)) + xlab('') +
      theme_logo() +
      theme(legend.position = 'none', axis.text.x = element_blank(),)
  } else{
    all <- rbind(Ref,Alt)
    num <- dim(Ref)[1]
    aln <- data.frame(
      letter = all,
      species=rep(c("Ref","Alt"), each=num),
      x=rep(1:num,2)
    )
    p2 <- ggplot(aln, aes(x, Type)) +
      geom_text(aes(label=letter,),check_overlap = TRUE) +
      scale_y_discrete(limits=c("Alt","Ref"))+
      scale_x_continuous(breaks=1:10, expand = c(0.13, 0)) + xlab('') +
      theme_logo() +
      theme(legend.position = 'none', axis.text.x = element_blank())  
  }
  snpEff$impact <- NA
  snpEff[which(snpEff$V6=="HIGH"),8] <- 4
  snpEff[which(snpEff$V6=="MODERATE"),8] <- 3
  snpEff[which(snpEff$V6=="LOW"),8] <- 2
  snpEff[which(snpEff$V6=="MODIFILER"),8] <- 1
  bp_data <- data.frame(
    x=snpEff$V2,
    Impact=snpEff$impact
  )
  #bp_data$x <- c(1:dim(bp_data)[1])
  p3 <- ggplot(bp_data, aes(x, Impact))+
    geom_bar(stat = "identity", fill="grey")+
    theme_logo()+
    #scale_x_discrete( expand = c(0.055, 0))+
    xlab("")+
    theme(axis.text.x = element_text(angle = 90,size=15))
  suppressMessages(require(cowplot))
  p <- plot_grid(p1,p2,p3,ncol = 1, align = "v")
  print(p)

dev.off()


##reviewer V2 -------------------------------------------------------------
library(ggplot2)
setwd("/Users/guoyafei/Desktop/NP/reviewer")
#data <- read.table("result4.txt", header=T, stringsAsFactors = F)
#data <- read.table("average_result.txt", header=T, stringsAsFactors = F)
data <- read.table("snp_lm.txt", header=T, stringsAsFactors = F)
sub <- data[which(data$pop1.1 != " NA"), c(1,2,61,62)]
sub1 <- as.data.frame(sub[-1,])
colnames(sub1) <- c("ID1","POS","r2","P")
sub1$POS <- as.character(sub1$POS)
sub1$r2 <- as.numeric(sub1$r2)
sub1$P <- as.numeric(sub1$P)
sub1$color <- -log10(sub1$P)
  
ggplot(data = sub1, 
       aes(x = POS, y = r2,color=color)) +
  geom_point() +
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ s(log(x)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 8))+
  scale_color_gradient(low = "cyan",high = "red")
  #scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB"))


snp <- as.data.frame(t(data[c(1,10),3:60]))
colnames(snp) <- c("Altitude","Genotype")
a <- ggplot(snp,aes(x=Altitude,y=Genotype))+geom_point()+stat_smooth(method=lm)
+geom_text(x=1000, y=0.9, label="r2=0.778",size=3.5) +ggtitle("33953684")+ylab("Allele Frequency")+theme_classic()+theme(plot.title = element_text(colour = "red",size = 10,face = "bold"))


### reviewer----------------


##R: library prepare
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
#display.brewer.all()
setwd("/Users/guoyafei/Desktop/NP/PPD")
annotation_col <- read.table("ppd_anno.txt",header=T,stringsAsFactors = F,sep="\t")
#rownames(annotation_col) <- c(paste("AB_",c(001:212),sep=""),paste("ABD_",c(0001:1196),sep=""),paste("D_",c(001:220),sep=""))
rownames(annotation_col) <- annotation_col[,1]
all <- read.table("ppd.txt",header=T,stringsAsFactors = F)
#colnames(all) <- c("ID","REF","ALT",paste("ABD_",c(0001:1143),sep=""),paste("AB_",c(001:212),sep=""),paste("ABD_",c(1144:1196),sep=""))
data <- all[,4:61]
### ploidy,common name,region ------>
anno <- annotation_col[colnames(data),2,drop=FALSE]
anno <- anno[rownames(data),]
#anno2 <- anno[which(anno$Common_name != "OtherHexaploids"),]
#anno2$type <- anno2$Common_name
#anno2[which(anno2$Common_name == "Polish_wheat" |anno2$Common_name == "Rivet_wheat" | anno2$Common_name == "Persian_wheat" |anno2$Common_name == "Khorasan_wheat"  |anno2$Common_name ==  "Durum"),4] <- "Freethreshing-Tetraploids"
#anno3 <- anno2[,c(1,3,4)]
#anno3$type <- factor(anno3$type,levels = c("Wild_emmer","Domesticated_emmer","Ispahanicum","Georgian_wheat","Freethreshing-Tetraploids","Spelt","Macha","Club_wheat","Tibetan_semi_wild", "Xinjiang_wheat", "Vavilovii","Indian_dwarf_wheat", "Yunan_wheat","Landrace","Cultivar"))
anno2 <- anno[order(anno$RDA_Region),,drop=FALSE]

anno$RDA_Region <- factor(anno$RDA_Region,levels = c("EU","WA","IA","EA","SH"))

### extract data & plot ------>
data2 <- data[,which(colnames(data) %in% rownames(anno2))]

seq <- read.table("file58_ppd_SH_taxa.txt",header=F,stringsAsFactors = F)
seq2 <- seq[which(seq$V1 %in% colnames(data2)),]
data3 <- data2[,seq2]

#ann_color = list(
#  Taxa = c(AABB="orange", AABBDD="blue"),
#  Common_name = c(Wild_emmer = "#8C510A", Domesticated_emmer = "#DFC27D", Freethreshing = "#F6E8C3", EU="#66C2A5", WA= "#FC8D62",North1="#8DA0CB", North2 ="#E78AC3",South="#A6D854", Tibet="#FFD92F",Other ="#B3B3B3"))
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
pdf("test.pdf",width=12,height = 8 )
pheatmap(data3, show_rownames=F, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE,  annotation_names_col = F)
dev.off()


test <- read.table("test.txt",header=T,stringsAsFactors = F)
test2 <- t(test[,-c(1:2)])
test3 <- as.data.frame(test2)
colnames(test3) <- c("a","b")


a <- ggplot(test3,aes(x=a,y=b))+geom_point()+stat_smooth(method=lm)+geom_text(x=1000, y=0.9, label="r2=0.778",size=3.5) +ggtitle("33953684")+ylab("Allele Frequency")+theme_classic()+theme(plot.title = element_text(colour = "red",size = 10,face = "bold"))

ggplot(data = test2, 
       aes(x = a, y = b)) +
  geom_point() +
  #geom_smooth(se = TRUE, method = "gam", formula = y ~ s(log(x)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 8))
  #scale_color_gradient(low = "cyan",high = "red")




