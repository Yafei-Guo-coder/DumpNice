##haplotype
##shell
cat $1 |while read vcfchr chr from to
do
bcftools filter /data1/publicData/wheat/genotype/VMap/VMap3.2/VMap3.2/chr${vcfchr}_VMap3.2.vcf.gz --regions $chr:$from-$to > $chr.$from-$to.vcf
plink --vcf $chr.$from-$to.vcf --recode vcf-iid --out $chr.$from-$to.vcf-geno
grep -v "##" $chr.$from-$to.vcf-geno.vcf | awk '{$1=null;$2=null;$6=null;$7=null;$8=null;$9=null;print $0}'  | sed 's/  //' | sed 's/    //' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' | sed 's/\.\/\./-1/g' |less -LS > $chr.$from-$to.vcf-geno.txt
done
##R: library prepare
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
#display.brewer.all()
setwd("/Users/guoyafei/Documents/02_VmapIII/13_Haplotype/h12/gene_ud10k/")
annotation_col <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3.info",header=T,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- c(paste("AB_",c(001:212),sep=""),paste("ABD_",c(0001:1196),sep=""),paste("D_",c(001:220),sep=""),paste("A_",c(0001:0091),sep=""))

#AB
fileAB <- c("1.9543967-9566916.vcf-geno","8.288358354-288379382.vcf-geno","9.443813483-443835345.vcf-geno","10.108988281-109010657.vcf-geno","13.7284362-7307912.vcf-geno","16.105950688-105974974.vcf-geno","19.2392444-2422218.vcf-geno","20.12130022-12151979.vcf-geno","21.14378436-14399403.vcf-geno","26.28659895-28683828.vcf-geno","28.145835585-145859850.vcf-geno","31.9362298-9388625.vcf-geno","39.109692167-109715332.vcf-geno")
input <- paste(fileAB,".txt",sep="")
output <- paste(fileAB,".pdf",sep="")

#D
fileD <- c("5.5129131-5150186.vcf-geno","5.5346764-5367678.vcf-geno","5.74681-95607.vcf-geno","11.5658870-5684210.vcf-geno","12.49941839-49965803.vcf-geno","29.3561184-3581678.vcf-geno","29.3581526-3601972.vcf-geno","29.3599672-3620121.vcf-geno","41.6473046-6495989.vcf-geno","42.3374067-3394830.vcf-geno")
input <- paste(fileD,".txt",sep="")
output <- paste(fileD,".pdf",sep="")

for (i in c(1:length(input))) {
all <- read.table(input[i],header=T,stringsAsFactors = F)
#AB
#colnames(all) <- c("ID","REF","ALT",paste("ABD_",c(0001:1143),sep=""),paste("AB_",c(001:212),sep=""),paste("ABD_",c(1144:1196),sep=""))
#data <- all[,4:1411]
#D
colnames(all) <- c("ID","REF","ALT",paste("ABD_",c(0001:1143),sep=""),paste("D_",c(001:220),sep=""),paste("ABD_",c(1144:1196),sep=""))
data <- all[,4:1419]

anno <- annotation_col[colnames(data),c(2,3,6),drop=FALSE]
anno2 <- anno[which(anno$Common_name != "OtherHexaploids"),]
anno2$type <- anno2$Common_name
anno2[which(anno2$Common_name == "Polish_wheat" |anno2$Common_name == "Rivet_wheat" | anno2$Common_name == "Persian_wheat" |anno2$Common_name == "Khorasan_wheat"  |anno2$Common_name ==  "Durum"),4] <- "Freethreshing-Tetraploids"
anno3 <- anno2[,c(1,3,4)]
anno3$type <- factor(anno3$type,levels = c("Strangulata","Wild_emmer","Domesticated_emmer","Ispahanicum","Georgian_wheat","Freethreshing-Tetraploids","Spelt","Macha","Club_wheat","Tibetan_semi_wild", "Xinjiang_wheat", "Vavilovii","Indian_dwarf_wheat", "Yunan_wheat","Landrace","Cultivar"))
anno3$Ctnt <- factor(anno3$Ctnt,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM"))
anno4 <- anno3[order(anno3$type,anno3$Ctnt),]
data2 <- data[,which(colnames(data) %in% rownames(anno4))]
data3 <- data[,rownames(anno4)]
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
pdf(output[i],width=12,height = 8 )
pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation_col = anno3, annotation_names_col = F)
dev.off()
}

### bam depth statistic ------>
setwd("/Users/guoyafei/Documents/02_VmapIII/13_Haplotype")
data <- read.table("all.TE.depth.txt", header=T, stringsAsFactors = F)
name <- read.table("bam-vcf-depth.txt",header=T,stringsAsFactors = F)
rownames(name) <- name[,1]

sub <- as.data.frame(t(data[,which(colnames(data) %in% name$bam_ID)]))
sub$vcf <- "NA"
sub$depth <- "NA"
sub[,222] <- name[rownames(sub),2]
sub[,223] <- name[rownames(sub),3]
rownames(sub) <- sub[,222]
sub2 <- sub[,-c(222,223)]
depth <-  sub[,223]
out <- t(sub2)

#plot
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
annotation_col <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3.info",header=T,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- c(paste("AB_",c(001:212),sep=""),paste("ABD_",c(0001:1196),sep=""),paste("D_",c(001:220),sep=""),paste("A_",c(0001:0091),sep=""))
anno <- annotation_col[which(rownames(annotation_col) %in% colnames(out)),c(2,3,6),drop=FALSE]

anno2 <- anno[which(anno$Common_name != "OtherHexaploids"),]
anno2$type <- anno2$Common_name
anno2[which(anno2$Common_name == "Polish_wheat" |anno2$Common_name == "Rivet_wheat" | anno2$Common_name == "Persian_wheat" |anno2$Common_name == "Khorasan_wheat"  |anno2$Common_name ==  "Durum"),4] <- "Freethreshing-Tetraploids"
anno3 <- anno2[,c(1,3,4)]
anno3$type <- factor(anno3$type,levels = c("Wild einkorn","Domesticated einkorn","Urartu","Wild_emmer","Domesticated_emmer","Ispahanicum","Georgian_wheat","Freethreshing-Tetraploids","Spelt","Macha","Club_wheat","Tibetan_semi_wild", "Xinjiang_wheat", "Vavilovii","Indian_dwarf_wheat", "Yunan_wheat","Landrace","Cultivar"))
anno3$Ctnt <- factor(anno3$Ctnt,levels = c("WA","EU","SA","CA","EA","AF","OA","Nth_AM","Sth_AM"))

anno4 <- anno3[order(anno3$type,anno3$Ctnt),]

### extract data & plot ------>
data2 <- out[,which(colnames(out) %in% rownames(anno4))]
data3 <- out[,rownames(anno4)]
cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
pdf(output[i],width=12,height = 8 )
pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, cluster_col = F, cluster_row = FALSE, annotation_col = anno3, annotation_names_col = F)
dev.off()





#extract sample
data <- as.data.frame(t(out))
data$mean <- as.numeric(apply(out,2,mean))
sub <- data[which(data$mean >=1),1:221]
colnames(sub) <- paste("pos",c(1:221),sep="")

library(readxl)
a <- read_excel("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx", sheet=1, na='NA')
b <- a[,c(1:12,14)]
rownames(b) <- c(paste("AB_",c(1:212),sep=""),paste("ABD_",c(1:1196),sep=""),paste("D_",c(1:220),sep=""))
c <- b[which(rownames(b) %in% rownames(sub)),]

write.table(sub,"read_depth2.txt", quote=F,sep="\t")
write.table(c,"read_depth2_info.txt", quote=F,sep="\t",row.names = F)



##########
##R: library prepare
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
#display.brewer.all()
setwd("/Users/guoyafei/Desktop/vmap1.1/btr1-B/")
setwd("/Users/guoyafei/Documents/02_VmapIII/13_Haplotype/positive-gene")
annotation_col <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3.info",header=T,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- c(paste("AB_",c(001:212),sep=""),paste("ABD_",c(0001:1196),sep=""),paste("D_",c(001:220),sep=""),paste("A_",c(0001:0091),sep=""))
file <- c("11.33947048-33957048.vcf-geno","17.55774958-55784958.vcf-geno","17.56276700-56286700.vcf-geno","20.129917259-129927259.vcf-geno","20.130279578-130289578.vcf-geno","21.30357277-30367277.vcf-geno","21.30856268-30866268.vcf-geno","23.18458838-18468838.vcf-geno","23.18776062-18786062.vcf-geno","24.58272633-58282633.vcf-geno","24.58329232-58339232.vcf-geno","26.196891739-196901739.vcf-geno","26.244926539-244936539.vcf-geno","28.122425011-122435011.vcf-geno","28.206714543-206724543.vcf-geno","30.15270579-15280579.vcf-geno","30.69806790-69816790.vcf-geno","7.36928684-36938684.vcf-geno","13.65864056-65874644.vcf-geno","15.88966298-88976838.vcf-geno","17.55774918-55785576.vcf-geno")
input <- paste(file,".geno3",sep="")
output <- paste(file,".pdf",sep="")

all <- read.table("/Users/guoyafei/Desktop/vmap1.1/btr1-B/btr1-B.geno3",header=T,stringsAsFactors = F)
for (i in c(1:length(input))) {
  all <- read.table(input[i],header=T,stringsAsFactors = F)
  all <- read.table(input,header=T,stringsAsFactors = F)
  colnames(all) <- c("ID","REF","ALT",paste("ABD_",c(0001:1143),sep=""),paste("AB_",c(001:212),sep=""),paste("ABD_",c(1144:1196),sep=""))
  data <- all[,4:1411]
  ### ploidy,common name,region ------>
  anno <- annotation_col[colnames(data),c(2,3,6),drop=FALSE]
  anno2 <- anno[which(anno$Common_name != "OtherHexaploids"),]
  anno2$type <- anno2$Common_name
  anno2[which(anno2$Common_name == "Polish_wheat" |anno2$Common_name == "Rivet_wheat" | anno2$Common_name == "Persian_wheat" |anno2$Common_name == "Khorasan_wheat"  |anno2$Common_name ==  "Durum"),4] <- "Freethreshing Tetraploids"
  bw <- c("Spelt","Macha","Club_wheat","Tibetan_semi_wild", "Xinjiang_wheat", "Vavilovii","Indian_dwarf_wheat", "Yunan_wheat","Landrace","Cultivar")
  anno2[which(anno2$Common_name %in% bw),4] <- "Bread Wheat"
  my <- c("Wild_emmer_S","Wild_emmer_N","Domesticated_emmer","Freethreshing Tetraploids","Bread Wheat")
  anno3 <- anno2[which(anno2$type %in% my),c(1,3,4)]
  
  anno3$type <- factor(anno3$type,levels = c("Wild_emmer_S","Wild_emmer_N","Domesticated_emmer","Freethreshing Tetraploids","Bread Wheat"))
  #anno3$type <- factor(anno3$type,levels = c("Wild_emmer","Domesticated_emmer","Ispahanicum","Georgian_wheat","Freethreshing-Tetraploids","Spelt","Macha","Club_wheat","Tibetan_semi_wild", "Xinjiang_wheat", "Vavilovii","Indian_dwarf_wheat", "Yunan_wheat","Landrace","Cultivar"))
  anno4 <- anno3[order(anno3$type,anno3$Ctnt),]
  type <- anno3[,3,drop=F]
  ### ploidy,common name ------>
  #anno <- annotation_col[colnames(data),c(1,2),drop=FALSE]
  #anno2 <- anno[which(anno$Common_name != "OtherHexaploids"),]
  #anno2$type <- anno2$Common_name
  #anno2[which(anno2$Common_name == "Polish_wheat" |anno2$Common_name == "Rivet_wheat" | anno2$Common_name == "Persian_wheat" |anno2$Common_name == "Khorasan_wheat"  |anno2$Common_name ==  "Durum"),3] <- "Freethreshing-Tetraploids"
  #anno3 <- anno2[,c(1,3)]
  #anno3$type <- factor(anno3$type,levels = c("Wild_emmer","Domesticated_emmer","Ispahanicum","Georgian_wheat","Freethreshing-Tetraploids","Spelt","Macha","Club_wheat","Tibetan_semi_wild", "Xinjiang_wheat", "Vavilovii","Indian_dwarf_wheat", "Yunan_wheat","Landrace","Cultivar"))
  #anno4 <- anno3[order(anno3$type),]
  ### ploidy,region ------>
  #anno <- annotation_col[colnames(data),c(1,5),drop=FALSE]
  #anno2 <- anno[which(anno$Ctnt != "NA"),]
  #anno3 <- anno2[,c(1,2)]
  #anno4 <- anno3[order(anno3$Taxa,anno3$Ctnt),]
  
  ### extract data & plot ------>
  data2 <- data[,which(colnames(data) %in% rownames(anno4))]
  data3 <- data[,rownames(anno4)]
  #ann_color = list(
  #  Taxa = c(AABB="orange", AABBDD="blue"),
  #  Common_name = c(Wild_emmer = "#8C510A", Domesticated_emmer = "#DFC27D", Freethreshing = "#F6E8C3", EU="#66C2A5", WA= "#FC8D62",North1="#8DA0CB", North2 ="#E78AC3",South="#A6D854", Tibet="#FFD92F",Other ="#B3B3B3"))
  cols <- c(brewer.pal(11, "Set3")[c(5)],"#DEEBF7","#FEB24C","#BD0026")
  pdf(output,width=12,height = 8 )
  pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_cols = TRUE,cluster_row = FALSE, annotation_col = type, annotation_names_col = F)
  dev.off()
}


depth <- read.table("/Users/guoyafei/Desktop/btr1/depth.txt", header=T,stringsAsFactors = F)
rownames(depth) <- paste("pos_",c(88971298:88976838),sep="")
sub <- depth[1:5000,]
pheatmap(sub, show_rownames=T, show_colnames=FALSE,cluster_col = F, fontsize = 1,cluster_row = FALSE)

depth <- read.table("/Users/guoyafei/Desktop/btr1/all.selec.depth.txt", header=T,stringsAsFactors = F)
rownames(depth) <- paste("pos_",c(88971298:88976838),sep="")
sub <- depth[1:5000,]
pheatmap(sub, show_rownames=T, show_colnames=FALSE,cluster_col = F, fontsize = 1,cluster_row = FALSE)


depth <- read.table("/Users/guoyafei/Desktop/btr1/all.depth.txt", header=T,stringsAsFactors = F)
pheatmap(depth, show_rownames=T,show_colnames=FALSE,  cluster_col = F, cluster_row = FALSE)




