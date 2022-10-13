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
#setwd("/Users/guoyafei/Documents/02_VmapIII/13_Haplotype")
annotation_col <- read.table("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/VMap3.info",header=T,stringsAsFactors = F,sep="\t")
rownames(annotation_col) <- c(paste("AB_",c(001:212),sep=""),paste("ABD_",c(0001:1196),sep=""),paste("D_",c(001:220),sep=""))

all <- read.table("8.297073952-297080886-geno.txt",header=T,stringsAsFactors = F)
colnames(all) <- c("ID","REF","ALT",paste("ABD_",c(0001:1143),sep=""),paste("AB_",c(001:212),sep=""),paste("ABD_",c(1144:1196),sep=""))
data <- all[,4:1411]
### ploidy,common name,region ------>
anno <- annotation_col[colnames(data),c(1,2,5),drop=FALSE]
anno2 <- anno[which(anno$Common_name != "OtherHexaploids"),]
anno2$type <- anno2$Common_name
anno2[which(anno2$Common_name == "Polish_wheat" |anno2$Common_name == "Rivet_wheat" | anno2$Common_name == "Persian_wheat" |anno2$Common_name == "Khorasan_wheat"  |anno2$Common_name ==  "Durum"),4] <- "Freethreshing-Tetraploids"
anno3 <- anno2[,c(1,3,4)]
anno3$type <- factor(anno3$type,levels = c("Wild_emmer","Domesticated_emmer","Ispahanicum","Georgian_wheat","Freethreshing-Tetraploids","Spelt","Macha","Club_wheat","Tibetan_semi_wild", "Xinjiang_wheat", "Vavilovii","Indian_dwarf_wheat", "Yunan_wheat","Landrace","Cultivar"))
anno4 <- anno3[order(anno3$type,anno3$Ctnt),]
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
pdf("test.pdf",width=12,height = 8 )
pheatmap(data3, show_rownames=FALSE, show_colnames=FALSE, color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation_col = anno3, annotation_names_col = F)
dev.off()
