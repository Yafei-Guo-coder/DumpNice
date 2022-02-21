cat gene_region.txt |while read chr from to ID Name
do
if [ $chr != "Chr" ];then
bcftools filter /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225/chr${chr}.E6_Landrace_locate.vcf.gz --regions ${chr}:${from}-${to} > VCF/${ID}-${Name}.vcf
fi
done
cd VCF

for i in `ls *vcf`
do
vcftools --vcf $i --012
mv out.012.pos ${i::-4}.pos
rm out.012*
done
#此处删除没有位点的vcf文件
for i in `ls *pos`
do
head -n 1 $i | awk '{print $1}' 
echo $i
done |xargs -n2 > chr-file.txt

awk 'NR==FNR {a[$2]=$1;b[$2]=$2} NR!=FNR {if($1 in b) print $0,a[$1]}' chr-num.txt chr-file.txt > chr-file-num.txt

cat chr-file-num.txt | while read chr file num
do
bcftools view -R ${file} -S ${num}.taxa.txt /data2/yafei/003_Project3/Vmap1.1/E6/VCF/chr${chr}.E6all.vcf.gz >${num}.${file::-4}.pos.vcf
done

for i in `ls AB*pos.vcf`
do
vcftools --vcf $i --012
paste out.012.pos <(cat out.012 | datamash transpose | sed '1d' ) > AB.${i::-4}.pos.vcf.txt
cat out.012.indv | datamash transpose | awk '{print "File\tChr\tPos\t"$0}' > indi.txt
rm out.012*
done
for i in `ls AB.*pos.vcf.txt`
do
awk '{print FILENAME"\t"$0}' $i  | sed 's/.pos.pos.vcf.txt//g' |  sed 's/AB.AB.//g' > AB.${i::-4}.pos.vcf.txt2
done
cat indi.txt AB.*pos.vcf.txt2 > AB.all.pos.txt

for i in `ls D*pos.vcf`
do
vcftools --vcf $i --012
paste out.012.pos <(cat out.012 | datamash transpose | sed '1d' ) > D.${i::-4}.pos.vcf.txt
cat out.012.indv | datamash transpose | awk '{print "File\tChr\tPos\t"$0}' > indi.txt
rm out.012*
done
for i in `ls D.*pos.vcf.txt`
do
awk '{print FILENAME"\t"$0}' $i  | sed 's/.pos.pos.vcf.txt//g'|  sed 's/D.D.//g'  > D.${i::-4}.pos.vcf.txt2
done
cat indi.txt D.*pos.vcf.txt2 > D.all.pos.txt

rm AB.AB*
rm D.D.*
  
#画单倍型图
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(RColorBrewer)
cols <- c("#225EA8","#DEEBF7","#FEB24C","#BD0026")
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/11_Haplotype")
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/10_Gene")
annotation <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/select_taxa4.txt", header=T, stringsAsFactors = F, sep="\t")
anno <- annotation[,9,drop=FALSE]

#--------------------------------------------------------------------------------------------------------------------------------------
#AB
data <- read.table("AB.all.pos.txt", header=T, stringsAsFactors = F)
ann_color = list(
  heatmap = c(Wild_emmer = "#8C510A", Domesticated_emmer = "#DFC27D", Freethreshing = "#F6E8C3", EU="#66C2A5", WA= "#FC8D62",CA="#8DA0CB", EA ="#E78AC3",SA="#A6D854",Tibet="#FFD92F"))
for(i in names(table(data$File))){
  sub <- data[which(data$File==i),4:253]
  file <- paste(i,".pdf",sep="")
  pdf(file)
  AB_anno <- anno[which(rownames(anno) %in% colnames(sub)),1,drop=F]
  sub <- sub[, rownames(AB_anno)]
  pheatmap(sub, show_rownames=FALSE, show_colnames=FALSE,color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = AB_anno, annotation_colors = ann_color,annotation_names_col = F,main=i)
  dev.off()
}
#--------------------------------------------------------------------------------------------------------------------------------------
#D
data <- read.table("D.all.pos.txt", header=T, stringsAsFactors = F)
data <- read.table("ppd_1.txt",header=T,stringsAsFactors = F)
data <- read.table("ppd_3.txt",header=T,stringsAsFactors = F)
ann_color = list(
  #Growing_Habit = c(Facultative = "yellow", Spring="orange", Winter="blue"),
  heatmap = c(Strangulata = "#8C510A", EU="#66C2A5", WA= "#FC8D62",CA="#8DA0CB", EA ="#E78AC3",SA="#A6D854",Tibet="#FFD92F"))
for(i in names(table(data$File))){
  sub <- data[which(data$File==i),4:253]
  file <- paste(i,".pdf",sep="")
  pdf(file)
  AB_anno <- anno[which(rownames(anno) %in% colnames(sub)),1,drop=F]
  sub <- sub[, rownames(AB_anno)]
  pheatmap(sub, show_rownames=FALSE, show_colnames=FALSE,color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = AB_anno, annotation_colors = ann_color,annotation_names_col = F,main=i)
  dev.off()
}
#单个基因画图
sub <- data
AB_anno <- anno[which(rownames(anno) %in% colnames(sub)),1,drop=F]
sub <- sub[, rownames(AB_anno)]
pheatmap(sub, show_rownames=FALSE, show_colnames=FALSE,color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = AB_anno, annotation_colors = ann_color,annotation_names_col = F)

anno <- annotation[,c(1,8),drop=FALSE]
sub <- data
AB_anno <- anno[which(anno$VMap3 %in% colnames(sub)),c(1,2),drop=F]
rownames(AB_anno) <- AB_anno$VMap3
AB_anno <- AB_anno[,2,drop=F]
sub <- sub[, rownames(AB_anno)]
pheatmap(sub, show_rownames=FALSE, show_colnames=FALSE,color = cols, legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_col = F, cluster_row = FALSE, annotation = AB_anno, annotation_colors = ann_color,annotation_names_col = F)

#annotation_colors = ann_colors
#接下来可以使用06_Qmatrix_PieMap.r来画在地图上的单倍型分布。
#plot haplotype heatmap
brewer.pal(9, "YlGnBu")[c(7)]
"#225EA8"
brewer.pal(9, "Blues")[c(2)]
"#DEEBF7"
brewer.pal(9, "YlOrRd")[c(4)]
"#FEB24C"
brewer.pal(9, "YlOrRd")[c(8)]
"#BD0026"





  
  



