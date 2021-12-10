java -Xmx200g -Xms512m -jar /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile chr36.maf0.01.recode.vcf -inputformat VCF -outputfile chr36.maf0.01.env -outputformat BAYENV -spid VCF_BAYENV.spid#!/bin/bash
#just a small bash script to calculate BFs for all SNPs from SNPFILE
#please copy this script into the same directory as bayenv and execute it there
#please see the Bayenv2 manual for details about usage
#make this script executable (chmod +x calc_bf.sh)
#Usage: ./calc_bf.sh <Name of your SNPSFILE> <Name of your ENVFILE> <Name of your MATFILE> <Nuber of populations> <Number of MCMC iterations> <Number of environmental factors>
SNPFILE=$1
ENVFILE=$2
MATFILE=$3
POPNUM=$4
ITNUM=$5
ENVNUM=$6

split -a 10 -l 2 $SNPFILE snp_batch

for f in $(ls snp_batch*)
do
./bayenv2 -i $f -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t
done

rm -f snp_batch*

#Working directory:
#204:yafei:/data1/home/yafei/003_Project3/Structure/E6_Landrace_locate_225/bayenv
#准备环境变量文件
awk 'NR==FNR{a[$1]=$2;b[$1]=$1}NR!=FNR{if($1 in b) print $0"\t"a[$1]}' pop.txt 225env.txt | sort -k24,24 > merge_env.txt
datamash groupby 24 mean 2 mean 3 mean 4 mean 5 mean 6 mean 7 mean 8 mean 9 mean 10 mean 11 mean 12 mean 13 mean 14 mean 15 mean 16 mean 17 mean 18 mean 19 mean 20 mean 21 mean 22 mean 23 < merge_env.txt | datamash transpose > format1.txt
#R
data <- read.table("format1.txt",header=T,stringsAsFactors=F)
m <- apply(data,1,mean)
s <- apply(data,1,sd)
sub <- (data-m)/s
write.table(sub,"format2.txt", sep="\t", quote=F,row.names=F)
#去掉第一行就是bayenv的输入文件。
#调整pop的顺序。
#awk '{print $1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$2"\t"$3"\t"$4"\t"$5}' format2.txt  > format3.txt
#shell
#datamash transpose < format2.txt > format3.txt

#准备基因型文件
vcftools --vcf chr36.E6_Landrace_locate.vcf --maf 0.01 --recode --recode-INFO-all --out chr36.maf0.01
java -Xmx200g -Xms512m -jar /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile chr36.maf0.01.recode.vcf -inputformat VCF -outputfile chr36.maf0.01.env -outputformat BAYENV -spid VCF_BAYENV.spid

#运行bayenv
#matrix estimation
#使用筛选过LD的VCF文件：50 10 0.2
java -Xmx200g -Xms512m -jar /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile LD/chr36.in.vcf -inputformat VCF -outputfile chr36.LD.env -outputformat BAYENV -spid VCF_BAYENV.spid
bayenv2 -i chr36.LD.env -s samplesize.txt -p 5 -k 100000 -r 83556 -o chr36.matrix,l

#环境变量相关性估计
#./calc_bf.sh SNPSFILE ENVIRONFILE MATRIXFILE NUMPOPS NUMITER NUMENVIRON
./calc_bf.sh chr36.maf0.01.env format2.txt chr36.matrix 5 100000 22
#bayenv2 -i rs316 -m hgdp_matrix_1 -e PCs.env -p 52 -k 1000 -n 4 -t -r 42 -o out_correlation


./calc_bf.sh chr1.1-20000001.envgenofile 13pop.env matrix/A.matrix 13 10000 22
bayenv2 -i test.txt -m matrix/A.matrix -e 13pop.env -p 13 -k 10000 -n 22 -t -r 42 -o out_correlation
./calc_bf.sh Aenvgenofile/chr8.80000001-100000001.envgenofile ENVBAY/13pop.env matrix/A.matrix 13 10000 22
#根据环境变量对样本进行聚类
#使用05_Sample_Cluster.r: working directory: /Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv
library(cluster)
library(factoextra)
#聚类
#根据经纬度给样本聚类
data <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv/225env.txt", header=F,stringsAsFactors = F)
colname <- c("elevation","temp1","temp2","temp3","temp4","temp5","temp6","temp7","temp8","temp9","temp10","temp11","prec1","prec2","prec3","prec4","prec5","prec6","prec7","prec8","Latitude","Logititude")
rownames(data) <- data[,1]
data2 <- data[!is.na(data$V6),-1]
data2 <- data2[which(data2$V23 > -40),]
colnames(data2) <- colname
#按列进行标准化并聚类
df = scale(data2,center = T,scale = T)
colnames(df) <- colname

------------------------------kmeans聚类---------
#确定应该分几个cluster
data2<- data2[,c(1:22)]
mydata <- data2
mydata <- df
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
     
#kmeans聚类，标准化后
data2<- data2[,c(1:22)]
km <- kmeans(df,13,iter.max = 10000)
km <- kmeans(data2, 11,iter.max = 5000) #用于画地图
fviz_cluster(km, data = df,
  #palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
  ellipse.type = "euclid",
  star.plot = T, 
  repel = TRUE,
  ggtheme = theme_minimal()
)
data2$type <- km$cluster
-------------------------------层次聚类---------
#求样本之间两两相似性，层次聚类
data2<- data2[,c(1:22)]
result <- dist(data2, method = "euclidean")
result2 <- dist(df, method = "euclidean")
#使用指定距离来计算数据矩阵行之间的距离
#euclidean：欧几里得距离
#maximum：最大距离
#manhattan：绝对距离
#canberra：堪培拉距离
#minkowski：闵可夫斯基距离
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
#Ward: 最小方差方法旨在寻找紧凑的球形簇的完整的联动方法找到相似集群。
#有两种不同的算法，"ward.D"（相当于只沃德选择"ward"不执行沃德（1963）聚类准则）& "ward.D2"实现了标准（Murtagh的及Legendre 2014）在集群更新之前对差异进行平方。注意agnes(*, method="ward")对应于hclust(*, "ward.D2").
#median和centroid在不导致单调的距离测量，或者等效产生的树状图可以具有所谓的倒置或颠倒。
data2$type <- cutree(result_hc, k=11)

------------------------计算类群环境变量---------
lat_mean <- tapply(data2[,21],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,22],data2$type,mean,na.rm = TRUE)
data2$cluster1 <- NA
data2$cluster2 <- NA
for(i in 1:219) {
  for(j in 1:8){
    if(data2[i,23] == j ){
      data2[i,24] <- as.numeric(lat_mean[j])
      data2[i,25] <- as.numeric(lon_mean[j])
    } 
  }
}
write.table(data2,"13_cluster.txt",sep="\t",row.names = T,quote=F)

-----------------------------地图上展示聚类结果--------
#new <- data.frame(cluster1=km$centers[,21], cluster2=km$centers[,22],size = km$size)
data <- read.table("13_cluster.txt",header=T,stringsAsFactors = F)
library(maps)
library(ggplot2)
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70") 
mp<-ggplot()+mapworld+ylim(-50,80)
mp_13<-mp+geom_point(aes(x=data2$Logititude, y=data2$Latitude,color = as.factor(data2$type)))+
  scale_size(range=c(1,1))+ 
  theme_classic()
  
data <- read.table("type.txt",header=T,stringsAsFactors = F)
mp_13<- mp+geom_point(aes(x=data$cluster2, y=data$cluster1,size=data$Type))+
  #scale_size(range=c(1,1))+ 
  theme_classic()
mp_13

#-------------------------------结果处理---------------------------------
Workinig directory:
  SNPsfile：203@xuebo /data1/home/xuebo/Projects/Speciation/BAYENV/genofile/Dlineage2
  结果文件：203@xuebo /data1/home/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageD2
  204@yafei /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY
120327.bayenv.top5.bed是top5的snp位点。
18203.merge50k.bayenv.top5.bed是合并距离小于50k的snp，并且把单个的snp去掉。
然后在xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%
step1:提取XPCLR和bayenv重叠的XPCLR区域
for i in `ls *_D.top5.bed`
do
bedtools intersect -a $i -b 18203.merge50k.bayenv.top5.bed  -wa > bayenv_over/$i
done
step2:定位已注释基因:gff
for i in `ls *bed`
do
bedtools intersect -a ../../gene_v1.1_Lulab.gff3 -b $i -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > ${i::-3}gff.gene
done

#定位已克隆基因:cloned gene 以及定位抗病基因:nlr gene
for i in `ls *gff.gene`
do
awk -F"ID=" '{print $2}' $i > ${i::-8}gene2
done
for i in `ls *.gene2`
do
grep -w -f $i ../../all_cloned_gene.txt > ${i::-5}cloned.gene
#grep -w -f $i ../nlr_gene.txt > nlr/${i::-5}cloned.gene
done
rm *gene2

ls *cloned.gene |xargs -n1 > cloned_gene.txt
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done | awk '{print $1}' | awk -F"_smooth" '{print $1}'|uniq > file_prefix.txt
#change file format to plot heatmap(A,B,D lineage seperate)
#A lineage
ls *A.top5.cloned.gene |xargs -n1 > A_cloned_gene.txt
sed 's/$/_smooth_A.top5.cloned.gene/' file_prefix.txt > A_file.txt
for i in `cat A_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > A_gene.txt

#B lineage
ls *B.top5.cloned.gene |xargs -n1 > B_cloned_gene.txt
sed 's/$/_smooth_B.top5.cloned.gene/' file_prefix.txt> B_file.txt
for i in `cat B_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > B_gene.txt

#D lineage 
ls *D.top5.cloned.gene |xargs -n1 > D_cloned_gene.txt
sed 's/$/_smooth_D.top5.cloned.gene/' file_prefix.txt > D_file.txt
for i in `cat D_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > D_gene.txt

#R
v <- "D"
for ( i in c(1)){
  filename <- paste(v[i],"_file.txt",sep="")
  genename <- paste(v[i],"_gene.txt",sep="")
  file <- read.table(filename,header=F,stringsAsFactors=F)
  gene <- read.table(genename,header=F,stringsAsFactors=F)
  files <- file[,1]
  genes <- gene[,1]
  x=vector()
  for (j in files){
    a<-paste(j,genes,sep=" ")
    x <- c(x,a)
  }
  out <- paste(v[i],"_file_gene_mode.txt",sep="")
  write.table(x,out,quote=F,row.names=F,col.names=F,sep="\t")
}

#shell
#A lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep A.top5.cloned.gene > A_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' A_postive_file_gene_mode.txt A_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > A_heatmap_format1.txt
#B lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep B.top5.cloned.gene > B_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' B_postive_file_gene_mode.txt B_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > B_heatmap_format1.txt
#D lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep D.top5.cloned.gene > D_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' D_postive_file_gene_mode.txt D_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > D_heatmap_format1.txt

step3:提取Sr45和Sr33的基因区上下游50k snp，重新做bayenv。









