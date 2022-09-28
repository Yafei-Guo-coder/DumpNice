#--------------------step1:运行bayenv---------------------
/data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene
java -Xmx200g -Xms512m -jar /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile chr36.maf0.01.recode.vcf -inputformat VCF -outputfile chr36.maf0.01.env -outputformat BAYENV -spid VCF_BAYENV.spid
#!/bin/bash
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
#204:yafei:/data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY
#准备环境变量文件
awk 'NR==FNR{a[$1]=$2;b[$1]=$1}NR!=FNR{if($1 in b) print $0"\t"a[$1]}' pop.txt 225env.txt | sort -k24,24 > merge_env.txt
datamash groupby 24 mean 2 mean 3 mean 4 mean 5 mean 6 mean 7 mean 8 mean 9 mean 10 mean 11 mean 12 mean 13 mean 14 mean 15 mean 16 mean 17 mean 18 mean 19 mean 20 mean 21 mean 22 mean 23 < merge_env.txt | datamash transpose > format1.txt
#R
data <- read.table("format1.txt",header=T,stringsAsFactors=F)
m <- apply(data,1,mean)
s <- apply(data,1,sd)
sub <- (data-m)/s
write.table(sub,"format2.txt", sep="\t", quote=F,row.names=F)
#去掉第一行就是bayenv的输入文件
#调整pop的顺序
#awk '{print $1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$2"\t"$3"\t"$4"\t"$5}' format2.txt  > format3.txt
#shell
#datamash transpose < format2.txt > format3.txt
#准备基因型文件
/data1/home/xuebo/software/PGDSpider/
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
data <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv/V1/225env.txt", header=F,stringsAsFactors = F)
colname <- c("elevation","temp1","temp2","temp3","temp4","temp5","temp6","temp7","temp8","temp9","temp10","temp11","prec1","prec2","prec3","prec4","prec5","prec6","prec7","prec8","Latitude","Logititude")
rownames(data) <- data[,1]
data2 <- data[!is.na(data$V6),-1]
data2 <- data2[which(data2$V23 > -40),]
colnames(data2) <- colname
#按列进行标准化并聚类
df = scale(data2,center = T,scale = T)
colnames(df) <- colname
#------------------------------kmeans聚类---------
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
km <- kmeans(data2, 13,iter.max = 10000) #用于画地图
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

----------------------------------------地图上展示聚类结果---------------------------------
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
#-----------------step2:提取所有克隆基因区上下游50k snp，重新做bayenv----------------
Working directory: 
  yafei@204: /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/Gene50k
awk '{print $1"\t"$2-50000"\t"$3+50000}' clone.gene.txt > gene.50k.txt
219个样本的VCF文件：/data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/VCF/chr1.lineage_Landrace_225_noAM_220_maf001_LD.vcf

1.区分A,B,D亚基因组
awk '{output="chr"$1".txt"; print $0 > output}' gene.50k.txt 
for i in {1..42}
do
bcftools view -R pos/chr${i}.txt /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225/chr${i}.E6_Landrace_locate.vcf.gz -o Chr/chr{i}.gene.50k.vcf &
done
2.从225个样本的vcf文件中提取出这些位点的snp
#yafei@204:/data1/home/yafei/003_Project3/Structure/E6_Landrace_locate_225/Lineage
bcftools view -R A.gene.50k.txt /data1/home/yafei/003_Project3/Structure/E6_Landrace_locate_225/Lineage/ABlineage.E6_Landrace_locate.vcf.gz -o A.gene.50k.vcf.gz -O z 
bcftools view -R B.gene.50k.txt /data1/home/yafei/003_Project3/Structure/E6_Landrace_locate_225/Lineage/ABlineage.E6_Landrace_locate.vcf.gz -o B.gene.50k.vcf.gz -O z 
bcftools view -R D.gene.50k.txt /data1/home/yafei/003_Project3/Structure/E6_Landrace_locate_225/Lineage/Dlineage.E6_Landrace_locate.vcf.gz -o D.gene.50k.vcf.gz -O z 
3.格式转换
#!/bin/bash
for chr in {1..42}
do
	java -jar -Xmx500g -Xms100g /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile Chr/chr${chr}.gene.50k.vcf -inputformat VCF -outputfile envgenofile/chr${chr}.envgenofile -outputformat BAYENV -spid spid/VCF_BAYENV_chr${chr}.spid
done
4.跑bayenv
./calc_bf.sh lineage/Alineage.envgenofile 13pop.env matrix/A.matrix 13 100000 22 Alineage_out
./calc_bf.sh lineage/Blineage.envgenofile 13pop.env matrix/B.matrix 13 100000 22 Blineage_out
./calc_bf.sh lineage/Dlineage.envgenofile 13pop.env matrix/D.matrix 13 100000 22 Dlineage_out
5.结果解析

--------------------------------------------------------重新统计---------------------------------------------------------------
204:/data2/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageA
204:/data2/xuebo/Projects/Speciation/BAYENV/genofile/Alineage2
这个里面是A lineage的结果和原始的计算文件
203:/data1/home/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageB
203:/data1/home/xuebo/Projects/Speciation/BAYENV/genofile/Blineage
这个里面是B lineage的结果和原始的计算文件
203:/data1/home/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageD2
203:/data1/home/xuebo/Projects/Speciation/BAYENV/genofile/Dlineage2
这个里面是D lineage的结果和原始的计算文件

203:/data1/home/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageA_New
203:/data1/home/xuebo/Projects/Speciation/BAYENV/genofile/AlineageNew
这是重新做的Alineage的情况
204:/data2/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageB_New
204:/data2/xuebo/Projects/Speciation/BAYENV/genofile/BlineageNew
这是重新做的Blineage的情况
204:/data2/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageD_New
204:/data2/xuebo/Projects/Speciation/BAYENV/genofile/DlineageNew
这是重新做的Dlineage的情况

204:/data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/Gene50k

----------------------------------------------------------结果文件--------------------------------------------------------------
yafei@204:/data1/home/yafei/003_Project3/Structure/bayenv/Result

总共分析的位点数：(三批数据)
A:444085; B:498827; D:482174
A:208479; B:222484; D:161206
A:1197; B:975; D:1090

总共得到的位点数：A:329977(16498) B:355505(17775) D:401549(20077)
A:235606; B:276343; D:321195;
A:93174; B:78187; D:79264;
A:1197; B:975; D:1090;

--------------------------------------基因上下游5k的结果统计-------------------------
204:/data1/home/yafei/003_Project3/bayenv/13pop/544Gene

-----------------------------------------20211221 新的结果--------------------------------------------
原始文件：xuebo
203:/data1/home/xuebo/Projects/Speciation/BAYENV/BAYENV_V2/bayenv2_out_lineageB
203:/data1/home/xuebo/Projects/Speciation/BAYENV/BAYENV_V2/bayenv2_out_lineageD
203:/data1/home/xuebo/Projects/Speciation/BAYENV/BAYENV_V2/genofile/B
203:/data1/home/xuebo/Projects/Speciation/BAYENV/BAYENV_V2/genofile/D
204:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V2/bayenv2_out_lineageA
204:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V2/genofile/A
204：/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225/Landrace_225_noAM_220_maf001/lineage/Alineage_Landrace_225_noAM_220_maf001_LD.vcf
204：/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225/Landrace_225_noAM_220_maf001/lineage/Blineage_Landrace_225_noAM_220_maf001_LD.vcf
204：/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225/Landrace_225_noAM_220_maf001/lineage/Dlineage_Landrace_225_noAM_220_maf001_LD.vcf

-----------------------------------------20211224 新的结果-------------------------
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/bayenv2_out_lineageA
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/bayenv2_out_lineageB
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/genofile/Alineage_bayenv_pop5_500K.vcf
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/genofile/Blineage_bayenv_pop5_500K.vcf
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/genofile/Alineage.pop5.envgenofile
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/genofile/Blineage.pop5.envgenofile

-----------------------------------------20220107 新的结果-------------------------
204:/data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2
#-----------------step3:统计富集情况--------
for i in {"A","B","D"}
do
awk 'NR==FNR{a[$4]=$1;b[$4]=$2;c[$4]=$3;d[$4]=$4}NR!=FNR{if($1 in d) print a[$1],b[$1],c[$1],$0}' ${i}.pos.bed ${i}_out.bf > ${i}_all.bf
cat /data1/home/yafei/003_Project3/bayenv/13pop/544Gene/${i}_all.bf ${i}_all.bf > Add544Gene_V2/${i}_all.bf
done

cd Add544Gene_V2

for chr in {"A","B","D"}
do
for i in {5..26}
do 
a=`sort -k${i},${i}g ${chr}_all.bf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$c}' c=$i | grep -v inf | wc -l`
let b=a/20
sort -k${i},${i}g ${chr}_all.bf | grep -v inf | tail -n ${b} | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$d}' d=$i > ${chr}bio${i}.top5.txt
done
done

mkdir 22bio Temp Prec All
mv *top5.txt 22bio

for i in {"A","B","D"}
do
for j in {6..16}
do
cat 22bio/${i}bio${j}.top5.txt >> ${i}.Temperature.bf
done
done 
for i in {"A","B","D"}
do
sort -k1,1n -k2,2n ${i}.Temperature.bf | awk '{print $1"\t"$2"\t"$3}' | uniq  > Temp/${i}temp.bayenv.top5.bed
done 
rm *Temperature.bf

for i in {"A","B","D"}
do
for j in {17..24}
do
cat 22bio/${i}bio${j}.top5.txt >> ${i}.Precipitation.bf
done
done 
for i in {"A","B","D"}
do
sort -k1,1n -k2,2n  ${i}.Precipitation.bf | awk '{print $1"\t"$2"\t"$3}' | uniq  > Prec/${i}prec.bayenv.top5.bed
done
rm *Precipitation.bf

for i in {"A","B","D"}
do
cat 22bio/${i}bio*.top5.txt >> ${i}.all.bf
done 
for i in {"A","B","D"}
do
sort -k1,1n -k2,2n  ${i}.all.bf | awk '{print $1"\t"$2"\t"$3}' | uniq  > All/${i}log.bayenv.top5.bed
done
rm *all.bf

#-----------------step4:分 region 和 site 定位基因----------------------
mkdir bayenv_xpclr_site
mkdir bayenv_xpclr_region

#all
for i in `ls *smooth_A.top5.bed`
do
bedtools intersect -b $i -a Aprec.bayenv.top5.bed  -wa |sort | uniq > bayenv_xpclr_site/$i
done
for i in `ls *smooth_B.top5.bed`
do
bedtools intersect -b $i -a Bprec.bayenv.top5.bed -wa|sort | uniq > bayenv_xpclr_site/$i
done
for i in `ls *smooth_D.top5.bed`
do
bedtools intersect -b $i -a Dprec.bayenv.top5.bed -wa |sort | uniq > bayenv_xpclr_site/$i
done

#all
for i in `ls *smooth_A.top5.bed`
do
bedtools intersect -a $i -b Aprec.bayenv.top5.bed  -wa |sort | uniq > bayenv_xpclr_region/$i
done
for i in `ls *smooth_B.top5.bed`
do
bedtools intersect -a $i -b Bprec.bayenv.top5.bed -wa|sort | uniq > bayenv_xpclr_region/$i
done
for i in `ls *smooth_D.top5.bed`
do
bedtools intersect -a $i -b Dprec.bayenv.top5.bed -wa |sort | uniq > bayenv_xpclr_region/$i
done
#site
cd /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2/Prec/bayenv_xpclr_site
mkdir nlr
for i in `ls *bed`
do
bedtools intersect -a /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2/gene_v1.1_Lulab.gff3 -b $i -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > ${i::-3}gff.gene
done

#定位已克隆基因:cloned gene 以及定位抗病基因:nlr gene
for i in `ls *gff.gene`
do
awk -F"ID=" '{print $2}' $i > ${i::-8}gene2
done
for i in `ls *.gene2`
do
grep -w -f $i /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2/cloned.gene.txt > ${i::-5}cloned.gene
grep -w -f $i /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2/nlr_gene.txt > nlr/${i::-5}cloned.gene
done
rm *gene2
#clone_gene
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
Rscript ../../gene_mode.r
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

Rscript ../../heatmap_format2.r
#nlr
cd nlr
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
Rscript ../../../gene_mode.r
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

Rscript ../../../heatmap_format2.r

#region
cd /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2/Prec/bayenv_xpclr_region
mkdir nlr
for i in `ls *bed`
do
bedtools intersect -a /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2/gene_v1.1_Lulab.gff3 -b $i -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > ${i::-3}gff.gene
done

#定位已克隆基因:cloned gene 以及定位抗病基因:nlr gene
for i in `ls *gff.gene`
do
awk -F"ID=" '{print $2}' $i > ${i::-8}gene2
done
for i in `ls *.gene2`
do
grep -w -f $i /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2/cloned.gene.txt > ${i::-5}cloned.gene
grep -w -f $i /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2/nlr_gene.txt > nlr/${i::-5}cloned.gene
done
rm *gene2
#clone_gene
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

Rscript ../../gene_mode.r
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

Rscript ../../heatmap_format2.r
#nlr
cd nlr
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
Rscript ../../../gene_mode.r
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

Rscript ../../../heatmap_format2.r
  
#-----------------step5:画曼哈顿图-----
本地：/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv/V4/manhuttan
服务器：204:/data1/home/yafei/003_Project3/bayenv/run1/Add544Gene_V2_Top5/all_bio

library(ggplot2)
library(gridExtra)
library(tidyverse)
centro <- read.table("centro.txt",header=F,stringsAsFactors=F,row.names=1)
names <- paste(rep("A",each=19),"bio",c(6:24),".all.txt",sep="")
Bio <- paste("BIO",c(1:19),sep="")
p<-list()

for (i in c(1:19)){
  cen <- centro[c(1,4,7,10,13,16,19),1,drop=F]
  data <- read.table(names[i],header=F,stringsAsFactors=F)
  data[which(data$V3>5),3] <- data[which(data$V3>5),3]/500
  data$V3 <- (data$V3-mean(data$V3))/sd(data$V3)
  #thresh <- data[order(data$V3),][nrow(data)*0.95,3]
  gwasResults <- data[sample(c(1:nrow(data)), 10000),]
  thresh <- gwasResults[order(gwasResults$V3),][nrow(gwasResults)*0.95,3]
  colnames(gwasResults) <- c("CHR", "BP","P")
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
  len <- don %>% group_by(CHR) %>% summarize(center=max(BPcum))
  asx <- c(0,len$center[1:6])
  cen$V2 <- cen$V2+asx
  p[[i]] <- ggplot(don, aes(x=BPcum, y=P)) +
    # Show all points
    #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.1) +
    geom_line(aes(color=as.factor(CHR)))+
    scale_color_manual(values = rep(c("grey","skyblue","grey","skyblue"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
    # Custom the theme:
    #geom_vline(xintercept = 33952048, color = 'skyblue', size = 0.5) + 
    #geom_vline(xintercept = 33956269, color = 'skyblue', size = 0.5) + 
    geom_hline(yintercept = thresh, color = 'red', size = 0.5) + 
    geom_point(data=cen,aes(V2,0),color="#BF812D")+
    #geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    scale_y_continuous(limits = c(-2,10))+
    theme_classic() +
    ylab("BF")+
    ggtitle(Bio[i])+
    xlab("A")+
    theme( 
      plot.title = element_text(color="black", size=15, face="bold"),
      legend.position="none",
      #panel.border = element_blank(),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=11),
      axis.text.y=element_text(size=11),
      axis.title.y=element_text(size = 14),
      axis.title.x=element_text(size = 14),
    )
}

names <- paste(rep(c("B"),each=19),"bio",c(6:24),".all.txt",sep="")
#p<-list()
for (i in c(1:19)){
  cen <- centro[c(2,5,8,11,14,17,20),1,drop=F]
  data <- read.table(names[i],header=F,stringsAsFactors=F)
  data[which(data$V3>100),3] <- data[which(data$V3>100),3]/4
  data$V3 <- (data$V3-mean(data$V3))/sd(data$V3)
  #thresh <- data[order(data$V3),][nrow(data)*0.95,3]
  gwasResults <- data[sample(c(1:nrow(data)), 10000),]
  thresh <- gwasResults[order(gwasResults$V3),][nrow(gwasResults)*0.95,3]
  colnames(gwasResults) <- c("CHR", "BP","P")
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
  len <- don %>% group_by(CHR) %>% summarize(center=max(BPcum))
  asx <- c(0,len$center[1:6])
  cen$V2 <- cen$V2+asx
  p[[i+19]] <- ggplot(don, aes(x=BPcum, y=P)) +
    # Show all points
    #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.1) +
    geom_line(aes(color=as.factor(CHR)))+
    scale_color_manual(values = rep(c("grey","#FDCDAC","grey","#FDCDAC"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
    # Custom the theme:
    #geom_vline(xintercept = 33952048, color = 'skyblue', size = 0.5) + 
    #geom_vline(xintercept = 33956269, color = 'skyblue', size = 0.5) + 
    geom_hline(yintercept = thresh, color = 'red', size = 0.5) + 
    geom_point(data=cen,aes(V2,0),color="#BF812D")+
    #geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    scale_y_continuous(limits = c(-2,5))+
    theme_classic() +
    #ggtitle(Bio[i+19])+
    ylab("BF")+
    xlab("B")+
    theme( 
      legend.position="none",
      #plot.title = element_text(color="black", size=15, face="bold"),
      #panel.border = element_blank(),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=11),
      axis.text.y=element_text(size=11),
      axis.title.y=element_text(size = 14),
      axis.title.x=element_text(size = 14),
    )
}


names <- paste(rep(c("D"),each=19),"bio",c(6:24),".all.txt",sep="")
#p<-list()
for (i in c(1:19)){
  cen <- centro[c(3,6,9,12,15,18,21),1,drop=F]
  data <- read.table(names[i],header=F,stringsAsFactors=F)
  data[which(data$V3>100),3] <- data[which(data$V3>100),3]/4
  data$V3 <- (data$V3-mean(data$V3))/sd(data$V3)
  #thresh <- data[order(data$V3),][nrow(data)*0.95,3]
  gwasResults <- data[sample(c(1:nrow(data)), 10000),]
  thresh <- gwasResults[order(gwasResults$V3),][nrow(gwasResults)*0.95,3]
  colnames(gwasResults) <- c("CHR", "BP","P")
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
  len <- don %>% group_by(CHR) %>% summarize(center=max(BPcum))
  asx <- c(0,len$center[1:6])
  cen$V2 <- cen$V2+asx  
  p[[i+38]] <- ggplot(don, aes(x=BPcum, y=P)) +
    # Show all points
    #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.1) +
    geom_line(aes(color=as.factor(CHR)))+
    scale_color_manual(values = rep(c("grey","#CCEBC5","grey","#CCEBC5"), 7)) +
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Add highlighted points
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    # Add label using ggrepel to avoid overlapping
    #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
    # Custom the theme:
    #geom_vline(xintercept = 33952048, color = 'skyblue', size = 0.5) + 
    #geom_vline(xintercept = 33956269, color = 'skyblue', size = 0.5) + 
    geom_hline(yintercept = thresh, color = 'red', size = 0.5) + 
    geom_point(data=cen,aes(V2,0),color="#BF812D")+
    #geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    scale_y_continuous(limits = c(-2,5))+
    theme_classic() +
    ylab("BF")+  
    xlab("D")+
    theme( 
      legend.position="none",
      #panel.border = element_blank(),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=11),
      axis.text.y=element_text(size=11),
      axis.title.y=element_text(size = 14),
      axis.title.x=element_text(size = 14),
    )
}

#plot_name <- c(paste("Temp",c(1:11),".pdf",sep=""),paste("Prec",c(1:8),".pdf",sep=""))
#count=0
#for (j in c(0:18)){
#  count=count+1
#  num <- vector()
#  for (i in seq(1,57,19)){
#    n <- i+j
#    num <- c(num,n)
#  }
#  pdf(plot_name[count],height = 5,width = 8)
#  grid.arrange(p[[num[1]]],p[[num[2]]],p[[num[3]]],nrow=3)
#  dev.off()
#}
pdf("Temp.pdf",height = 20, width = 20)
grid.arrange(p[[1]],p[[20]],p[[39]],p[[2]],p[[21]],p[[40]],p[[3]],p[[22]],p[[41]],p[[4]],p[[23]],p[[42]],p[[5]],p[[24]],p[[43]],p[[6]],p[[25]],p[[44]],p[[7]],p[[26]],p[[45]],p[[8]],p[[27]],p[[46]],p[[9]],p[[28]],p[[47]],p[[10]],p[[29]],p[[48]],p[[11]],p[[30]],p[[49]],nrow=11)
dev.off()

pdf("Prec.pdf",height = 20, width = 20)
grid.arrange(p[[12]],p[[31]],p[[50]],p[[13]],p[[32]],p[[51]],p[[14]],p[[33]],p[[52]],p[[15]],p[[34]],p[[53]],p[[16]],p[[35]],p[[54]],p[[17]],p[[36]],p[[55]],p[[18]],p[[37]],p[[56]],p[[19]],p[[38]],p[[57]],nrow=8)
dev.off()

#只计算xp-clr top5的bayenv值
#/data1/home/yafei/003_Project3/bayenv/run1/XP-CLR/Top5/Alineage.top5.bed2，Blineage.top5.bed2，Dlineage.top5.bed2
#/data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225/Lineages/Alineage.E6_Landrace_locate.vcf.gz, Blineage.E6_Landrace_locate.vcf.gz, Dlineage.E6_Landrace_locate.vcf.gz
#/data2/yafei/003_Project3/bayenv
#xuebo:204:/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225/Lineages
#/data1/home/yafei/010_DataSet/gene_v1.1_withupdown300k_Lulab.bed

bedtools intersect -a Alineage.E6_Landrace_locate.vcf.gz -b Alineage.top5.bed2 -header | bgzip -c > Alineage.top5.vcf.gz &
bedtools intersect -a Blineage.E6_Landrace_locate.vcf.gz -b Blineage.top5.bed2 -header | bgzip -c > Blineage.top5.vcf.gz &
bedtools intersect -a Dlineage.E6_Landrace_locate.vcf.gz -b Dlineage.top5.bed2 -header | bgzip -c > Dlineage.top5.vcf.gz &

#A:2705021 B:2161251 D:2498868
#过基因及上下游300k
  bedtools intersect -a Alineage.top5.vcf.gz -b /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.bed -header | bgzip -c > Alineage.gene.top5.vcf.gz &
  bedtools intersect -a Blineage.top5.vcf.gz -b /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.bed -header | bgzip -c > Blineage.gene.top5.vcf.gz &
  bedtools intersect -a Dlineage.top5.vcf.gz -b /data1/home/yafei/010_DataSet/gene_v1.1_Lulab.bed -header | bgzip -c > Dlineage.gene.top5.vcf.gz &
#过LD(无)
  plink --vcf Alineage.gene5k.top5.vcf.gz --indep-pairwise 50 10 0.8 --out Alineage.filterLD --allow-extra-chr --autosome-num 42 &
  plink --vcf Alineage.gene5k.top5.vcf.gz --make-bed --extract Alineage.filterLD.prune.in --recode vcf-iid --out Alineage.filter.prune.in --autosome-num 42 &
  plink --vcf Blineage.gene5k.top5.vcf.gz --indep-pairwise 50 10 0.8 --out Blineage.filterLD --allow-extra-chr --autosome-num 42 &
  plink --vcf Blineage.gene5k.top5.vcf.gz --make-bed --extract Blineage.filterLD.prune.in --recode vcf-iid --out Blineage.filter.prune.in --autosome-num 42 &
  plink --vcf Dlineage.gene5k.top5.vcf.gz --indep-pairwise 50 10 0.8 --out Dlineage.filterLD --allow-extra-chr --autosome-num 42 &
  plink --vcf Dlineage.gene5k.top5.vcf.gz --make-bed --extract Dlineage.filterLD.prune.in --recode vcf-iid --out Dlineage.filter.prune.in --autosome-num 42 &
  
#/data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/V1/Gene50k/spid
  for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
java -jar -Xmx500g -Xms100g /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile ../Alineage.filter.prune.in.vcf -inputformat VCF -outputfile envgenofile/chr${chr}.envgenofile -outputformat BAYENV -spid /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/V1/Gene50k/spid/VCF_BAYENV_chr${chr}.spid
done
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
java -jar -Xmx500g -Xms100g /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile ../Blineage.filter.prune.in.vcf -inputformat VCF -outputfile envgenofile/chr${chr}.envgenofile -outputformat BAYENV -spid /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/V1/Gene50k/spid/VCF_BAYENV_chr${chr}.spid
done
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
java -jar -Xmx500g -Xms100g /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile ../Dlineage.filter.prune.in.vcf -inputformat VCF -outputfile envgenofile/chr${chr}.envgenofile -outputformat BAYENV -spid /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/V1/Gene50k/spid/VCF_BAYENV_chr${chr}.spid
done

#/data2/yafei/003_Project3/bayenv/spid/envgenofile/chr1.envgenofile 
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
./calc_bf.sh /data2/yafei/003_Project3/bayenv/spid/envgenofile/chr${chr}/chr${chr}.envgenofile /data2/yafei/003_Project3/bayenv/13pop.env /data2/yafei/003_Project3/bayenv/matrix/matrixA 13 100000 22 run1/chr${chr}_out
done

for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
./calc_bf.sh /data2/yafei/003_Project3/bayenv/spid/envgenofile/chr${chr}/chr${chr}.envgenofile /data2/yafei/003_Project3/bayenv/13pop.env /data2/yafei/003_Project3/bayenv/matrix/matrixB 13 100000 22 run1/chr${chr}_out
done

for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
./calc_bf.sh /data2/yafei/003_Project3/bayenv/spid/envgenofile/chr${chr}/chr${chr}.envgenofile /data2/yafei/003_Project3/bayenv/13pop.env /data2/yafei/003_Project3/bayenv/matrix/matrixD 13 100000 22 run1/chr${chr}_out
done
  