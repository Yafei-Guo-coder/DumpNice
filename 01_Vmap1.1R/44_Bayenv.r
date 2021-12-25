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
#-----------------------------------------结果处理-----------------------------------------
结果文件路径:
  筛选了LD的vcf文件：/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225/Landrace_225_noAM_220_maf001/lineage(Blineage_Landrace_225_noAM_220_maf001_LD.vcf)
  SNPsfileD：203@xuebo /data1/home/xuebo/Projects/Speciation/BAYENV/genofile/Dlineage2
  SNPsfileB：203@xuebo /data1/home/xuebo/Projects/Speciation/BAYENV/genofile/Blineage
  SNPsfileA：204@xuebo /data2/xuebo/Projects/Speciation/BAYENV/genofile/Alineage2
  D结果文件：203@xuebo /data1/home/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageD2
  B结果文件：203@xuebo /data1/home/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageB
  A结果文件：204@xuebo /data2/xuebo/Projects/Speciation/BAYENV/bayenv2_out_lineageA
结果处理:
  处理路径：204@yafei /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY
  把结果合并在A_out.bf，B_out.bf，D_out.bf文件中。
  提取vcf文件的位点(A:444085; B:498827; D:482174)，给bf文件进行位点注释(A:235606; B:276343; D:321195)。
total: 235606	276343	321195
Top001: 2356	2763	3211
Top005: 11780	13817	16059 
005Snp: 251773 303553 352583
uniqSnp: 30304 102996 120327

for i in {"A","B","D"}
do
cd /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/${i}
awk 'NR==FNR{a[$1]=$1;b[$1]=$0;c[$1]=$2"-"$3}NR!=FNR{print b[$1]"\t"c[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' ${i}.pos2.txt ${i}_out.bf > ${i}_out.bf2
bash ${i}get005.sh
cat *top5.txt | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n | uniq  > ${i}.bayenv.top5.bed
bedtools merge -d 50000 -i ${i}.bayenv.top5.bed |awk '{if($3-$2 != 1) print $0}' > ${i}.merge50k.bayenv.top5.bed
done

204@yafei /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY
30304.bayenv.top5.bed是Alineage的top5的snp位点。
3628.merge50k.bayenv.top5.bed是合并Alineage的距离小于50k的snp，并且把单个的snp去掉。
102996.bayenv.top5.bed是Blineage的top5的snp位点。
14470.merge50k.bayenv.top5.bed是合并Blineage的距离小于50k的snp，并且把单个的snp去掉。
120327.bayenv.top5.bed是Dlineage的top5的snp位点。
14359.merge50k.bayenv.top5.bed是合并Dlineage的距离小于50k的snp，并且把单个的snp去掉。

比对基因:
204@xuebo /data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%

----------------------------step1:提取XPCLR和bayenv重叠的XPCLR区域-------------------------
for i in `ls *_A.shuf.top5.bed`
do
bedtools intersect -a $i -b A.merge10k.bayenv.top5.bed -wa |sort | uniq > bayenv_over/$i
done
for i in `ls *_B.shuf.top5.bed`
do
bedtools intersect -a $i -b B.merge10k.bayenv.top5.bed -wa |sort | uniq> bayenv_over/$i
done
for i in `ls *_D.shuf.top5.bed`
do
bedtools intersect -a $i -b D.merge10k.bayenv.top5.bed -wa |sort | uniq> bayenv_over/$i
done

-----------------step3:提取所有克隆基因区上下游50k snp，重新做bayenv----------------
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
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' Alineage_out.bf > Alineage_out.bf2
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' Blineage_out.bf > Blineage_out.bf2
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' Dlineage_out.bf > Dlineage_out.bf2
zcat ../Alineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > A.pos2.txt
zcat ../Blineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > B.pos2.txt
zcat ../Dlineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > D.pos2.txt
for i in {"A","B","D"}
do
awk 'NR==FNR{a[$1]=$1;b[$1]=$0;c[$1]=$2"-"$3}NR!=FNR{print b[$1]"\t"c[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' ${i}.pos2.txt ${i}lineage_out.bf2 > ${i}lineage_out.bf2
cat /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/${i}/${i}_out.bf2 ${i}lineage_out.bf2 > ${i}_out.bf2
bash ${i}get005.sh
cat ${i}bio*top5.txt | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n | uniq  > ${i}.bayenv.top5.bed
bedtools merge -d 50000 -i ${i}.bayenv.top5.bed |awk '{if($3-$2 != 1) print $0}' > ${i}.merge50k.bayenv.top5.bed
done

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
结果合并解析
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' A3_out.bf > A3_out.bf2
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' B3_out.bf > B3_out.bf2
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' D3_out.bf > D3_out.bf2
zcat ../Alineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > A.pos2.txt
zcat ../Blineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > B.pos2.txt
zcat ../Dlineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > D.pos2.txt

for i in {"A","B","D"}
do
for j in {1..3}
do
awk 'NR==FNR{a[$1]=$1;b[$1]=$0;c[$1]=$2"-"$3}NR!=FNR{print b[$1]"\t"c[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' POS/${i}${j}.pos.txt ${i}${j}_out.bf > ${i}${j}_out.bf2
done
done

cat A1_out.bf2 A2_out.bf2 A3_out.bf2 > A_all.bf (329977)
cat B1_out.bf2 B2_out.bf2 B3_out.bf2 > B_all.bf (355505)
cat D1_out.bf2 D2_out.bf2 D3_out.bf2 > D_all.bf (401549)

5%: A: 16498 B: 17775 D: 20077

total: A: 354783 B: 385077 D: 435991
uniq: A: 38572 B: 99203 D: 97149

cat *top5.txt | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n | uniq  > ${i}.bayenv.top5.bed
bedtools merge -d 50000 -i ${i}.bayenv.top5.bed |awk '{if($3-$2 != 1) print $0}' > ${i}.merge50k.bayenv.top5.bed

--------------------计算富集程度-----------

sed '1d' EU_South_smooth_D.txt|shuf -n 4476 >shuf005/EU_South_smooth_D.top5.bed
sed '1d' North2_South_smooth_D.txt|shuf -n 5078 >shuf005/North2_South_smooth_D.top5.bed
sed '1d' Strang_WA_smooth_D.txt|shuf -n 5404 >shuf005/Strang_WA_smooth_D.top5.bed
sed '1d' Tibet_South_smooth_D.txt|shuf -n 6658 >shuf005/Tibet_South_smooth_D.top5.bed
sed '1d' WA_EU_smooth_D.txt|shuf -n 6204 >shuf005/WA_EU_smooth_D.top5.bed
sed '1d' WA_South_smooth_D.txt|shuf -n 4296 >shuf005/WA_South_smooth_D.top5.bed

sed '1d' EU_South_smooth_A.txt|shuf -n 5779 >shuf005/EU_South_smooth_A.top5.bed
sed '1d' North2_South_smooth_A.txt|shuf -n 6122 >shuf005/North2_South_smooth_A.top5.bed
sed '1d' Tibet_South_smooth_A.txt|shuf -n 6289 >shuf005/Tibet_South_smooth_A.top5.bed
sed '1d' WA_EU_smooth_A.txt|shuf -n 5795 >shuf005/WA_EU_smooth_A.top5.bed
sed '1d' WA_South_smooth_A.txt|shuf -n 5886 >shuf005/WA_South_smooth_A.top5.bed

sed '1d' EU_South_smooth_B.txt|shuf -n 7706 >shuf005/EU_South_smooth_B.top5.bed
sed '1d' North2_South_smooth_B.txt|shuf -n 7089 >shuf005/North2_South_smooth_B.top5.bed
sed '1d' Tibet_South_smooth_B.txt|shuf -n 6107 >shuf005/Tibet_South_smooth_B.top5.bed
sed '1d' WA_EU_smooth_B.txt|shuf -n 7447 >shuf005/WA_EU_smooth_B.top5.bed
sed '1d' WA_South_smooth_B.txt|shuf -n 6980 >shuf005/WA_South_smooth_B.top5.bed

for i in `ls *_A.top5.bed`
do
bedtools intersect -a $i -b Alog.merge50k.bayenv.shuf.1.5.bed -wa |sort | uniq > bayenv_over4/region_shuf_$i
done
for i in `ls *_B.top5.bed`
do
bedtools intersect -a $i -b Blog.merge50k.bayenv.shuf.1.5.bed -wa |sort | uniq> bayenv_over4/region_shuf_$i
done
for i in `ls *_D.top5.bed`
do
bedtools intersect -a $i -b Dlog.merge50k.bayenv.shuf.1.5.bed -wa |sort | uniq> bayenv_over4/region_shuf_$i
done

shuf -n 38572 A_all.bf | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n > A.bayenv.shuf2.top5.bed
bedtools merge -d 50000 -i A.bayenv.shuf2.top5.bed |awk '{if($3-$2 != 1) print $0}' > A.merge50k.bayenv.shuf2.top5.bed

shuf -n 99203 B_all.bf | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n > B.bayenv.shuf2.top5.bed
bedtools merge -d 50000 -i B.bayenv.shuf2.top5.bed |awk '{if($3-$2 != 1) print $0}' > B.merge50k.bayenv.shuf2.top5.bed

shuf -n 97149 D_all.bf | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n > D.bayenv.shuf2.top5.bed
bedtools merge -d 50000 -i D.bayenv.shuf2.top5.bed |awk '{if($3-$2 != 1) print $0}' > D.merge50k.bayenv.shuf2.top5.bed

--------------------提取log-----
bash get005.sh 
cat Alog_bio*txt | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n | uniq  > Alog.bayenv.1.5.bed
cat Blog_bio*txt | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n | uniq  > Blog.bayenv.1.5.bed
cat Dlog_bio*txt | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n | uniq  > Dlog.bayenv.1.5.bed

bedtools merge -d 50000 -i Alog.bayenv.1.5.bed |awk '{if($3-$2 != 1) print $0}' > Alog.merge50k.bayenv.1.5.bed
bedtools merge -d 50000 -i Blog.bayenv.1.5.bed |awk '{if($3-$2 != 1) print $0}' > Blog.merge50k.bayenv.1.5.bed
bedtools merge -d 50000 -i Dlog.bayenv.1.5.bed |awk '{if($3-$2 != 1) print $0}' > Dlog.merge50k.bayenv.1.5.bed

scp *1.5.bed xuebo@159.226.116.204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%

wc -l *1.5.bed
 
shuf -n 62529 A_all.bf | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n > Alog.bayenv.shuf.1.5.bed
shuf -n 16005 B_all.bf | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n  > Blog.bayenv.shuf.1.5.bed
shuf -n 19991 D_all.bf | awk '{print $2"\t"$3-1"\t"$3}'| sort -k1,1n -k2,2n  > Dlog.bayenv.shuf.1.5.bed
bedtools merge -d 50000 -i Alog.bayenv.shuf.1.5.bed |awk '{if($3-$2 != 1) print $0}' > Alog.merge50k.bayenv.shuf.1.5.bed
bedtools merge -d 50000 -i Blog.bayenv.shuf.1.5.bed |awk '{if($3-$2 != 1) print $0}' > Blog.merge50k.bayenv.shuf.1.5.bed
bedtools merge -d 50000 -i Dlog.bayenv.shuf.1.5.bed |awk '{if($3-$2 != 1) print $0}' > Dlog.merge50k.bayenv.shuf.1.5.bed

scp *1.5.bed xuebo@159.226.116.204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%


for i in `ls *smooth_A.top1.bed`
do
bedtools intersect  -a Alog.bayenv.1.5.bed -b  $i -wa |sort | uniq > bayenv_over/$i
done
for i in `ls *smooth_B.top1.bed`
do
bedtools intersect  -a Blog.bayenv.1.5.bed -b  $i -wa |sort | uniq > bayenv_over/$i
done
for i in `ls *smooth_D.top1.bed`
do
bedtools intersect -a Dlog.bayenv.1.5.bed -b $i -wa |sort | uniq > bayenv_over/$i
done

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

180893 A_out.bf uniq -> 180882
204653 B_out.bf uniq -> 202982
197614 D_out.bf uniq -> 197523

444085 A_maf001_LD.vcf
498827 B_maf001_LD.vcf
482174 D_maf001_LD.vcf

awk '{print $1"\t"$2-1"\t"$2"\t"NR}' A_maf001_LD.txt > A_maf001_LD.bed
awk '{print $1"\t"$2-1"\t"$2"\t"NR}' B_maf001_LD.txt > B_maf001_LD.bed
awk '{print $1"\t"$2-1"\t"$2"\t"NR}' D_maf001_LD.txt > D_maf001_LD.bed

awk 'NR==FNR{a[$4]=$1;b[$4]=$2;c[$4]=$3;d[$4]=$4}NR!=FNR{if($1 in d) print a[$1],b[$1],c[$1],$0}' A.pos.bed A_out.bf > A_all.bf
awk 'NR==FNR{a[$4]=$1;b[$4]=$2;c[$4]=$3;d[$4]=$4}NR!=FNR{if($1 in d) print a[$1],b[$1],c[$1],$0}' B.pos.bed B_out.bf > B_all.bf
awk 'NR==FNR{a[$4]=$1;b[$4]=$2;c[$4]=$3;d[$4]=$4}NR!=FNR{if($1 in d) print a[$1],b[$1],c[$1],$0}' D.pos.bed D_out.bf > D_all.bf

AB 密度 1.905844e-06
bash get005.sh
cat Abio*.top5.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n | uniq  > Alog.bayenv.005.bed
cat Bbio*.top5.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n | uniq  > Blog.bayenv.005.bed
cat Dbio*.top5.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n | uniq  > Dlog.bayenv.005.bed

180893 A_all.bf 0.05 -> 9045 密度 -> 14199/4934891648 -> 1.832867e-06 2.877267e-06
204653 B_all.bf 0.05 -> 10233 密度 -> 45868/5180314468 -> 1.975363e-06 8.854289e-06
197614 D_all.bf 0.05 -> 9881 密度 -> 50085/3951074735 ->2.500839e-06 1.26763e-05

135848 A_all.bf 0.05 -> 6792
98705 B_all.bf 0.05 -> 4935

11498 Alog.bayenv.005.bed 背景密度 -> 11498/4934891648 -> 2.32994e-06
24221 Blog.bayenv.005.bed 背景密度 -> 24221/5180314468 -> 4.675585e-06

-----------------------------------------20211224 新的结果--------------------------------------------
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/bayenv2_out_lineageA
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/bayenv2_out_lineageB
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/genofile/Alineage_bayenv_pop5_500K.vcf
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/genofile/Blineage_bayenv_pop5_500K.vcf
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/genofile/Alineage.pop5.envgenofile
204@xuebo:/data2/xuebo/Projects/Speciation/BAYENV/BAYENV_V3/genofile/Blineage.pop5.envgenofile

---------------------------bayenv_xpclr_overlap 计算富集程度----------------------
204@yafei:/data1/home/yafei/003_Project3/bayenv/V3_500k
for i in `ls *smooth_A.top5.bed`
do
bedtools intersect -a Alog.bayenv.005.bed -b  $i -wa |sort | uniq > bayenv_xpclr2/$i
done
for i in `ls *smooth_B.top5.bed`
do
bedtools intersect -a Blog.bayenv.005.bed -b  $i -wa |sort | uniq > bayenv_xpclr2/$i
done
for i in `ls *smooth_D.top5.bed`
do
bedtools intersect -a Dlog.bayenv.005.bed -b $i -wa |sort | uniq > bayenv_xpclr/$i
done

#--------------------------------------定位基因-----------------------------------
204@yafei:/data1/home/yafei/003_Project3/bayenv/V3_500k/gene
for i in `ls *smooth_A.top5.bed`
do
bedtools intersect -a $i -b Alog.bayenv.005.bed -wa |sort | uniq > gene/$i
done
for i in `ls *smooth_B.top5.bed`
do
bedtools intersect -a $i -b Blog.bayenv.005.bed -wa |sort | uniq > gene/$i
done
for i in `ls *smooth_D.top5.bed`
do
bedtools intersect -a $i -b Dlog.bayenv.005.bed -wa |sort | uniq > gene/$i
done

定位已注释基因:gff
for i in `ls *bed`
do
bedtools intersect -a gene_v1.1_Lulab.gff3 -b $i -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > ${i::-3}gff.gene
done
#定位已克隆基因:cloned gene 以及定位抗病基因:nlr gene
for i in `ls *gff.gene`
do
awk -F"ID=" '{print $2}' $i > ${i::-8}gene2
done
for i in `ls *.gene2`
do
grep -w -f $i all_cloned_gene.txt > ${i::-5}cloned.gene
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
v <- c("A","B")
for ( i in c(1:2)){
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

library(reshape)
name <- read.table("file_prefix.txt",header=F,stringsAsFactors=F)
colnames(name) <- "file"
v <- c("A","B","D")
num <- c(1,0.6,0.2)
for ( i in c(1:3)){
  filename <- paste(v[i],"_heatmap_format1.txt",sep="")
  file <- read.table(filename,header=F,stringsAsFactors=F)
  sub <- file[,c(1,3,4)]
  colnames(sub) <- c("file","gene","type")
  sub[which(sub$type==1),3] <- num[i]
  sub[which(sub$type==0),3] <- num[i]-0.2
  mt <- melt(sub,id=c("file","gene"))
  st <- cast(mt,file~variable+gene)
  name<-merge(name,st,by="file") 
}
write.table(name,"heatmap_format2.txt",quote=F,row.names=F,col.names=T,sep="\t")

#-----------------------------------bayenv_overlap--------------------------------
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv")
A <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv/21.Alog.bayenv.005.uniq.txt",header=T,stringsAsFactors = F)
B <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv/21.Blog.bayenv.005.uniq.txt",header=T,stringsAsFactors = F)
D <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/04_bayenv/21.Dlog.bayenv.005.uniq.txt",header=T,stringsAsFactors = F)
ggplot(D, aes(x = Value)) +
  geom_histogram(binwidth = 1, fill = "lightblue", colour = "black")+
  theme_classic()

require(RIdeogram)
density <- rbind(A,B,D)
wheat_karyotype <- read.table("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/13_Plots/01_Density/wheat_karyotype.txt", header=T, stringsAsFactors = F)
ideogram(karyotype = wheat_karyotype, overlaid = density)
convertSVG("chromosome.svg", device = "pdf")

-----------------step3:提取所有克隆基因区上下游500k snp，重新做bayenv----------------
Working directory: 
/data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/V1/Gene50k


yafei@204: /data1/home/yafei/003_Project3/bayenv/Gene_50k
bedtools merge -i <(awk '{print $1"\t"$2-50000"\t"$3+50000}' clone.gene.txt) > gene.50k.bed

1.区分A,B,D亚基因组
awk '{output="chr"$1".txt"; print $0 > output}' gene.50k.bed2 

2.从225个样本的vcf文件中提取出这些位点的snp
#yafei@204:/data1/home/yafei/003_Project3/Structure/E6_Landrace_locate_225/Lineage
bcftools view -R ../pos/A.pos.bed Alineage_145.vcf.gz -o A.gene.500k.vcf.gz -O z &
bcftools view -R ../pos/B.pos.bed Blineage_145.vcf.gz -o B.gene.500k.vcf.gz -O z &
bcftools view -R ../pos/D.pos.bed Dlineage_145.vcf.gz -o D.gene.500k.vcf.gz -O z &

3.格式转换
#!/bin/bash
for chr in {"A","B","D"}
do
java -jar -Xmx500g -Xms100g /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile Chr/chr${chr}.gene.50k.vcf -inputformat VCF -outputfile envgenofile/chr${chr}.envgenofile -outputformat BAYENV -spid spid/VCF_BAYENV_chr${chr}.spid
done
4.跑bayenv
./calc_bf.sh lineage/Alineage.envgenofile 13pop.env matrix/A.matrix 13 100000 22 Alineage_out
./calc_bf.sh lineage/Blineage.envgenofile 13pop.env matrix/B.matrix 13 100000 22 Blineage_out
./calc_bf.sh lineage/Dlineage.envgenofile 13pop.env matrix/D.matrix 13 100000 22 Dlineage_out
5.结果解析
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' Alineage_out.bf > Alineage_out.bf2
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' Blineage_out.bf > Blineage_out.bf2
awk '{print NR"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' Dlineage_out.bf > Dlineage_out.bf2
zcat ../Alineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > A.pos2.txt
zcat ../Blineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > B.pos2.txt
zcat ../Dlineage.gene.50k.vcf.gz | awk '{print $1"\t"$2}' |grep -v "#" | awk '{print NR"\t"$1"\t"$2}' > D.pos2.txt
for i in {"A","B","D"}
do
awk 'NR==FNR{a[$1]=$1;b[$1]=$0;c[$1]=$2"-"$3}NR!=FNR{print b[$1]"\t"c[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' ${i}.pos2.txt ${i}lineage_out.bf2 > ${i}lineage_out.bf2
cat /data1/home/yafei/003_Project3/Structure/bayenv/ENVBAY/${i}/${i}_out.bf2 ${i}lineage_out.bf2 > ${i}_out.bf2
bash ${i}get005.sh
cat ${i}bio*top5.txt | awk '{print $2"\t"$3-1"\t"$3}' | sort -k1,1n -k2,2n | uniq  > ${i}.bayenv.top5.bed
bedtools merge -d 50000 -i ${i}.bayenv.top5.bed |awk '{if($3-$2 != 1) print $0}' > ${i}.merge50k.bayenv.top5.bed
done



