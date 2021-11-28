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
#204:yafei:/data1/home/yafei/003_Project3/Structure/E6_Landrace_locate_225/bayenv
#准备环境变量文件
awk 'NR==FNR{a[$1]=$2;b[$1]=$1}NR!=FNR{if($1 in b) print $0,a[$1]}' pop.txt 225env.txt | sort -k24,24 > merge_env.txt
datamash groupby 24 mean 2 mean 3 mean 4 mean 5 mean 6 mean 7 mean 8 mean 9 mean 10 mean 11 mean 12 mean 13 mean 14 mean 15 mean 16 mean 17 mean 18 mean 19 mean 20 mean 21 mean 22 mean 23 < merge_env.txt | datamash transpose > format1.txt
#R
data <- read.table("format1.txt",header=T,stringsAsFactors=F)
m <- apply(data,1,mean)
s <- apply(data,1,sd)
sub <- (data-m)/s
write.table(sub,"format2.txt", sep="\t", quote=F,row.names=F)
#shell
datamash transpose < format2.txt > format3.txt

#准备基因型文件
vcftools --vcf chr36.E6_Landrace_locate.vcf --maf 0.01 --recode --recode-INFO-all --out chr36.maf0.01
java -Xmx200g -Xms512m -jar /data1/home/yafei/008_Software/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile chr36.maf0.01.recode.vcf -inputformat VCF -outputfile chr36.maf0.01.env -outputformat BAYENV -spid VCF_BAYENV.spid

#运行bayenv
#matrix estimation
bayenv2 -i hgdp_no_X_37_freqs -s samplesize.txt -p 52 -k 1000 -r 83556 -o out_matrix
bayenv2 -i rs316 -m hgdp_matrix_1 -e PCs.env -p 52 -k 1000 -n 4 -t -r 42 -o out_correlation
calc_bf.sh hgdp_no_X_37_freqs PCs.env hgdp_matrix_1 52 1000 4

MATRIXFILE 单个协方差矩阵，矩阵中的条目由制表符分隔。 这可以通过从程序的输出中复制和粘贴最后打印的矩阵来获得。 协方差矩阵示例如下： hgdp 矩阵 1. 如果您担心协方差矩阵绘制的变化，那么您可以考虑在此处使用迭代中的平均协方差矩阵作为输入


