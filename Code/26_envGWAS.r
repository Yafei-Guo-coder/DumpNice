#444126 Alineage_Landrace_225_noAM_220_maf001_LD.vcf
#498868 Blineage_Landrace_225_noAM_220_maf001_LD.vcf
#482215 Dlineage_Landrace_225_noAM_220_maf001_LD.vcf

#工作目录
#服务器:/data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_220/maf005/env1/Result
#本地:/Users/guoyafei/Documents/01_Migration/02_Environment/09_GWAS

library(dplyr)
library(ggridges)
library(RColorBrewer)
library(ggrepel)
require(gridExtra)
env1 和matrix1是单个亚基因组做群体结构，分开做GWAS
env2 和matrix2是三个亚基因组合起来做一个群体结构，做GWAS

#准备环境变量文件
for i in {2..23}
do
a=`expr $i - 1`
for j in `cat 220.txt`
do
grep -w $j 225env.txt
done  | awk '{print $1"\t"$'"$i"'}' | sed "1i <Trait>\tenv$a" > env${a}.txt
done

1.按MAF>0.05和缺失率<0.1过滤
plink --vcf Alineage_Landrace_225_noAM_220_maf001_LD.vcf --maf 0.05 --geno 0.1 --recode vcf-iid --out Alineage_Landrace_220_maf005 --allow-extra-chr --double-id --autosome-num 42 &
plink --vcf Blineage_Landrace_225_noAM_220_maf001_LD.vcf --maf 0.05 --geno 0.1 --recode vcf-iid --out Blineage_Landrace_220_maf005 --allow-extra-chr --double-id --autosome-num 42 &
plink --vcf Dlineage_Landrace_225_noAM_220_maf001_LD.vcf --maf 0.05 --geno 0.1 --recode vcf-iid --out Dlineage_Landrace_220_maf005 --allow-extra-chr --double-id --autosome-num 42 &

plink --vcf Alineage_Landrace_220_maf005.vcf --recode 12 --out Alineage_Landrace_220_maf005 --autosome-num 42 &
plink --vcf Blineage_Landrace_220_maf005.vcf --recode 12 --out Blineage_Landrace_220_maf005 --autosome-num 42 &
plink --vcf Dlineage_Landrace_220_maf005.vcf --recode 12 --out Dlineage_Landrace_220_maf005 --autosome-num 42 &

2.计算structure
admixture --cv Alineage_Landrace_220_maf005.ped 5 >> Alog.txt &
admixture --cv Blineage_Landrace_220_maf005.ped 5 >> Blog.txt &
admixture --cv Dlineage_Landrace_220_maf005.ped 5 >> Dlog.txt &
grep "CV error" log.txt > k_1to13

A:10 B:8 D:8
PCA 分析：
plink --vcf Alineage_Landrace_220_maf005.vcf --pca 3 header tabs -out Alineage.maf0.05 --double-id --autosome-num 42 &
plink --vcf Blineage_Landrace_220_maf005.vcf --pca 3 header tabs -out Blineage.maf0.05 --double-id --autosome-num 42 &
plink --vcf Dlineage_Landrace_220_maf005.vcf --pca 3 header tabs -out Dlineage.maf0.05 --double-id --autosome-num 42 &

3.亲缘关系矩阵
run_pipeline.pl -Xmx200g -importGuess Alineage_Landrace_220_maf005.vcf -KinshipPlugin -method Centered_IBS -endPlugin -export A_kinship.txt -exportType SqrMatrix &
run_pipeline.pl -Xmx200g -importGuess Blineage_Landrace_220_maf005.vcf -KinshipPlugin -method Centered_IBS -endPlugin -export B_kinship.txt -exportType SqrMatrix &
run_pipeline.pl -Xmx200g -importGuess Dlineage_Landrace_220_maf005.vcf -KinshipPlugin -method Centered_IBS -endPlugin -export D_kinship.txt -exportType SqrMatrix &
  run_pipeline.pl -Xmx200g -importGuess ABDlineage_Landrace_220_maf005.vcf -KinshipPlugin -method Centered_IBS -endPlugin -export ABD_kinship.txt -exportType SqrMatrix &
-mlmVarCompEst P3D

4.VCF转格式
run_pipeline.pl -Xmx200g -fork1 -vcf ../Alineage_Landrace_220_maf005.vcf -export Alineage_maf005 -exportType Hapmap -runfork1 &
run_pipeline.pl -Xmx200g -fork1 -vcf ../Blineage_Landrace_220_maf005.vcf -export Blineage_maf005 -exportType Hapmap -runfork1 &
run_pipeline.pl -Xmx200g -fork1 -vcf ../Dlineage_Landrace_220_maf005.vcf -export Dlineage_maf005 -exportType Hapmap -runfork1 &
  
run_pipeline.pl -Xmx200g -SortGenotypeFilePlugin -inputFile Alineage_maf005.hmp.txt  -outputFile Alineage_maf005_sort -fileType Hapmap &
run_pipeline.pl -Xmx200g -SortGenotypeFilePlugin -inputFile Blineage_maf005.hmp.txt  -outputFile Blineage_maf005_sort -fileType Hapmap &
run_pipeline.pl -Xmx200g -SortGenotypeFilePlugin -inputFile Dlineage_maf005.hmp.txt  -outputFile Dlineage_maf005_sort -fileType Hapmap &
  
thread_num=15
tempfifo="my_temp_fifo"
mkfifo ${tempfifo}
exec 6<>${tempfifo}
rm -f ${tempfifo}
for ((i=1;i<=${thread_num};i++))
  do
{
  echo
}
done >&6
for i in {8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7}
do
{
  read -u6
  {
    run_pipeline.pl -fork1 -h ABDlineage_maf005.hmp.txt -fork2 -r env${i}.txt -fork3 -q ABD.Q8.matrix -excludeLastTrait -fork4 -k ABD_kinship.txt -combine5  -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export env${i}_A_mlm -runfork1 -runfork2 -runfork3
    echo "" >&6
  } &
}
done
wait
exec 6>&-
  
for i in `ls *mlm2.txt`; do sed '1d' $i | grep -v "NaN" | awk '{print $2"\t"$3"\t"$4"\t"$7}' | sed '1i SNP\tCHR\tBP\tP'> Result/${i::-9}.txt ; done
for i in *txt
do
sed '1d' $i | shuf -n 5000 | sed '1i SNP\tCHR\tBP\tP' > ${i::-4}.qq.txt
done

5.本地画图
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/09_GWAS/LD01/QQ")
#QQ plot-----
#pdf("qq2.pdf")
v <- c("AB","D")
for ( i in c(1:20)){
  p <- list()
  name1 <- paste("env",i,".pdf",sep="")
  for ( chr in c(1:2)){
    name <- paste("env",i,": ",v[chr],sep="")
    filename <- paste("env",i,"_",v[chr],".qq.txt",sep="")
    data <- read.table(filename,header=T,stringsAsFactors = F)
    colnames(data) <- c("SNP", "CHR", "BP","P")
    qq_dat <- data.frame(obs=-log10(sort(data$P,decreasing=FALSE)),
                         exp=-log10( ppoints(length(data$P))))
    p[[chr]] <- ggplot(data=qq_dat,aes(exp,obs))+
      geom_point(alpha=0.7,color="#7F7F7FFF")+
      geom_abline(color="#D62728FF")+
      xlab("Expected -log10(P-value)")+
      ylab("Observed -log10(P-value)")+
      scale_x_continuous(limits = c(0,7))+
      scale_y_continuous(limits = c(0,7))+
      #ggtitle(name)+
      theme(
        plot.title = element_text(color="red", size=20, face="bold.italic"),
        axis.title = element_text(size=12,face="bold"),
        axis.text = element_text(face="bold",size=8,color = "black"),
        #axis.line = element_line(size=0.8,color="black"),
        axis.ticks= element_line(size=0.8,colour = "black"),
        panel.grid =element_blank(),
        panel.border = element_rect(fill=NA,size = 0.8),
        panel.background = element_blank())
  }
  pdf(name1,width = 6,height = 3)
  grid.arrange(p[[1]],p[[2]],nrow=1)
  dev.off()
}

#--------------manhuttun-------------
library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/09_GWAS/LD01")
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/09_GWAS/LD01/Manhuttan" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F)})

#threshold
library(gdata)
thresh <- read.xls("thresh.xlsx",sheet=1,row.name=1,na.strings=c("NA","#DIV/0!"))

#for (i in c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58)){
for (i in seq(1,40,by=2)){
  filename <- strsplit( names(data)[i], "_")[[1]][1]
  p<-list()
  for (j in c(0,1)){
    a <- i+j
    gwasResults2 <- data[[a]]
    colnames(gwasResults2) <- c("SNP", "CHR", "BP","P")
    gwasResults3 <- gwasResults2[which(gwasResults2$P < 0.01),]
    other <- gwasResults2[which(gwasResults2$P > 0.01 & gwasResults2$P <0.2),]
    #other2 <- other[sample(nrow(other), 20000), ]
    gwasResults <- rbind(gwasResults3,other)
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
    p[[j+1]] <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
      # Show all points
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      # Add highlighted points
      #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
      # Add label using ggrepel to avoid overlapping
      #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size = 15),
        axis.title.x=element_text(size = 15),
      )+
      scale_y_continuous(limits = c(0,7))+
      geom_point(data=point,aes(x=BPcum,y=-log10(P)),color="red")
      #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
      #geom_hline(yintercept = -log10(thresh[a,1]), colour="red",linetype=2, size=1)+
      #geom_hline(yintercept = -log10(thresh[a,2]), colour="blue",linetype=2, size=1)
  }
  pdf(paste(filename,".pdf",sep=""),height = 9,width = 9)
  grid.arrange(p[[1]],p[[2]],nrow=2)
  dev.off()
}
----------------跑20个随机shuf的情况-----------
for j in {51..100}
do
mkdir run${j}
for i in env*txt
do 
sed '1d' $i | shuf | sed "1i <Trait>\t${i::-4}" | awk '{print $2}' > test
paste names test -d "\t" > run${j}/$i
done
done

for j in {51..100}
do
cp gwas.sh run${j}
done

-----------------整理shuf的结果------------------
{1,12,13,14,15,16,17,18,19,20,21,22}

for i in {2,3,4,5,6,7,8,9,10,11}
do
cd run${i}
mkdir Result
for j in `ls *mlm2.txt`; do sed '1d' $j | grep -v "NaN" | awk '{print $2"\t"$3"\t"$4"\t"$7}' | sed '1i SNP\tCHR\tBP\tP'> Result/${j::-9}.txt ; done
cd Result
for m in *txt; do sed '1d' $m | shuf -n 5000 | sed '1i SNP\tCHR\tBP\tP'  > ${m::-4}.qq.txt; done
done

for i in {2,3,4,5,6,7,8,9,10,11}
do
{
  cd /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_220/maf005/PCA/shufTrait/run${i}/Result/
    for j in {"A","B","D"}
  do
  {
    for m in *${j}.txt; do sort -k4,4g $m | head -n 2 | tail -n 1 ; echo $i; done 
  }
  done
}
done


for i in {23..50}
do
cd run${i}
mkdir Result
for j in `ls *mlm2.txt`; do sed '1d' $j | grep -v "NaN" | awk '{print $2"\t"$3"\t"$4"\t"$7}' | sed '1i SNP\tCHR\tBP\tP'> Result/${j::-9}.txt ; done
cd Result
for m in *txt; do sed '1d' $m | shuf -n 5000 | sed '1i SNP\tCHR\tBP\tP'  > ${m::-4}.qq.txt; done
done


for i in {23..50}
do
{
  cd /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_220/maf005/PCA/shufTrait/run${i}/Result/
    for j in {"A","B","D"}
  do
  {
    for m in *${j}.txt; do sort -k4,4g $m | head -n 2 | tail -n 1 ; echo $i; done 
  }
  done
}
done

for i in {1,12,13,14,15,16,17,18,19,20,21,22}
do
cd run${i}
mkdir Result
for j in `ls *mlm2.txt`; do sed '1d' $j | grep -v "NaN" | awk '{print $2"\t"$3"\t"$4"\t"$7}' | sed '1i SNP\tCHR\tBP\tP'> Result/${j::-9}.txt ; done
cd Result
for m in *txt; do sed '1d' $m | shuf -n 5000 | sed '1i SNP\tCHR\tBP\tP'  > ${m::-4}.qq.txt; done
done

#thresh005
for i in {1..50}
do
{
  cd /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_220/maf005/PCA/shufTrait/run${i}/Result/
  for m in *A.txt; do sort -k4,4g $m | head -n 7519 | tail -n 1 | awk '{print $4}'; echo $m; done 
  for n in *B.txt; do sort -k4,4g $n | head -n 7381 | tail -n 1 | awk '{print $4}'; echo $n; done 
  for q in *D.txt; do sort -k4,4g $q | head -n 5720 | tail -n 1 | awk '{print $4}'; echo $q; done 
}
done
#thresh001
for i in {1..50}
do
{
  cd /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_220/maf005/PCA/shufTrait/run${i}/Result/
  for m in *A.txt; do sort -k4,4g $m | head -n 1504 | tail -n 1 | awk '{print $4}'; echo $m; done 
  for n in *B.txt; do sort -k4,4g $n | head -n 1476 | tail -n 1 | awk '{print $4}'; echo $n; done 
  for q in *D.txt; do sort -k4,4g $q | head -n 1144 | tail -n 1 | awk '{print $4}'; echo $q; done 
}
done

setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/shuf")
pdf("density.pdf")
v <- c("A","B","D")
for ( i in c(1:20)){
  for ( chr in c(1:3)){
    name <- paste("env",i,": ",v[chr],sep="")
    filename <- paste("env",i,"_",v[chr],".qq.txt",sep="")
    data <- read.table(filename,header=T,stringsAsFactors = F)
    data<- read.table("env2_B.qq.txt", header=T,stringsAsFactors = F)
    a <- as.numeric(quantile(data$P,0.001))
    b <- as.numeric(quantile(data$P,0.01))
    p <- ggplot(data, aes(x=P)) + 
      geom_density()+
      scale_x_log10()+
      geom_vline(xintercept = a, color = 'gray', size = 0.5) + 
      geom_vline(xintercept = b, color = 'skyblue', size = 0.5) 
    print(pd_qq)
  }
}
dev.off()

#合并去掉单个位置（位点富集以及区域富集）

for i in {"A","B","D"}
do
for j in {1..20}
do
bedtools merge -i Top0001.env${j}_${i}.txt -d 1000000 -c 1 -o count | awk '{if($4>1) print $0}' > Top0001.env${j}_${i}.region
bedtools intersect -a Top0001.env${j}_${i}.txt -b Top0001.env${j}_${i}.region -wa > Top0001.env${j}_${i}.site
done
cat Top0001.env${j}_${i}.site | sort | uniq > Merge.pos/${i}.site.bed
done




