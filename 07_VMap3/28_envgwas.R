#environmental gwas
#203:/data1/home/yafei/004_Vmap3/gwas
#203:/data1/home/yafei/004_Vmap3/gwas/genotype

1.按MAF>0.05和缺失率<0.1过滤
for i in {001..042}
do
plink --bfile /data4/home/yafei/plink_VCF/maf005r202/chr${i}.filter.prune.in --keep group.txt --maf 0.05 --geno 0.1 --recode vcf-iid --out chr${i}_GWAS --allow-extra-chr --double-id --autosome-num 42 &
done

2.合并VCF
vcf-concat chr001_GWAS.vcf chr002_GWAS.vcf chr003_GWAS.vcf chr004_GWAS.vcf chr007_GWAS.vcf chr008_GWAS.vcf chr009_GWAS.vcf chr010_GWAS.vcf chr013_GWAS.vcf chr014_GWAS.vcf chr015_GWAS.vcf chr016_GWAS.vcf chr019_GWAS.vcf chr020_GWAS.vcf chr021_GWAS.vcf chr022_GWAS.vcf chr025_GWAS.vcf chr026_GWAS.vcf chr027_GWAS.vcf chr028_GWAS.vcf chr031_GWAS.vcf chr032_GWAS.vcf chr033_GWAS.vcf chr034_GWAS.vcf chr037_GWAS.vcf chr038_GWAS.vcf chr039_GWAS.vcf chr040_GWAS.vcf | bgzip -c > ABlineage_gwas.vcf.gz 
vcf-concat chr005_GWAS.vcf chr006_GWAS.vcf chr011_GWAS.vcf chr012_GWAS.vcf chr017_GWAS.vcf chr018_GWAS.vcf chr023_GWAS.vcf chr024_GWAS.vcf chr029_GWAS.vcf chr030_GWAS.vcf chr035_GWAS.vcf chr036_GWAS.vcf chr041_GWAS.vcf chr042_GWAS.vcf | bgzip -c > Dlineage_gwas.vcf.gz 

3.PCA 分析：
plink --vcf ABlineage_gwas.vcf.gz --pca 3 header tabs -out ABlineage_pca --double-id --autosome-num 42 &
plink --vcf Dlineage_gwas.vcf.gz --pca 3 header tabs -out Dlineage_pca --double-id --autosome-num 42 &

4.亲缘关系矩阵
run_pipeline.pl -Xmx200g -importGuess ABlineage_gwas.vcf.gz -KinshipPlugin -method Centered_IBS -endPlugin -export ABlineage_kinship.txt -exportType SqrMatrix &
run_pipeline.pl -Xmx200g -importGuess Dlineage_gwas.vcf.gz -KinshipPlugin -method Centered_IBS -endPlugin -export Dlineage_kinship.txt -exportType SqrMatrix &

5.VCF转格式
run_pipeline.pl -Xmx200g -fork1 -vcf ABlineage_gwas.vcf.gz -export ABlineage -exportType Hapmap -runfork1 
run_pipeline.pl -Xmx200g -fork1 -vcf Dlineage_gwas.vcf.gz -export Dlineage -exportType Hapmap -runfork1 

run_pipeline.pl -Xmx200g -SortGenotypeFilePlugin -inputFile ABlineage.hmp.txt  -outputFile ABlineage_sort -fileType Hapmap &
run_pipeline.pl -Xmx200g -SortGenotypeFilePlugin -inputFile Dlineage.hmp.txt  -outputFile Dlineage_sort -fileType Hapmap &

6.广义线性回归
run_pipeline.pl -fork1 -h ABlineage_sort.hmp.txt -fork2 -r env1.txt -fork3 -q ABlineage_pca.qmatrix -excludeLastTrait -fork4 -k ABlineage_kinship.txt -combine5  -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export env${i}_A_mlm -runfork1 -runfork2 -runfork3

7.结果处理
for i in `ls *mlm2.txt`; do sed '1d' $i | grep -v "NaN" | awk '{print $2"\t"$3"\t"$4"\t"$7}' | sed '1i SNP\tCHR\tBP\tP'> Result/${i::-9}.txt ; done
for i in *txt
do
sed '1d' $i | shuf -n 5000 | sed '1i SNP\tCHR\tBP\tP' > ${i::-4}.qq.txt
done

8.本地画图
#-----qq----
library(ggplot2)
setwd("/Users/guoyafei/Documents/02_Vmap3/23_gwas/pdf")
path <- "/Users/guoyafei/Documents/02_Vmap3/23_gwas/plot"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
data_all <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F)})

for ( i in c(1:7)){
  p <- list()
  n1 <- strsplit( names(data_all)[i], "_")[[1]][1]
  n2 <- strsplit( names(data_all)[i], "_")[[1]][2]
  filename <- paste(n1,n2,"qq",sep="_")
    data <- data_all[[i]]
    colnames(data) <- c("Trait", "SNP", "CHR", "BP", "df", "F", "P")
    qq_dat <- data.frame(obs=-log10(sort(data$P,decreasing=FALSE)),
                         exp=-log10( ppoints(length(data$P))))
    p[[i]] <- ggplot(data=qq_dat,aes(exp,obs))+
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
    pdf(paste(filename,".pdf",sep=""),height = 3,width = 3)
    print(p[[i]])
    dev.off()
  }

#--------------manhuttun-----
library(qqman)
library(tidyverse)
setwd("/Users/guoyafei/Documents/02_Vmap3/23_gwas/pdf")
path <- "/Users/guoyafei/Documents/02_Vmap3/23_gwas/plot"
fileNames <- dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F)})
#threshold
#library(gdata)
#thresh <- read.xls("thresh.xlsx",sheet=1,row.name=1,na.strings=c("NA","#DIV/0!"))
for (i in seq(2)){
  n1 <- strsplit( names(data)[i], "_")[[1]][1]
  n2 <- strsplit( names(data)[i], "_")[[1]][2]
  filename <- paste(n1,n2,sep="_")
  p<-list()
    gwasResults <- data[[i]]
    colnames(gwasResults) <- c("Trait", "SNP", "CHR", "BP", "df", "F", "P")
    don <- gwasResults %>% 
      group_by(CHR) %>% 
      summarise(chr_len=max(BP)) %>% 
      mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
      select(-chr_len) %>%
      left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
      arrange(CHR, BP) %>%
      mutate( BPcum=BP+tot) 
    axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
    p[[i]] <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) ) +   
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10),
      )+
      scale_y_continuous(limits = c(0,7))+
      geom_hline(yintercept = -log10(0.0699), colour="red",linetype=2, size=1)
    pdf(paste(filename,".pdf",sep=""),height = 3,width = 6.5)
    print(p[[i]])
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


#31个高低海拔样本的环境GWAS分析
#服务器:SH1&SH2:
204:/data2/yafei/003_Project3/Vmap1.1/E6/GWAS
203:/data1/home/yafei/003_Project3/GWAS
#58个高低海拔样本的环境GWAS分析
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/09_GWAS/58taxa")
data <- read.table("D_mlm2.txt",header=T,stringsAsFactors = F)

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
    gwasResults3 <- gwasResults2[which(gwasResults2$P < 0.25),]
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
    p[[j+1]] <- ggplot(don, aes(x=BPcum, y=F)) +
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
      geom_point(data=point,aes(x=BPcum,y=-log10(P)),color="red")+
    #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
    geom_hline(yintercept = -log10(0.07043), colour="red",linetype=2, size=1)
    #geom_hline(yintercept = -log10(thresh[a,2]), colour="blue",linetype=2, size=1)
  }
  pdf(paste(filename,".pdf",sep=""),height = 9,width = 9)
  grid.arrange(p[[1]],p[[2]],nrow=2)
  dev.off()
}




  