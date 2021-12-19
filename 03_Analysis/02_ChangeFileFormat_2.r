#-----------------------------------------------------------提取WW Reference Allele.sh---------------------------------------------------------------------------------------------------------
#使用bedtools提取WW的SNP位点的Reference Allele
for chr in {1..42}
do
bedtools getfasta -fi /data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa -bed chr${chr}.bed | sed 's/.*-//' |sed -n '{N;s/\n/ /p}' > WW_site/chr${chr}_WW.site
done
#测试bedtools是否提取到了正确的位点及正确的Reference Allele
#工作目录：203:/data1/home/yafei/SNP/WW
bedtools getfasta -fi /data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa -bed chr2.bed  | sed 's/.*-//' |xargs -n 2 > chr2_WW.site
zcat /data1/home/yafei/Project3/Vmap1.1/chr02.all.vcf.gz | awk '{if (NF>20) print $2,$4}' > /data1/home/yafei/SNP/WW/chr2_ref.site
cat chr2_WW_1.site /data1/home/yafei/SNP/WW/chr2_ref.site | sort |uniq -c|grep '^2' |head
awk '{print $4}' chr2.tab > test.tab
cat test.tab /data1/home/yafei/SNP/WW/chr2_ref.site |awk '{print $1}' | sort |uniq -c|awk '{if ($1==2) print $0}' |wc  #642

#-----------------------------------------------------转换WW.tab文件为vcf文件.jar(Tab2vcf.jar)---------------------------------------------------------------------------------------------------
#输入文件
#/data1/home/yafei/SNP/WWWW_site/chr${chr}_WW.site: 该染色体WW检测到的SNP位点以及该位点的Reference Allele。格式(snpPosition Ref_allele)，没有表头。
#/data1/home/yafei/SNP/WW/chr${chr}.tab: 该染色体原始文件，表头为(psId\tsnpId\tchromosome\tsnpPosition\tsnpQuality\tERGE23958\tERGE24682...)
#/data1/home/yafei/SNP/WW/chr${chr}.convers: 该染色体 A, B allele的对应碱基文件。表头为(psId\tsnpId\tchromosome\tsnpPosition\tA\tB)
#输出文件
#/data1/home/yafei/SNP/WW/outWWvcf_E2/chr${chr}.vcf
#检查vcf文件的正确性
#awk '{if(NF>100) print $3"\t"$4"\t"$5"\t"$10}' chr36.vcf |head
#head WW_site/chr36_WW.site
#head chr36.bed
#head chr36.convers
#awk '{print $4"\t"$6}' chr36.tab |head
#MAF分布

#-----------------------------------------------------以WW的位点为库，对Vmap1进行HapScanner.sh---------------------------------------------------------------------------------------------------
#生成HapScanner的pos以及posAllele文件
#203服务器：/data1/home/yafei/SNP/WW/outWWvcf/
for chr in {1..42}
do
awk '{if(NF > 20) print $1"\t"$2}' chr${chr}.vcf | sed '1d' > pos/chr${chr}_pos.txt
awk '{if(NF > 20) print $1"\t"$2"\t"$4"\t"$5}' chr${chr}.vcf | sed '1c Chr\tPos\tRef\tAlt' > posAllele/chr${chr}_posAllele.txt
done
#运行HapScanner
#204服务器xuebo:/data2/xuebo/Projects/Speciation/More_accessions/hapScan/para_file/bash**.sh
#HapScanner的参数文件para_file的生成：para.jar

#---------------------------------------------WW: 合并vcf文件成Alineage.vcf/Blineage.vcf/Dlineage.vcf并做基本统计.sh---------------------------------------------------------------------------------------------------
#合并lineage文件
nohup vcf-concat chr1.vcf chr2.vcf chr7.vcf chr8.vcf chr13.vcf chr14.vcf chr19.vcf chr20.vcf chr25.vcf chr26.vcf chr31.vcf chr32.vcf chr37.vcf chr38.vcf > Alineage.vcf &
  nohup vcf-concat chr3.vcf chr4.vcf chr9.vcf chr10.vcf chr15.vcf chr16.vcf chr21.vcf chr22.vcf chr27.vcf chr28.vcf chr33.vcf chr34.vcf chr39.vcf chr40.vcf > Blineage.vcf &
  nohup vcf-concat chr5.vcf chr6.vcf chr11.vcf chr12.vcf chr17.vcf chr18.vcf chr23.vcf chr24.vcf chr29.vcf chr30.vcf chr35.vcf chr36.vcf chr41.vcf chr42.vcf > Dlineage.vcf &
  
  bgzip -c Alineage.vcf > Alineage.vcf.gz &
  bgzip -c Blineage.vcf > Blineage.vcf.gz &
  bgzip -c Dlineage.vcf > Dlineage.vcf.gz &
  
  #统计
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data1/home/yafei/SNP/WW/outWWvcf/Alineage.vcf.gz --out Alineage_siteQCfile.txt --out2 Alineage_taxaQCfile.txt > nohupA 2>& 1 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data1/home/yafei/SNP/WW/outWWvcf/Blineage.vcf.gz --out Blineage_siteQCfile.txt --out2 Blineage_taxaQCfile.txt > nohupB 2>& 1 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data1/home/yafei/SNP/WW/outWWvcf/Dlineage.vcf.gz --out Dlineage_siteQCfile.txt --out2 Dlineage_taxaQCfile.txt > nohupD 2>& 1 

#---------------------------------------------WW: 42个vcf文件做基本统计并合并成AA/BB/DD/AABB/AABBDD不同家系统计情况.sh---------------------------------------------------------------------------------------------------
for chr in {1..42}
do
bgzip -c chr${chr}.vcf > chr${chr}.vcf.gz
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data1/home/yafei/SNP/WW/outWWvcf/chr${chr}.vcf.gz --out chr${chr}_siteQCfile.txt --out2 chr${chr}_taxaQCfile.txt &
  done

#---------------------------------------------------------WW: 检查原始WW文件的42个vcf文件基本统计情况.sh & .R---------------------------------------------------------------------------------------------------
for chr in {1..42}
do
bgzip -c chr${chr}.vcf > chr${chr}.vcf.gz
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data1/home/yafei/SNP/WW/allWWvcf/chr${chr}.vcf.gz --out chr${chr}_siteQCfile.txt --out2 chr${chr}_taxaQCfile.txt &
  done

for chr in {1..42}
do
cat chr${chr}_siteQCfile.txt >> siteQCfile.txt
done
sed -i  '/^Chr/d' siteQCfile.txt
sed -i  '1i Chr\tPos\tHeterozygousProportion\tMissingRate\tMaf' siteQCfile.txt

#R: 检查原始文件的Maf和Missing Rate分布
library(ggplot2)
data <- read.table("siteQCfile.txt",header=T,stringsAsFactors=F)

for(i in c(1,3,5,7,9,11,13,17,19,21,23,25,27,29,31,33,35,37,39,41)){
  data[which(data$Chr == i  | data$Chr== i+1 ),1]  <- (i+1)/2
}
#labels <- c(1 = "1A", 2 = "1B", 3 = "1D", 4 = "2A", 5 = "2B", 6 = "2D",7 = "3A", 8 = "3B", 9 = "3D", 10 = "4A", 11 = "4B", 12 = "4D",13 = "5A", 14 = "5B", 15 = "5D",16 = "6A", 17 = "6B", 18 = "6D",19 = "7A", 20 = "7B", 21 = "7D")
chr <- c("1A", "1B", "1D", "2A", "2B",  "2D","3A",  "3B", "3D",  "4A", "4B", "4D","5A", "5B", "5D","6A",  "6B", "6D","7A",  "7B", "7D")
for(i in 1:21){
  data[data$Chr == i,1] <- chr[i]
}
pdf("AllSite_Maf.pdf",width=20,height=8)
ggplot(data, aes(x=Maf))+ geom_histogram() + facet_grid(. ~ Chr) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
pdf("AllSite_MissRate.pdf",width=20,height=8)
ggplot(data, aes(x=MissingRate))+ geom_histogram() + facet_grid(. ~ Chr) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#sh: 检查40k位点缺失的原因
#原因是：之前合并三个文件，WW里的位点如果在其他两个里面存在（1.1）并且变异含有Ref Allele（1.2）则保留，或者WW里的位点在其他两个里面不存在（2），则保留，不考虑是否含有Ref Allele。
#所以通过是否含有ref allele的条件的过滤，过滤掉了40K的位点。