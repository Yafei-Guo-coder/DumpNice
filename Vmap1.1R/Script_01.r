#------------------------------------------------------------------转换WW原始文件.R-------------------------------------------------------------------------------------------------------------
#读取文件：203服务器
setwd("/data1/home/yafei/SNP/WW")
all <- read.table("Balfourier_et_al_Wheat_Phylogeography_DataS2.tab",header=T,stringsAsFactors=F)
names <- read.table("741names.txt",header=F,stringsAsFactors=F)
convers <- read.table("tabw280kAlleleConversionTable.tab",header=T,stringsAsFactors=F)
geno <- all[,which(colnames(all) %in% names[,1])]
header <- all[,1:5]
data <- rbind(header,geno)
#write.table(alldata,"741.tab",row.names=F,quote=F,sep="\t")
#data <- read.table("741.tab",header=T,stringsAsFactors=F)
#修改染色体号和snp位点
chr <- paste(rep(paste("chr",1:7,sep=""),each=3),rep(c("A","B","D"),times=7), sep = "")
pos <- c(471304005,438720154,452179604, 462376173, 453218924,462216879,454103970,448155269,476235359,452555092,451014251,451004620,453230519,451372872,451901030,452440856,452077197,450509124,450046986,453822637,453812268)
for(i in 1:21){
  data[which(data$chromosome == chr[i]),3]  <- i*2-1
  data[which(data$chromosome == i*2-1 & data$snpPosition > pos[i]),3] <- i*2
  data[which(data$chromosome == i*2 & data$snpPosition > pos[i]),4] <- data[which(data$chromosome == i*2),4] - pos[i]
  orderd <- data[order(data[,3], data[,4]), ]
}
sub <- merge(data[,c(1:4)], convers,by.x="psId",by.y="probesetId",all.x=TRUE)
sub_orderd <- sub[order(sub[,3], sub[,4]), ]
#写42个(.convers)文件
varName <- paste("chr",1:42,".convers",sep="")
for(i in 1:42){
  write.table(sub_orderd[which(sub_orderd$chromosome == i),], varName[i],row.names=F,quote=F,sep="\t")
}
#写42个(.tab)文件
varName <- paste("chr",1:42,".tab",sep="")
for(i in 1:42){
  write.table(orderd[which(orderd$chromosome == i),], varName[i],row.names=F,quote=F,sep="\t")
}
#写42个(.bed)文件
varName <- paste("chr",1:42,".bed",sep="")
for(i in 1:42){
  var <- cbind(orderd[which(orderd$chromosome == i),3],orderd[which(orderd$chromosome == i),4],orderd[which(orderd$chromosome == i),4]+1)
  write.table(var, varName[i],row.names=F,col.names=F,quote=F,sep="\t")
}

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
head posAllele_E1/chr1_posAllele.txt 
head chr1.bed
head chr1.convers 
head WW_site/chr1_WW.site 
less -LS Trancing.SNPs.vcf
less -LS Exome.vcf
less -LS output_World.vcf
#原因是：之前合并三个文件，WW里的位点如果在其他两个里面存在（1.1）并且变异含有Ref Allele（1.2）则保留，或者WW里的位点在其他两个里面不存在（2），则保留，不考虑是否含有Ref Allele。
#所以通过是否含有ref allele的条件的过滤，过滤掉了40K的位点。

#--------------------------------------------合并WW和Vmap数据：(AA/BB/DD/AABB)/AABBDD/Landrace/ABDlineages并做基本统计.sh------------------------------------------------------------------------------------------------------

#把hascan结果文件从204:/data2/xuebo/Projects/Speciation/More_accessions/hapScan/out/ 移动到203:/data2/yafei/Project3/Vmap1.1/Hapscan/
#文件位置做一些调整
#for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40,5,6,11,12,17,18,23,24,29,30,35,36,41,42}; do mv chr${i}/VCF/ ./; done
#for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40,5,6,11,12,17,18,23,24,29,30,35,36,41,42}; do rm -r chr${i}; done
#cp -r outWWvcf/ /data2/yafei/Project3/WW

#目前工作目录：203: /data2/yafei/Project3/
#对Vmap1.1和WW的所有vcf文件按照染色体进行合并，然后挑选样本，做maf过滤，提取出分离位点，再输出到不同的品系中（AA/BB/DD/AABB/AABBDD/Landrace）。再将其按照染色体合并成不同的lineage（A lineage，B lineage，D lineage）。

#压缩样本并建立索引
for i in `ls`
do 
bgzip -c ${i} > ${i}.gz
tabix -p vcf ${i}.gz
done


#merge_maf.sh
#合并样本
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
bcftools merge /data2/yafei/Project3/Vmap1.1/Hapscan/AA/chr0${chr}.vcf.gz /data2/yafei/Project3/Vmap1.1/Hapscan/AABB/chr0${chr}.vcf.gz /data2/yafei/Project3/Vmap1.1/Hapscan/AABBDD/chr0${chr}.vcf.gz /data2/yafei/Project3/WW_E2/chr${chr}.vcf.gz -o /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf 
done
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
bcftools merge /data2/yafei/Project3/Vmap1.1/Hapscan/BB/chr0${chr}.vcf.gz /data2/yafei/Project3/Vmap1.1/Hapscan/AABB/chr0${chr}.vcf.gz /data2/yafei/Project3/Vmap1.1/Hapscan/AABBDD/chr0${chr}.vcf.gz /data2/yafei/Project3/WW_E2/chr${chr}.vcf.gz -o /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf 
done
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
bcftools merge /data2/yafei/Project3/Vmap1.1/Hapscan/DD/chr0${chr}.vcf.gz /data2/yafei/Project3/Vmap1.1/Hapscan/AABBDD/chr0${chr}.vcf.gz /data2/yafei/Project3/WW_E2/chr${chr}.vcf.gz -o /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf 
done

#过滤maf提取样本
#AA
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/AA_taxa.txt --maf 0.0000001 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/AA/chr${chr}.mafAA.vcf.gz 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_maf/AA/chr${chr}.mafAA.vcf.gz --out /data2/yafei/Project3/V+W_E2/E2_maf/AA/chr${chr}_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_maf/AA/chr${chr}_taxaQCfile.txt > nohupA 2>& 1 
done
#BB
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/chr${chr}.all.vcf.gz --keep /data2/yafei/Project3/group/SS_taxa.txt --max-missing 0.1 --maf 0.0000001 --recode --recode-INFO-all --stdout | bgzip -c > /data1/home/yafei/Project3/Maf/mafBB/chr${chr}.mafBB.vcf.gz 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data1/home/yafei/Project3/Maf/mafBB/chr${chr}.mafBB.vcf.gz --out /data1/home/yafei/Project3/Maf/mafBB/chr${chr}_siteQCfile.txt --out2 /data1/home/yafei/Project3/Maf/mafBB/chr${chr}_taxaQCfile.txt > nohupB 2>& 1 
done
#DD
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/DD_taxa.txt --maf 0.0000001 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/DD/chr${chr}.mafDD.vcf.gz 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_maf/DD/chr${chr}.mafDD.vcf.gz --out /data2/yafei/Project3/V+W_E2/E2_maf/DD/chr${chr}_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_maf/DD/chr${chr}_taxaQCfile.txt > nohupD 2>& 1 
done
#AABB
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/AABB_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/AABB/chr${chr}.mafAABB.vcf.gz 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_maf/AABB/chr${chr}.mafAABB.vcf.gz --out /data2/yafei/Project3/V+W_E2/E2_maf/AABB/chr${chr}_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_maf/AABB/chr${chr}_taxaQCfile.txt > nohupAB 2>& 1 
done
#AABBDD
for chr in {1..42}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/AABBDD_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/AABBDD/chr${chr}.mafAABBDD.vcf.gz 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_maf/AABBDD/chr${chr}.mafAABBDD.vcf.gz --out /data2/yafei/Project3/V+W_E2/E2_maf/AABBDD/chr${chr}_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_maf/AABBDD/chr${chr}_taxaQCfile.txt > nohupABD 2>& 1 
done
#Landrace
for chr in {1..42}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/subspecies/sub_Landrace_new.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/Landrace/chr${chr}.mafLandrace.vcf.gz 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_maf/Landrace/chr${chr}.mafLandrace.vcf.gz --out /data2/yafei/Project3/V+W_E2/E2_maf/Landrace/chr${chr}_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_maf/Landrace/chr${chr}_taxaQCfile.txt > nohupLand 2>& 1 
done
#nohup sh merge_maf.sh &


#合并lineage文件
vcf-concat chr1.all.vcf chr2.all.vcf chr7.all.vcf chr8.all.vcf chr13.all.vcf chr14.all.vcf chr19.all.vcf chr20.all.vcf chr25.all.vcf chr26.all.vcf chr31.all.vcf chr32.all.vcf chr37.all.vcf chr38.all.vcf > lineage/Alineage.vcf 
vcf-concat chr3.all.vcf chr4.all.vcf chr9.all.vcf chr10.all.vcf chr15.all.vcf chr16.all.vcf chr21.all.vcf chr22.all.vcf chr27.all.vcf chr28.all.vcf chr33.all.vcf chr34.all.vcf chr39.all.vcf chr40.all.vcf > lineage/Blineage.vcf 
vcf-concat chr5.all.vcf chr6.all.vcf chr11.all.vcf chr12.all.vcf chr17.all.vcf chr18.all.vcf chr23.all.vcf chr24.all.vcf chr29.all.vcf chr30.all.vcf chr35.all.vcf chr36.all.vcf chr41.all.vcf chr42.all.vcf > lineage/Dlineage.vcf 

#过滤maf
vcftools --vcf Alineage.vcf --maf 0.0000001 --recode --stdout | bgzip -c > mafAlineage.vcf.gz 
vcftools --vcf Blineage.vcf --maf 0.0000001 --recode --stdout | bgzip -c > mafBlineage.vcf.gz
vcftools --vcf Dlineage.vcf --maf 0.0000001 --recode --stdout | bgzip -c > mafDlineage.vcf.gz

#统计
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_all/lineage/mafAlineage.vcf. --out /data2/yafei/Project3/V+W_E2/E2_all/lineage/Alineage_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_all/lineage/Alineage_taxaQCfile.txt > nohupA 2>& 1 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_all/lineage/mafBlineage.vcf. --out /data2/yafei/Project3/V+W_E2/E2_all/lineage/Blineage_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_all/lineage/Blineage_taxaQCfile.txt > nohupB 2>& 1 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_all/lineage/mafDlineage.vcf. --out /data2/yafei/Project3/V+W_E2/E2_all/lineage/Dlineage_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_all/lineage/Dlineage_taxaQCfile.txt > nohupD 2>& 1 

#---------------------------------------------------plot A/B/Dlineage Taxa/Site MissRate/Maf.R---------------------------------------------------------------------------------------------------
library(ggplot2)
#ABDlineage site
A <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Alineage_siteQCfile.txt",header=T,stringsAsFactors=F)
B <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Blineage_siteQCfile.txt",header=T,stringsAsFactors=F)
D <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Dlineage_siteQCfile.txt",header=T,stringsAsFactors=F)
A$Lineage <- "A-lineage"
B$Lineage <- "B-lineage"
D$Lineage <- "D-lineage"
#colnames(A) <- c("Lineage","A","B","MissingRate","Maf")
#colnames(B) <- c("Lineage","A","B","MissingRate","Maf")
#colnames(D) <- c("Lineage","A","B","MissingRate","Maf")
site <- rbind(A,B,D)
pdf("LineageSite_Maf.pdf",width=20,height=8)
ggplot(site, aes(x=Maf))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
pdf("LineageSite_MissRate.pdf",width=20,height=8)
ggplot(site, aes(x=MissingRate))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#ABDlineage taxa
A <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Alineage_taxaQCfile.txt",header=T,stringsAsFactors=F)
B <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Blineage_taxaQCfile.txt",header=T,stringsAsFactors=F)
D <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Dlineage_taxaQCfile.txt",header=T,stringsAsFactors=F)
A$Taxa <- "A-lineage"
B$Taxa <- "B-lineage"
D$Taxa <- "D-lineage"
taxa <- rbind(A,B,D)
colnames(taxa)[1] <- "Lineage"
pdf("LineageMind_MissRate.pdf",width=20,height=8)
ggplot(taxa, aes(x=MissRate))+ geom_histogram() + labs(x="Individual Missing Rate", y="Count")+facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#-----------------------------------------------plot AA/BB/DD/AABB/AABBDD/Landrace Taxa/Site MissRate/Maf.R------------------------------------------------------------------------------------------------------------
#AA
nohup cat chr1_siteQCfile.txt chr2_siteQCfile.txt chr7_siteQCfile.txt chr8_siteQCfile.txt chr13_siteQCfile.txt  chr14_siteQCfile.txt chr19_siteQCfile.txt  chr20_siteQCfile.txt  chr25_siteQCfile.txt chr26_siteQCfile.txt  chr31_siteQCfile.txt  chr32_siteQCfile.txt chr37_siteQCfile.txt chr38_siteQCfile.txt > AA_siteQCfile.txt &
  #AABB
  nohup cat chr10_siteQCfile.txt  chr16_siteQCfile.txt  chr21_siteQCfile.txt  chr27_siteQCfile.txt  chr32_siteQCfile.txt  chr38_siteQCfile.txt  chr4_siteQCfile.txt chr13_siteQCfile.txt  chr19_siteQCfile.txt  chr22_siteQCfile.txt  chr28_siteQCfile.txt  chr33_siteQCfile.txt  chr39_siteQCfile.txt  chr7_siteQCfile.txt chr14_siteQCfile.txt  chr1_siteQCfile.txt   chr25_siteQCfile.txt  chr2_siteQCfile.txt   chr34_siteQCfile.txt  chr3_siteQCfile.txt   chr8_siteQCfile.txt chr15_siteQCfile.txt  chr20_siteQCfile.txt  chr26_siteQCfile.txt  chr31_siteQCfile.txt  chr37_siteQCfile.txt  chr40_siteQCfile.txt  chr9_siteQCfile.txt > AABB_siteQCfile.txt &
  #AABBDD
  nohup cat chr10_siteQCfile.txt  chr16_siteQCfile.txt  chr21_siteQCfile.txt  chr27_siteQCfile.txt  chr32_siteQCfile.txt  chr38_siteQCfile.txt  chr4_siteQCfile.txt chr11_siteQCfile.txt  chr17_siteQCfile.txt  chr22_siteQCfile.txt  chr28_siteQCfile.txt  chr33_siteQCfile.txt  chr39_siteQCfile.txt  chr5_siteQCfile.txt chr12_siteQCfile.txt  chr18_siteQCfile.txt  chr23_siteQCfile.txt  chr29_siteQCfile.txt  chr34_siteQCfile.txt  chr3_siteQCfile.txt   chr6_siteQCfile.txt chr13_siteQCfile.txt  chr19_siteQCfile.txt  chr24_siteQCfile.txt  chr2_siteQCfile.txt   chr35_siteQCfile.txt  chr40_siteQCfile.txt  chr7_siteQCfile.txt chr14_siteQCfile.txt  chr1_siteQCfile.txt   chr25_siteQCfile.txt  chr30_siteQCfile.txt  chr36_siteQCfile.txt  chr41_siteQCfile.txt  chr8_siteQCfile.txt chr15_siteQCfile.txt  chr20_siteQCfile.txt  chr26_siteQCfile.txt  chr31_siteQCfile.txt  chr37_siteQCfile.txt  chr42_siteQCfile.txt  chr9_siteQCfile.txt > AABBDD_siteQCfile.txt &
  #DD
  nohup cat chr11_siteQCfile.txt  chr17_siteQCfile.txt  chr23_siteQCfile.txt  chr29_siteQCfile.txt  chr35_siteQCfile.txt  chr41_siteQCfile.txt  chr5_siteQCfile.txt chr12_siteQCfile.txt  chr18_siteQCfile.txt  chr24_siteQCfile.txt  chr30_siteQCfile.txt  chr36_siteQCfile.txt  chr42_siteQCfile.txt  chr6_siteQCfile.txt  > DD_siteQCfile.txt &
  #Landrace
  nohup cat chr10_siteQCfile.txt  chr16_siteQCfile.txt  chr21_siteQCfile.txt  chr27_siteQCfile.txt  chr32_siteQCfile.txt  chr38_siteQCfile.txt  chr4_siteQCfile.txt chr11_siteQCfile.txt  chr17_siteQCfile.txt  chr22_siteQCfile.txt  chr28_siteQCfile.txt  chr33_siteQCfile.txt  chr39_siteQCfile.txt  chr5_siteQCfile.txt chr12_siteQCfile.txt  chr18_siteQCfile.txt  chr23_siteQCfile.txt  chr29_siteQCfile.txt  chr34_siteQCfile.txt  chr3_siteQCfile.txt   chr6_siteQCfile.txt chr13_siteQCfile.txt  chr19_siteQCfile.txt  chr24_siteQCfile.txt  chr2_siteQCfile.txt   chr35_siteQCfile.txt  chr40_siteQCfile.txt  chr7_siteQCfile.txt chr14_siteQCfile.txt  chr1_siteQCfile.txt   chr25_siteQCfile.txt  chr30_siteQCfile.txt  chr36_siteQCfile.txt  chr41_siteQCfile.txt  chr8_siteQCfile.txt chr15_siteQCfile.txt  chr20_siteQCfile.txt  chr26_siteQCfile.txt  chr31_siteQCfile.txt  chr37_siteQCfile.txt  chr42_siteQCfile.txt  chr9_siteQCfile.txt Landrace_siteQCfile.txt &
  #SS
  nohup cat chr10_siteQCfile.txt  chr16_siteQCfile.txt  chr22_siteQCfile.txt  chr28_siteQCfile.txt  chr34_siteQCfile.txt  chr3_siteQCfile.txt   chr4_siteQCfile.txt chr15_siteQCfile.txt  chr21_siteQCfile.txt  chr27_siteQCfile.txt  chr33_siteQCfile.txt  chr39_siteQCfile.txt  chr40_siteQCfile.txt  chr9_siteQCfile.txt > SS_siteQCfile.txt &
  
  #site_miss & site_maf
  AA <- read.table("/data1/home/yafei/Project3/Maf/mafAA/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
AA_maf <- AA[,c(1,5)]
AA_miss <- AA[,c(1,4)]
AA_maf$V1 <- "AA"
AA_miss$V1 <- "AA"

SS <- read.table("/data1/home/yafei/Project3/Maf/mafSS/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
SS_maf <- SS[,c(1,5)]
SS_miss <- SS[,c(1,4)]
SS_maf$V1 <- "SS"
SS_miss$V1 <- "SS"

DD <- read.table("/data1/home/yafei/Project3/Maf/mafDD/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
DD_maf <- DD[,c(1,5)]
DD_miss <- DD[,c(1,4)]
DD_maf$V1 <- "DD"
DD_miss$V1 <- "DD"

AABB <- read.table("/data1/home/yafei/Project3/Maf/mafAABB/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
AABB_maf <- AABB[,c(1,5)]
AABB_miss <- AABB[,c(1,4)]
AABB_maf$V1 <- "AABB"
AABB_miss$V1 <- "AABB"

AABBDD <- read.table("/data1/home/yafei/Project3/Maf/mafAABBDD/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
AABBDD_maf <- AABBDD[,c(1,5)]
AABBDD_miss <- AABBDD[,c(1,4)]
AABBDD_maf$V1 <- "AABBDD"
AABBDD_miss$V1 <- "AABBDD"

Landrace <- read.table("/data1/home/yafei/Project3/Maf/mafLand/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
Landrace_maf <- Landrace[,c(1,5)]
Landrace_miss <- Landrace[,c(1,4)]
Landrace_maf$V1 <- "Landrace"
Landrace_miss$V1 <- "Landrace"

maf <- rbind(AA_maf,SS_maf,DD_maf,AABB_maf,AABBDD_maf)
colnames(maf) <- c("Lineage","Maf")
maf$Maf <- as.numeric(maf$Maf)

miss <- rbind(AA_miss,SS_miss,DD_miss,AABB_miss,AABBDD_miss)
colnames(miss) <- c("Lineage","MissRate")
miss$MissRate <- as.numeric(miss$MissRate)

library(ggplot2)
pdf("site_Maf.pdf",width=20,height=8)
ggplot(maf, aes(x=Maf))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=25),panel.border = element_blank())
dev.off()

pdf("site_MissRate.pdf",width=20,height=8)
ggplot(miss, aes(x=MissRate))+ geom_histogram() + labs(x="Site Missing Rate", y="count")+facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=25),panel.border = element_blank())
dev.off()

**mind_miss
AA <- read.table("/data1/home/yafei/Project3/Maf/mafAA/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
SS <- read.table("/data1/home/yafei/Project3/Maf/mafSS/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
DD <- read.table("/data1/home/yafei/Project3/Maf/mafDD/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
AABB <- read.table("/data1/home/yafei/Project3/Maf/mafAABB/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
AABBDD <- read.table("/data1/home/yafei/Project3/Maf/mafAABBDD/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
Landrace <- read.table("/data1/home/yafei/Project3/Maf/mafLand/all_taxaQCfile.txt",header=F,stringsAsFactors=F)

AA$V1 <- "AA"
SS$V1 <- "SS"
DD$V1 <- "DD"
AABB$V1 <- "AABB"
AABBDD$V1 <- "AABBDD"
Landrace$V1 <- "Landrace"

mind<- rbind(AA,SS,DD,AABB,AABBDD,Landrace)
colnames(mind) <- c("Lineage","Het","MissingRate")
mind$MissingRate <- as.numeric(mind$MissingRate)

**应该画密度
pdf("Mind_MissRate.pdf",width=20,height=8)
ggplot(mind, aes(x=MissingRate,y=..density..))+ geom_histogram() + labs(x="Individual Missing Rate", y="Frequency")+facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#------------------------------------------------------------------------计算FST.sh--------------------------------------------------------------------------------------
#FST_group: /data2/yafei/Project3/FST_group/subgroups/ 共27个groups。
#工作目录：/data2/yafei/Project3/FST_group/VmapData/
#AABB
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/AABB_AABBDD_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/AABB_AABBDD/chr${chr}.mafAABB_AABBDD.vcf.gz 
done
vcf-concat *mafAABB_AABBDD.vcf.gz  > mafAABB_AABBDD.vcf 
bgzip -c mafAABB_AABBDD.vcf > mafAABB_AABBDD.vcf.gz

#pop1:	A_Wild_einkorn	
#pop2:	A_Domesticated_einkorn	
#pop3:	A_Urartu
#pop4:	B_Speltoides
#pop5:	AABB_Wild_emmer
#pop6:	AABB_Domesticated_emmer
#pop7:	AABB_Georgian_wheat
#pop8:	AABB_Ispahanicum
#pop9:	AABB_Rivet_wheat
#pop10:	AABB_Polish_wheat
#pop11:	AABB_Persian_wheat
#pop12:	AABB_Khorasan_wheat
#pop13:	AABB_Durum
#pop14:	AABBDD_Spelt
#pop15:	AABBDD_Macha
#pop16:	AABBDD_Club_wheat
#pop17:	AABBDD_Indian_dwarf_wheat
#pop18:	AABBDD_Yunan_wheat
#pop19:	AABBDD_Xinjiang_wheat
#pop20:	AABBDD_Tibetan_semi_wild
#pop21:	AABBDD_Vavilovii
#pop22:	AABBDD_Landrace_new
#pop23:	AABBDD_Cultivar
#pop24:	AABBDD_Synthetic
#pop25:	DD_Strangulata
#pop26:	DD_Meyeri
#pop27:	DD_Anathera

#A: /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz
#B: /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz
#D: /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz
#AABB: /data2/yafei/Project3/V+W_E2/E2_maf/AABB_AABBDD/mafAABB_AABBDD.vcf.gz
#AABBDD: /data2/yafei/Project3/V+W_E2/E2_maf/AABBDD/AABBDD.vcf.gz

#分lineage计算fst: /data2/yafei/Project3/FST_group/lineage/
#Alineage
#pop1:
for i in {2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/pop1.txt --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/A/"p_1_"$i 
done

#pop2:
for i in {3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/pop2.txt --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/A/"p_2_"$i &
  done
#pop3:
for i in {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/pop3.txt --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/A/"p_3_"$i &
  done

#pop5-24:

thread_num=20
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

for i in {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
do
for((j=i+1;j<=24;j++))
  do
{
  read -u6
  {
    vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$j".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/A/"p_"$i"_"$j 
    echo "" >&6
  } & 
} 
done 
done
wait

exec 6>&-
  
  
  #Blineage
  #pop5-24:
  
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

for i in {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
do
for((j=i+1;j<=24;j++))
  do
{
  read -u6
  {
    vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$j".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/B/"p_"$i"_"$j > Bnohup$i 2>& 1
    echo "" >&6
  } & 
} 
done 
done

wait

exec 6>&-
  
  #Dlineage
  #pop14-27:
  thread_num=60
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
for i in {14,15,16,17,18,19,20,21,22,23,24,25,26}
do
for((j=i+1;j<=27;j++))
  do
{
  read -u6
  {
    vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$j".txt" --out /data2/yafei/Project3/FST_group/VmapData/D/site_fst/"p_"$i"_"$j >> Dnohup$i 2>& 1
    echo "" >&6
  } & 
} 
done 
done

wait
exec 6>&-
  
  #grep "out p" nohup.out | awk '{print $2}'
  #grep "weight" nohup.out | awk '{print $7}'
  
  #A/B/D.fst 修改输出文件名
  for i in `ls *windowed.weir.fst`
do
awk 'BEGIN {count = 0} {if($5<10){sum+=$5; count++}} END {print FILENAME,"Average = ", sum/count}' $i >> A_weight.fst
done

#求中位数
for i in `ls *windowed.weir.fst`
do
echo $i
awk '{print $5}' $i |sed '1d' |sort |awk -f medium.awk
done |xargs -n 2 > D_medium.fst

#AB/ABD.fst
for i in `cat ABfiles`
do
echo $i 
cat /data2/yafei/Project3/FST_group/VmapData/A/$i  /data2/yafei/Project3/FST_group/VmapData/B/$i | sed '/WEIGHTED/d'|awk 'BEGIN {count = 0} {if($5<10){sum+=$5; count++}} END {print "Average = ", sum/count}' >> AB_weight.fst 
done 

#画热图 修改文件名，数据集名，输出文件名
#yafei@203:/data2/yafei/Project3/FST_group/VmapData/heatmap

library("corrplot")
data <- read.table("D_weight.fst_ibs",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$Dlineage
data<- data[,-1]
data<- as.matrix(data)
pdf("fst_ibs_D.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()

#----------------------------------------------------------计算Vmap和WW的AABBDD的分离位点并进行统计画图 .sh & .R--------------------------------------------------------------------------------------
#203: /data2/yafei/Project3/V+W/E2_maf/Vmap
#WW
#AABBDD
for chr in {1..42}
do
vcftools --vcf /data2/yafei/Project3/V+W/E2_all_ChangeName/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/WW_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W/E2_maf/Vmap/chr${chr}.WWmafAABBDD.vcf.gz 
done

vcf-concat *WWmafAABBDD.vcf.gz  > WWmafAABBDD.vcf 
bgzip -c WWmafAABBDD.vcf > WWmafAABBDD.vcf.gz

#Vmap
#AABBDD
for chr in {1..42}
do
vcftools --vcf /data2/yafei/Project3/V+W/E2_all_ChangeName/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/AABBDD_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W/E2_maf/Vmap/chr${chr}.VmapmafAABBDD.vcf.gz 
done
vcf-concat *VmapmafAABBDD.vcf.gz  > VmapmafAABBDD.vcf 
bgzip -c VmapmafAABBDD.vcf > VmapmafAABBDD.vcf.gz

#统计
nohup java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W/E2_maf/Vmap/VmapmafAABBDD.vcf.gz --out /data2/yafei/Project3/V+W/E2_maf/Vmap/Vmap_siteQCfile.txt --out2 /data2/yafei/Project3/V+W/E2_maf/Vmap/Vmap_taxaQCfile.txt > nohupVmap 2>& 1 &
  nohup java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W/E2_maf/Vmap/WWmafAABBDD.vcf.gz --out /data2/yafei/Project3/V+W/E2_maf/Vmap/WW_siteQCfile.txt --out2 /data2/yafei/Project3/V+W/E2_maf/Vmap/WW_taxaQCfile.txt > nohupWW 2>& 1 &
  
  #画图
  library(ggplot2)
#site frequency
WW <- read.table("/data2/yafei/Project3/V+W/E2_maf/Vmap/WW_siteQCfile.txt",header=T,stringsAsFactors=F)
Vmap <- read.table("/data2/yafei/Project3/V+W/E2_maf/Vmap/Vmap_siteQCfile.txt",header=T,stringsAsFactors=F)

WW$Lineage <- "WW"
Vmap$Lineage <- "Vmap"

site <- rbind(WW,Vmap)
pdf("Site_Maf.pdf",width=20,height=8)
ggplot(site, aes(x=Maf))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
pdf("Site_MissRate.pdf",width=20,height=8)
ggplot(site, aes(x=MissingRate))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#indi site
Vmap$ID <- paste(Vmap$Chr, Vmap$Pos,sep="-")
WW$ID <- paste(WW$Chr, WW$Pos,sep="-")
WW_a <- WW[which(WW$ID %in% Vmap$ID),]
Vmap_a <- Vmap[which(Vmap$ID %in% WW$ID),]
WW_a$Index <- c(1:60969)
Vmap_a$Index <- c(1:60969)
pdf("WW_Maf.pdf",width=20,height=8)
ggplot(data = WW_a, mapping = aes(x = Index, y = Maf)) + geom_bar(stat = 'identity')+theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
pdf("Vmap_Maf.pdf",width=20,height=8)
ggplot(data = Vmap_a, mapping = aes(x = Index, y = Maf)) + geom_bar(stat = 'identity')+theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

WW_a$Lineage <- "WW"
Vmap_a$Lineage <- "Vmap"
all <- rbind(WW_a,Vmap_a)
pdf("All_Maf.pdf",width=20,height=8)
ggplot(data = all, mapping = aes(x = Index, y = Maf)) + geom_bar(stat = 'identity')+facet_grid(Lineage ~ .)+theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#taxa
WW <- read.table("/data2/yafei/Project3/V+W/E2_maf/Vmap/WW_taxaQCfile.txt",header=T,stringsAsFactors=F)
Vmap <- read.table("/data2/yafei/Project3/V+W/E2_maf/Vmap/Vmap_taxaQCfile.txt",header=T,stringsAsFactors=F)
WW$Lineage <- "WW"
Vmap$Lineage <- "Vmap"

taxa <- rbind(WW,Vmap)
pdf("Mind_MissRate.pdf",width=20,height=8)
ggplot(taxa, aes(x=MissRate))+ geom_histogram() + labs(x="Individual Missing Rate", y="Count")+facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#-------------------------------------------------------------计算 Tajima's D .sh--------------------------------------------------------------------------------------
#调色板：display.brewer.pal(n = 8, name = 'Dark2')
#分lineage计算Tajima's D: /data2/yafei/Project3/FST_group/lineage/
#工作目录：203:/data2/yafei/Project3/TajimasD/VmapData
#Alineage
for i in {1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --TajimaD 1000000 --out /data2/yafei/Project3/TajimasD/VmapData/A/"pop"$i"_1M.group" &
  done

#Blineage
for i in {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --TajimaD 1000000 --out /data2/yafei/Project3/TajimasD/VmapData/B/"pop"$i"_1M.group" &
  done

#Dlineage
for i in {14,15,16,17,18,19,20,21,22,23,24,25,26,27}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --TajimaD 1000000 --out /data2/yafei/Project3/TajimasD/VmapData/D/"pop"$i"_1M.group" &
  done

#ABlineage

for i in `cat ABfiles`
do
cat /data2/yafei/Project3/TajimasD/VmapData/A/$i  /data2/yafei/Project3/TajimasD/VmapData/B/$i | sed '/CHROM/d' | sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > /data2/yafei/Project3/TajimasD/VmapData/AB/$i
done

#ABDlineage
for i in `cat ABDfiles`
do
cat /data2/yafei/Project3/TajimasD/VmapData/A/$i  /data2/yafei/Project3/TajimasD/VmapData/B/$i  /data2/yafei/Project3/TajimasD/VmapData/D/$i| sed '/CHROM/d' | sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > /data2/yafei/Project3/TajimasD/VmapData/ABD/$i
done

#重新计算新疆的TajimaD
for chr in {1..42}
do
nohup vcftools --gzvcf /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/chr0${chr}.vcf --keep /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/pop19.txt --TajimaD 1000000 --out /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/xinjiang_Tajima/${chr}_pop19_1M.group &
  done

nohup cat 1_pop19_1M.group.Tajima.D 2_pop19_1M.group.Tajima.D 7_pop19_1M.group.Tajima.D 8_pop19_1M.group.Tajima.D 13_pop19_1M.group.Tajima.D 14_pop19_1M.group.Tajima.D 19_pop19_1M.group.Tajima.D 20_pop19_1M.group.Tajima.D 25_pop19_1M.group.Tajima.D 26_pop19_1M.group.Tajima.D 31_pop19_1M.group.Tajima.D 32_pop19_1M.group.Tajima.D 37_pop19_1M.group.Tajima.D 38_pop19_1M.group.Tajima.D | sed '/CHROM/d' |sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > A/pop19_1M.group &
  nohup cat 3_pop19_1M.group.Tajima.D 4_pop19_1M.group.Tajima.D 9_pop19_1M.group.Tajima.D 10_pop19_1M.group.Tajima.D 15_pop19_1M.group.Tajima.D 16_pop19_1M.group.Tajima.D 21_pop19_1M.group.Tajima.D 22_pop19_1M.group.Tajima.D 27_pop19_1M.group.Tajima.D 28_pop19_1M.group.Tajima.D 33_pop19_1M.group.Tajima.D 34_pop19_1M.group.Tajima.D 39_pop19_1M.group.Tajima.D 40_pop19_1M.group.Tajima.D | sed '/CHROM/d' |sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > B/pop19_1M.group &
  nohup cat 5_pop19_1M.group.Tajima.D 6_pop19_1M.group.Tajima.D 11_pop19_1M.group.Tajima.D 12_pop19_1M.group.Tajima.D 17_pop19_1M.group.Tajima.D 18_pop19_1M.group.Tajima.D 23_pop19_1M.group.Tajima.D 24_pop19_1M.group.Tajima.D 29_pop19_1M.group.Tajima.D 30_pop19_1M.group.Tajima.D 35_pop19_1M.group.Tajima.D 36_pop19_1M.group.Tajima.D 41_pop19_1M.group.Tajima.D 42_pop19_1M.group.Tajima.D | sed '/CHROM/d' |sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > D/pop19_1M.group &
  #转移到yafei203对应位置
  
  #改名字
  library(ggplot2)
nameA <- c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar")
nameB <- c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar")
nameD <- c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera")
nameAB <- c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar")
nameABD <- c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera")

#修改 路径和name以及输出文件名
path <- "/data2/yafei/Project3/TajimasD/VmapData/ABD"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 
length(data)

#A: 4, 13
#B: 2, 11
#D: 11
#AB: 10
#ABD:
for(i in 1:length(data)){
  data[[i]]$Name <- nameABD[i]
  if (i<20) {
    data[[i]]$Fill <- "pink"
  } else if (i<20) {
    data[[i]]$Fill <- "lightblue"
  } else {
    data[[i]]$Fill <- "orange"
  }
}
all <- data.frame()
for(i in 1:length(data)){
  all <- rbind(all, data[[i]])
}
head(all)
tail(all)
#绘图：修改输出文件名
library(ggplot2)
P2<- ggplot(all, aes(x=Name, y=TajimaD,fill=Fill),na.rm=TRUE) + 
  #geom_violin() + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5, alpha = 0.8)+ #绘制箱线图
  scale_fill_manual(values = c("#6000b4", "#1f97ff", "#ffce6e"))+
  #Alineage
  #scale_x_discrete(limits= c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Blineage
  #scale_x_discrete(limits= c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Dlineage
  #scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera"))+
  #ABlineage
  #scale_x_discrete(limits= c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #ABDlineage
  scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  theme_bw()+ #背景变为白色
  ylab("TajimaD")+xlab("Lineages")+
  theme(axis.text.x=element_text(angle=80,hjust = 1,colour="black",family="Times",size=200),
        axis.text.y=element_text(family="Times",size=300,face="plain"),
        axis.title.y=element_text(family="Times",size = 300,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")

pdf(file = "ABD_TajimaD.pdf",width=200,height=100) #结果保存
print(P2)
dev.off()

#-------------------------------------------------------------计算 pi .sh--------------------------------------------------------------------------------------

#分lineage计算pi: 
#Alineage
for i in {1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do 
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --site-pi --out /data2/yafei/Project3/Pi/A/"pop"$i"1M.group" 
done

#Blineage
for i in {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --site-pi --out /data2/yafei/Project3/Pi/B/"pop"$i"1M.group"
done

#Dlineage
for i in {14,15,16,17,18,19,20,21,22,23,24,25,26,27}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --site-pi --out /data2/yafei/Project3/Pi/D/"pop"$i"1M.group"
done

#分染色体计算Windows pi

#Alineage:{1,2,7,8,13,14,19,20,25,26,31,32,37,38}
#Blineage:{3,4,9,10,15,16,21,22,27,28,33,34,39,40}
#Dlineage:{5,6,11,12,17,18,23,24,29,30,35,36,41,42}

for i in `ls *group.sites.pi`
do
for j in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
awk '{if($1=='"$j"') print $0}'  ${i} >> chr/${j}_${i}
done
done
#/data2/yafei/Project3/Pi/A/chr
#修改文件表头
for i in `ls *group.sites.pi`
do
sed -i '1i CHROM\tPOS\tPI' ${i}
done

#由于gcc版本的原因，转移到204上做
#/data1/home/yafei/Project3/Pi/A/chr
#计算windows pi
#Alineage name
#Blineage name
#Dlineage name

for i in `cat D.name`
do
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/home/yafei/Project3/Pi/D/chr/${chr}_${i} --file2 /data1/home/yafei/Project3/Pi/SNPbed_1M/chr${chr}_1M.bed --out /data1/home/yafei/Project3/Pi/D/1M_window/${chr}_${i} &
  done
done

#把chr Windows pi合并成 群体 pi

for i in `cat D.name`
do
cat D/1M_window/ *_${i} >> D/by_chr/${i}
done

cat A.name B.name |sort |uniq -c |grep "2 " |awk '{print $2}' > ABfiles
cat A.name B.name D.name|sort |uniq -c |grep "3 " |awk '{print $2}' > ABDfiles

#ABlineage
for i in `cat ABfiles`
do
cat /data1/home/yafei/Project3/Pi/A/by_chr/$i  /data1/home/yafei/Project3/Pi/B/by_chr/$i  > /data1/home/yafei/Project3/Pi/AB/$i
done
#ABDlineage
for i in `cat ABDfiles`
do
cat /data1/home/yafei/Project3/Pi/A/by_chr/$i  /data1/home/yafei/Project3/Pi/B/by_chr/$i  /data1/home/yafei/Project3/Pi/D/by_chr/$i > /data1/home/yafei/Project3/Pi/ABD/$i
done

#重新计算新疆的PI值
#目录xuebo@204:/data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf  pop19.txt

nohup vcf-concat chr01.vcf chr02.vcf chr07.vcf chr08.vcf chr013.vcf chr014.vcf chr019.vcf chr020.vcf chr025.vcf chr026.vcf chr031.vcf chr032.vcf chr037.vcf chr038.vcf > Alineage.vcf &
  nohup vcf-concat chr03.vcf chr04.vcf chr09.vcf chr010.vcf chr015.vcf chr016.vcf chr021.vcf chr022.vcf chr027.vcf chr028.vcf chr033.vcf chr034.vcf chr039.vcf chr040.vcf > Blineage.vcf &
  nohup vcf-concat chr05.vcf chr06.vcf chr011.vcf chr012.vcf chr017.vcf chr018.vcf chr023.vcf chr024.vcf chr029.vcf chr030.vcf chr035.vcf chr036.vcf chr041.vcf chr042.vcf > Dlineage.vcf &
  bgzip -c ${i} > ${i}.gz
nohup tabix -p vcf Alineage.vcf.gz &
  nohup tabix -p vcf Blineage.vcf.gz &
  nohup tabix -p vcf Dlineage.vcf.gz &
  #分染色体计算PI
  for chr in {1..42}
do
nohup vcftools --vcf /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/chr0${chr}.vcf --keep /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/pop19.txt --site-pi --out /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/xinjiang/${chr}_pop191M.group.sites.pi > nohupVmap 2>& 1 &
  done
#转移到yafei204对应的位置:/data1/home/yafei/Project3/Pi/xinjiang

for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/home/yafei/Project3/Pi/xinjiang/${chr}_pop191M.group.sites.pi.sites.pi --file2 /data1/home/yafei/Project3/Pi/SNPbed_1M/chr${chr}_1M.bed --out Awindow/${chr}_pop191M.group.sites.pi &
  done


for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/home/yafei/Project3/Pi/xinjiang/${chr}_pop191M.group.sites.pi.sites.pi --file2 /data1/home/yafei/Project3/Pi/SNPbed_1M/chr${chr}_1M.bed --out Bwindow/${chr}_pop191M.group.sites.pi &
  done

for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/home/yafei/Project3/Pi/xinjiang/${chr}_pop191M.group.sites.pi.sites.pi --file2 /data1/home/yafei/Project3/Pi/SNPbed_1M/chr${chr}_1M.bed --out Dwindow/${chr}_pop191M.group.sites.pi &
  done

#把Windows pi合并成染色体

for i in `cat D.name`
do
cat D/1M_window/ *_${i} >> D/by_chr/${i}
done

#ggplot2画图
#改名字
library(ggplot2)
nameA <- c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic")
nameB <- c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic")
nameD <- c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic","DD_Strangulata", "DD_Meyeri", "DD_Anathera")
nameAB <- c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic")
nameABD <- c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic")

#修改 路径和name以及输出文件名
path <- "/data1/home/yafei/Project3/Pi/D/by_chr"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 

length(data)

#A: 4, 13
#B: 2, 11
#D: 11
#AB: 10
#ABD:
for(i in 1:length(data)){
  data[[i]]$Name <- nameD[i]
  if (i<11) {
    data[[i]]$Fill <- "pink"
  } else if (i<20) {
    data[[i]]$Fill <- "lightblue"
  } else {
    data[[i]]$Fill <- "orange"
  }
}
all <- data.frame()
for(i in 1:length(data)){
  all <- rbind(all, data[[i]])
}
head(all)
tail(all)

#绘图：修改横坐标顺序和输出文件名
library(ggplot2)
P2<- ggplot(all, aes(x=Name, y=V4,fill=Fill),na.rm=TRUE) + 
  #geom_violin() + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  scale_fill_manual(values = c("#6000b4", "#1f97ff", "#ffce6e"))+
  #Alineage
  #scale_x_discrete(limits= c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Blineage
  #scale_x_discrete(limits= c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Dlineage
  scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera"))+
  #ABlineage
  #scale_x_discrete(limits= c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #ABDlineage
  #scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  theme_bw()+ #背景变为白色
  #B:0.01 D:0.004 AB:0.01
  ylab("Pi")+xlab("Lineages")+scale_y_continuous(limits=c(0,0.004)) + 
  theme(axis.text.x=element_text(angle=80,hjust = 1,colour="black",family="Times",size=200),
        axis.text.y=element_text(family="Times",size=300,face="plain"),
        axis.title.y=element_text(family="Times",size = 300,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())

pdf(file = "D_Pi.pdf",width=200,height=100) #结果保存
print(P2)
dev.off()


#-------------------------------------------------------------------shell--------------------------------------------------------------------------------------

#Linux Shell中使用awk完成两个文件的关联Join
#外关联
awk 'NR==FNR{a[$2]=$0;}NR!=FNR{print $0,a[$2]}' b.txt a.txt
#内关联
#方法1
awk -F',' 'NR==FNR{a[$1]=$2;}NR!=FNR && a[$1] {print $0,a[$1]}' b.txt a.txt
#方法2
awk -F',' 'NR==FNR{a[$1]=$2;}NR!=FNR && $1 in a {print $0,a[$1]}' b.txt a.txt

awk 'FNR==NR{a[$1];next}($1 in a){next} {print}' a b 

#-------------------------------------------------------------------CS IBS--------------------------------------------------------------------------------------
#提取CS的pos和posAllele文件做hapscanner。工作目录：203:/data1/home/yafei/Project3/Vmap1.1
for chr in {1..42}
do
zcat chr${chr}.all.vcf.gz |awk '{if(NF>100) print $1"\t"$2}' |sed '1d' > CS_pos/chr${chr}_pos.txt &
  zcat chr${chr}.all.vcf.gz |awk '{if(NF>100) print $1"\t"$2"\t"$4"\t"$5}' |sed '1c Chr\tPos\tRef\tAlt' > CS_posAllele/chr${chr}_posAllele.txt  &
  done

#生成参数文件 工作目录：204:/data2/xuebo/Projects/Speciation/More_accessions/hapScan/
nohup java -jar Yafei_Guo.jar &
  
  #进行CS的hapscanner扫描
  nohup bash bashCS.sh &
  #生成的结果文件 工作目录：203:/data2/yafei/Project3/CS_Vmap/CS/
  #合并lineage文件
  nohup vcf-concat chr01.vcf chr02.vcf chr07.vcf chr08.vcf chr013.vcf chr014.vcf chr019.vcf chr020.vcf chr025.vcf chr026.vcf chr031.vcf chr032.vcf chr037.vcf chr038.vcf | bgzip -c > AlineageCS.vcf.gz &
  nohup vcf-concat chr03.vcf chr04.vcf chr09.vcf chr010.vcf chr015.vcf chr016.vcf chr021.vcf chr022.vcf chr027.vcf chr028.vcf chr033.vcf chr034.vcf chr039.vcf chr040.vcf | bgzip -c > BlineageCS.vcf.gz &
  nohup vcf-concat chr05.vcf chr06.vcf chr011.vcf chr012.vcf chr017.vcf chr018.vcf chr023.vcf chr024.vcf chr029.vcf chr030.vcf chr035.vcf chr036.vcf chr041.vcf chr042.vcf | bgzip -c > DlineageCs.vcf.gz &
  nohup bcftools merge /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz /data2/yafei/Project3/CS_Vmap/CS/AlineageCS.vcf.gz -o Alineage_withCS.vcf &
  nohup bcftools merge /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz /data2/yafei/Project3/CS_Vmap/CS/BlineageCS.vcf.gz -o Blineage_withCS.vcf &
  nohup bcftools merge /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz /data2/yafei/Project3/CS_Vmap/CS/DlineageCS.vcf.gz -o Dlineage_withCS.vcf &
  
  #提取barely的vcf
  这个文件夹是分染色体的: /data2/xuebo/Projects/Speciation/tree/withBarley_segregate
这个是合起来的: /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage

#计算IBS:barely
nohup java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage/Alineage_withBarley.vcf.gz --out Alineage_withBarley.all.ibs.txt > logA.txt &
  nohup java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage/Blineage_withBarley.vcf.gz --out Blineage_withBarley.all.ibs.txt > logB.txt &
  nohup java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage/Dlineage_withBarley.vcf.gz --out Dlineage_withBarley.all.ibs.txt > logD.txt &
  
  #计算IBS:CS
  nohup java -jar -Xmx100g /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 Alineage_withCS.vcf.gz --out Alineage_withCS.all.ibs.txt > logA.txt &
  nohup java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 Blineage_withCS.vcf.gz --out Blineage_withCS.all.ibs.txt > logB.txt &
  nohup java -jar -Xmx100g /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 Dlineage_withCS.vcf.gz --out Dlineage_withCS.all.ibs.txt > logD.txt &
  
  #IBS画图 分不同的lineage
  #工作目录：yafei@203:/data2/yafei/Project3/CS_Vmap
  #读取group名的文件
  path <- "/data2/yafei/Project3/group/subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)

#修改 输入文件名
A <- read.table("Alineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
A$Lineage <- NA
A$Fill <- NA
nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_Durum","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_Khorasan_wheat","AABBDD_Landrace","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_Polish_wheat","AABB_Rivet_wheat","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
#A:487  485-487
#B:406  404-406
#D:306  304-306
for(i in 1:length(data)){
  A[which(A$Dxy %in% data[[i]][,1]),486] <- nameABD[i]
  if(i %in% c(3,4,9,12,13,14,19,22,23,25,28,29)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "pink"
  }else if(i %in% c(5,24,26)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "lightblue"
  }else if(i %in% c(6,7,8,10,11,16,17,18,27)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "orange"
  }else if(i %in% c(1,20)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "red"
  }else if(i %in% c(2,15,21)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "yellow"
  }
}
sub <- A[,485:487]

library(ggplot2)
P2<- ggplot(sub, aes(x=Lineage, y=CS, fill=Fill),na.rm=TRUE) + 
  #geom_violin() + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  scale_fill_manual(values = c("#6000b4", "#1f97ff", "#ffce6e"))+
  #Alineage
  #scale_x_discrete(limits= c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Blineage
  scale_x_discrete(limits= c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Dlineage
  #scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera"))+
  #ABlineage
  #scale_x_discrete(limits= c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #ABDlineage
  #scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  theme_bw()+ #背景变为白色
  #A:0.4 B:0.4 D:0.7
  ylab("IBS_CS")+xlab("Lineages")+scale_y_continuous(limits=c(0,0.4)) + 
  theme(axis.text.x=element_text(angle=80, hjust = 1, colour="black", family="Times", size=200),
        axis.text.y=element_text(family="Times",size=300,face="plain"),
        axis.title.y=element_text(family="Times",size = 300,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())

pdf(file = "B_IBS_CS.pdf",width=200,height=100) #结果保存
print(P2)
dev.off()

#画热图文件的产生
library("corrplot")
library("reshape")
path <- "/data2/yafei/Project3/group/subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)
D <- read.table("Blineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
rownames(D) <- D$Dxy
D<- D[,-1]
D<- as.matrix(D)
all <- as.data.frame(melt(D))
#colnames(all) <- c("name1","name2","value")
all$X1 <- as.character(all$X1)
all$X2 <- as.character(all$X2)
nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_Durum","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_Khorasan_wheat","AABBDD_Landrace","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_Polish_wheat","AABB_Rivet_wheat","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
for(i in 1:length(data)){
  all[which(all$X1 %in% data[[i]][,1]),1] <- nameABD[i]
}
for(i in 1:length(data)){
  all[which(all$X2 %in% data[[i]][,1]),2] <- nameABD[i]
}
colnames(all) <- c("id1","id2","value")
all$variable <- "var"
sub <- all[-which(all$value == 0),]

casted <- cast(sub,id1+id2~variable,mean)
ibs <- cast(casted, id1~id2)
ibs[is.na(ibs)] <- 0
names <- ibs$id1
write.table(ibs,"AB_ibs.txt",quote=F) 

#热图  
library("corrplot")
ibs <- read.table("D_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

#A:0.5 B:0.5 D:0.7
pdf("IBS_D_CS_heat.pdf",width=100,height=100)
#corrplot(data,method = "color",col = col3(10),tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=2,tl.cex=2,cl.cex=2,cl.lim = c(0, 1))
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=6,number.digits=3,tl.cex=10,cl.cex=12,cl.lim = c(0, 0.7))
dev.off()


#-------------------------------------------------------------------shell--------------------------------------------------------------------------------------
#提取vcf的特定区域
#yafei@203:/data2/yafei/Project3/make_tree
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
bedtools intersect -a /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr${chr}.withBarley.vcf.gz -b merge_A.bed -header > chr${chr}.withBarley.vcf &
  done

for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
bedtools intersect -a /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr${chr}.withBarley.vcf.gz -b merge_B.bed -header > chr${chr}.withBarley.vcf &
  done

for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
bedtools intersect -a /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr${chr}.withBarley.vcf.gz -b merge_D.bed -header > chr${chr}.withBarley.vcf &
  done

for i in `ls *vcf`
do 
bgzip -c ${i} > ${i}.gz &
  done

#-------------------------------------------------------------------计算reads depth.sh--------------------------------------------------------------------------------------
#工作目录：204:/data2/xuebo/Projects/Speciation/More_accessions/hapScan/
#位点depth
samtools depth -b bed_file sample.bam > sample.depth  #bed 用来指定统计区间，运行后输出指定区间每一个碱基的测序深度。
#区间depth
bedtools makewindows -g genome.txt -w 10000000 -s 2000000 > windows.bed
#bedtools makewindows用来自动生成划窗区间。-g genome.txt是要划分的基因组，格式为两列：染色体、染色体长度；-w 10000000为窗口大小为10M；-s 1000000为步长为1M，即窗口在染色体上每次向右平移1M的距离；windows.bed为输出文件，格式为三列：染色体、区间开始位点、区间结束位点。
bedtools coverage -a windows.bed -b xxx.sort.bam > xxx.depth.txt      
#bedtools coverage对划分好的每个滑动窗口进行reads数（depth)的统计。-a windows为上一步划分好的区间；-b xxx.sort.bam为测序数据mapping到参考基因组的比对文件；xxx.depth.txt为统计结果的输出文件，格式为7列：染色体、区间起始位点、区间结束位点、该区间内的reads数、该区间内的碱基数、区间大小、该区间的平均覆盖度。       
samtools bedcov bed_file samplename.bam > sample.bedcov
#输出的文件中计算了bed文件每一个区间的碱基总数，这里并不是reads的条数

#计算：Homozygous alternative allele
java -jar Yafei_Guo.jar(IBS)

#计算个体的杂合率和缺失率
#AABB
for i in {2..149}
do
sed -n "$i~149p" all_taxaQCfile.txt | awk '{sum1+=$2;sum2+=$3} END{print $i"\t"sum1/28"\t"sum2/28}'; done >> taxaQCfile.txt
#AABBDD
for i in {2..245}
do
sed -n "$i~245p" all_taxaQCfile.txt | awk '{sum1+=$2;sum2+=$3} END{print $i"\t"sum1/42"\t"sum2/42}'; done >> taxaQCfile.txt
#AA
for i in {2..92}
do
sed -n "$i~92p" all_taxaQCfile.txt | awk '{sum1+=$2;sum2+=$3} END{print $i"\t"sum1/14"\t"sum2/14}'; done >> taxaQCfile.txt
#BB
for i in {2..11}
do
sed -n "$i~11p" all_taxaQCfile.txt | awk '{sum1+=$2;sum2+=$3} END{print $i"\t"sum1/14"\t"sum2/14}'; done >> taxaQCfile.txt
#DD
for i in {2..59}
do
sed -n "$i~59p" all_taxaQCfile.txt | awk '{sum1+=$2;sum2+=$3} END{print $i"\t"sum1/14"\t"sum2/14}'; done >> taxaQCfile.txt

#计算平均Pi值
#工作目录：yafei@204:/data1/home/yafei/Project3/Pi/AB/by_chr
for i in `ls`; do awk '{print $4*1000000000}' $i|awk '{sum+=$1} END {printf("%.6f\n",sum/4942000000000)}'; done
#计算标准误.R
path <- "/data1/home/yafei/Project3/Pi/D/by_chr"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
all <- vector()
std <- function(x) sd(x,na.rm=T)/sqrt(length(x))
for(i in 1:length(data)){
  all[i] <-  std(data[[i]][,4])
}

#计算平均TajimaD值以及标准误
#工作目录：yafei@203:/data2/yafei/Project3/TajimasD/VmapData/A
path <- "/data2/yafei/Project3/TajimasD/VmapData/D"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 
SE <- vector()
Mean <- vector()
std <- function(x) sd(x,na.rm=T)/sqrt(length(x))
for(i in 1:length(data)){
  Mean[i] <- mean(data[[i]][,4],na.rm=T)
  SE[i] <-  std(data[[i]][,4])
}

cat species_2675.fasta |xargs |sed 's/>/\n>/g'  |sed 's/ /\t/1' |sed 's/ //g'  |awk '{print NFR, $2}' |less -LS

awk 'ORS=NR%2?" ":"\n"{print}' 

cat species_1.fasta species_2.fasta species_3.fasta species_4.fasta species_5.fasta species_6.fasta species_7.fasta species_8.fasta species_9.fasta species_10.fasta |awk 'ORS=NR%20000?" ":"\n"{print}' |sed 's/>/\n>/g' |sed 's/ /\t/' |sed 's/ //g'| sort |sed 's/>//'|sed 's/LandraceA/LandrA/' |sed 's/LandraceB/LandrB/' |sed 's/LandraceD/LandrD/'|sed 's/SpeltoidB/SpeltB/'|sed 's/Wildeink/Wildek/' |sed 's/WildemmerA/EmmerA/'|sed 's/WildemmerB/EmmerB/' |less -LS 

for i in {0..9}; do  for j in {1..13}; do echo $i; done; done |xargs -n 1 > row_num
paste row_num temp0 -d " " | sort -k2,2 -k1n,1 | awk '{$1=null;print $0}' |sed 's/ //' |awk '{if(FNR%10==1) {print "1\t"$1"\n"$0} else {print $0}}' > temp1

for i in {0..12}; do  for j in {1..10}; do echo $i; done; done |xargs -n 1 > row_num2
awk '{if(FNR%10==1) {print "\n"$0} else {print $0}}' row_num2 > row_num3
paste row_num3 temp1 -d "" |
  
  #批量杀死程序
  ps aux|grep Volca|tail -n 20 | awk '{print $2}' > id
for i in `cat id`; do kill -9 $i; done

#冰川数据
setwd("/Users/guoyafei/Downloads/stack/")
data <- read.table("LR04stack.txt",header=T,stringsAsFactors = F,sep="\t")
ggplot(data, aes(x=data$Time..ka., y=data$Benthic.d18O..per.mil.)) + 
  geom_line() +
  scale_x_log10()+scale_y_reverse() +
  geom_hline(aes(yintercept=4), colour="#990000", linetype="dashed")

#######################################################################  bwa  ##########################################################################

zcat file.fq.gz | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > file_sorted.fastq

for ID in {"CRR061683","CRR061684"}
do
ref=/data1/home/xinyue/ref/abd_iwgscV1.fa.gz
in1=/data1/home/xinyue/data/cleandata/${ID}_f1_paired.fq.gz
in2=/data1/home/xinyue/data/cleandata/${ID}_r2_paired.fq.gz
in="$in1 $in2"
out=/data1/home/xinyue/data/bamdata
echo -e "#!/bin/bash \n
bwa mem -t 32 -R '@RG\tID:$ID\tSM:$ID\tLB:SL\tPL:IL' $ref $in | samtools view -S -b - > $out/$ID.bam " > NGS_bwa.$ID.sh
sh NGS_bwa.$ID.sh
done
## sort
samtools sort -@ 2 -O bam -o $out/$ID.sorted.bam $out/$ID.bam
## MarkDuplicates and index
java -jar /data1/home/xuebo/software/picard.jar MarkDuplicates INPUT=$out/$ID.sorted.bam OUTPUT=$out/$ID.rm.bam METRICS_FILE=$out/$ID.metrics.txt
## index
samtools index $out/$ID.rm.bam

## split by chromosome
for chr in {0..42}
do
samtools view -b $out/$ID.rm.bam $chr > $out/$ID.chr$chr.bam
done
index bam file
for chr in {0..42}
do
samtools index $out/$ID.chr$chr.bam
done

samtools view test.bam | awk '{print $5"\t"$9}'| sed '1i mapping-qulity\tmate-length' > test
awk '{print $5"\t"$9}' standard.bam| sed '1i mapping-qulity\tmate-length' > standard

cat 687_f1_test.fq | paste - - - - | sort -k1,1 -t " "  > 687_f1_test.sorted
cat 688_f1_test.fq | paste - - - - | sort -k1,1 -t " "  > 688_f1_test.sorted
cat 691_f1_test.fq | paste - - - - | sort -k1,1 -t " "  > 691_f1_test.sorted

cat 687_r2_test.fq | paste - - - - | sort -k1,1 -t " "  > 687_r2_test.sorted
cat 688_r2_test.fq | paste - - - - | sort -k1,1 -t " "  > 688_r2_test.sorted
cat 691_r2_test.fq | paste - - - - | sort -k1,1 -t " "  > 691_r2_test.sorted

#cat sub_f1_687.fq | paste - - - - | sort -k1,1 -t " "  > 687_f1_test.sorted
#cat sub_r2_687.fq | paste - - - - | sort -k1,1 -t " "  > 687_r2_test.sorted
#grep "@" sub_f1_687.fq > 687_f1.id
#grep "@" sub_r2_687.fq > 687_r2.id
#cat 687_f1.id 687_r2.id |sort |uniq -c |grep "1 " | awk '{print $2}' > 687_uniqid

grep "@" 687_f1_test.fq > 687_f1.id
grep "@" 687_r2_test.fq > 687_r2.id
cat 687_f1.id 687_r2.id |sort |uniq -c |grep "1 " | awk '{print $2}' > 687_uniqid
grep "@" 688_f1_test.fq > 688_f1.id
grep "@" 688_r2_test.fq > 688_r2.id
cat 688_f1.id 688_r2.id |sort |uniq -c |grep "1 " | awk '{print $2}' > 688_uniqid
grep "@" 691_f1_test.fq > 691_f1.id
grep "@" 691_r2_test.fq > 691_r2.id
cat 691_f1.id 691_r2.id |sort |uniq -c |grep "1 " | awk '{print $2}' > 691_uniqid

rm *.id

grep -v -f 687_uniqid 687_f1_test.sorted | tr "\t" "\n" > 687_f1_test.sorted.fq
grep -v -f 688_uniqid 688_f1_test.sorted | tr "\t" "\n" > 688_f1_test.sorted.fq
grep -v -f 691_uniqid 691_f1_test.sorted | tr "\t" "\n" > 691_f1_test.sorted.fq

grep -v -f 687_uniqid 687_r2_test.sorted | tr "\t" "\n" > 687_r2_test.sorted.fq
grep -v -f 688_uniqid 688_r2_test.sorted | tr "\t" "\n" > 688_r2_test.sorted.fq
grep -v -f 691_uniqid 691_r2_test.sorted | tr "\t" "\n" > 691_r2_test.sorted.fq

bwa mem -t 8 /data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz 687_f1_test.sorted.fq  687_r2_test.sorted.fq  | samtools view -S -b -> 687_test.bam
bwa mem -t 8 /data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz 688_f1_test.sorted.fq  688_r2_test.sorted.fq  | samtools view -S -b -> 688_test.bam
bwa mem -t 8 /data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz 691_f1_test.sorted.fq  691_r2_test.sorted.fq  | samtools view -S -b -> 691_test.bam
bwa mem -t 8 /data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz 704_f1_test.fq  704_r2_test.fq  | samtools view -S -b -> 704_test.bam
bwa mem -t 8 /data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz 705_f1_test.fq  705_r2_test.fq  | samtools view -S -b -> 705_test.bam
bwa mem -t 8 /data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz 706_f1_test.fq  706_r2_test.fq  | samtools view -S -b -> 706_test.bam

samtools view 687_test.bam | awk '{print $5"\t"$9"\t687_test"}' | tail -n 2500 | sed '1i mapping-qulity\tmate-length\tfile' > 687_test

tail -n 5000 687_test > 687_test2
samtools view 688_test.bam | awk '{print $5"\t"$9"\t688_test"}' > 688_test
samtools view 691_test.bam | awk '{print $5"\t"$9"\t691_test"}' > 691_test
samtools view 704_test.bam | awk '{print $5"\t"$9"\t704_test"}' > 704_test
samtools view 705_test.bam | awk '{print $5"\t"$9"\t705_test"}'> 705_test
samtools view 706_test.bam | awk '{print $5"\t"$9"\t706_test"}' > 706_test

cat 687_test2  688_test 691_test 704_test 705_test 706_test| sed '1i mapping-qulity\tmate-length\tfile' > test

data <- read.table("/Users/guoyafei/Documents/Lulab/Project-4-VmapIII/test",header=T,stringsAsFactors = F)
data2 <- data[-which(data$mate.length > 1000),]
data3 <- data2[-which(data2$mate.length < -1000),]
ggplot(data3,aes(x=mapping.qulity) )+  geom_histogram(binwidth = 1) +facet_grid(file~.)+theme_classic()

nohup zcat CRR061687_r2.filtered.fq.gz \
| paste - - - - \
| sort -k1,1 -S 500G \
| tr '\t' '\n' \
| gzip > CRR061687_r2.filtered_sorted.fq.gz &
  
nohup seqkit sort -n CRR061687_r2.filtered.fq.gz | gzip -c > CRR061687_r2.filtered_sorted.fq.gz &
  
bcftools merge chr0${chr}.vcf.gz chr0${chr}.vcf.gz chr0${chr}.vcf.gz chr${chr}.vcf.gz -o chr${chr}.all.vcf

bcftools filter 1000Genomes.vcf.gz --regions 9:4700000-4800000 > 4700000-4800000.vcf


#E6:xuebo@204:/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225

#!/bin/bash
sleep 4800s
for i in {1..42}
do
  WGS --model vcf --type toXPCLR --group1 ../groupTaxa/EA.txt --rec ../../slidewindow_recomrate_updown_20M.txt --chr $i --file /data2/xuebo/Projects/Speciation/E6/chr${i}.Land.vcf.gz --out groupEAChr${i} &
done

#!/bin/bash
for i in {1..42}
do
	XPCLR -xpclr ../groupEA/groupEAChr${i}.geno ../groupWA/groupWAChr${i}.geno ../groupWA/groupWAChr${i}.snp EA_WA_10kchr${i} -w1 0.005 500 10000 $i -p1 0.95 &
done


#做迁徙和环境适应性路径
204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6
203:yafei:/data1/home/yafei/008_Software/snpEff/data2

/data2/xuebo/Projects/Speciation/xpclr/Selection_V2

#统计VCF文件的snp的密度
#hg19.bed
chr1 248956422
chr2 242193529
chr3 198295559
....

bedtools makewindows -g hg19.bed -w 1000000 > windows.bed
bedtools coverage -a windows.bed -b test.vcf -counts > coverage.txt


