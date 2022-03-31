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
