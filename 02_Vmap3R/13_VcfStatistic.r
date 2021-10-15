#过滤maf提取样本
#AA
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/AA_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/AA/chr${chr}.mafAA.vcf.gz 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W_E2/E2_maf/AA/chr${chr}.mafAA.vcf.gz --out /data2/yafei/Project3/V+W_E2/E2_maf/AA/chr${chr}_siteQCfile.txt --out2 /data2/yafei/Project3/V+W_E2/E2_maf/AA/chr${chr}_taxaQCfile.txt > nohupA 2>& 1 
done
#BB
for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/chr${chr}.all.vcf.gz --keep /data2/yafei/Project3/group/SS_taxa.txt --max-missing 0.1 --maf 0.0000001 --recode --stdout | bgzip -c > /data1/home/yafei/Project3/Maf/mafBB/chr${chr}.mafBB.vcf.gz 
java -jar /data1/home/yafei/C36_checkQuality.jar --file /data1/home/yafei/Project3/Maf/mafBB/chr${chr}.mafBB.vcf.gz --out /data1/home/yafei/Project3/Maf/mafBB/chr${chr}_siteQCfile.txt --out2 /data1/home/yafei/Project3/Maf/mafBB/chr${chr}_taxaQCfile.txt > nohupB 2>& 1 
done
#DD
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/DD_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/DD/chr${chr}.mafDD.vcf.gz 
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

