discoal 16 20 550 -Pt 175.204699 1752.046985 -Pre 1925.251684 5775.755051 -x 0.5

python diploSHIC.py fvecSim diploid out.txt out.diploid.fvec --totalPhysLen 55000 --maskFileName exampleApplication/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.accessible.fa.gz --chrArmsForMasking 3R

for f in exampleApplication/*.msOut.gz; do python2 diploSHIC.py fvecSim diploid $f $f.diploid.fvec --totalPhysLen 55000 --maskFileName exampleApplication/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.accessible.fa.gz --chrArmsForMasking 3R; done

vcf-concat *vcf.gz.snp.gz | bgzip -c > Alineage.vcf.gz
vcf-concat *vcf.gz.snp.gz | bgzip -c > Blineage.vcf.gz

mv chr001.mac.bialle.vcf.gz.snp.gz chr002.mac.bialle.vcf.gz.snp.gz chr007.mac.bialle.vcf.gz.snp.gz chr008.mac.bialle.vcf.gz.snp.gz chr013.mac.bialle.vcf.gz.snp.gz chr014.mac.bialle.vcf.gz.snp.gz chr019.mac.bialle.vcf.gz.snp.gz chr020.mac.bialle.vcf.gz.snp.gz chr025.mac.bialle.vcf.gz.snp.gz chr026.mac.bialle.vcf.gz.snp.gz chr031.mac.bialle.vcf.gz.snp.gz chr032.mac.bialle.vcf.gz.snp.gz chr037.mac.bialle.vcf.gz.snp.gz chr038.mac.bialle.vcf.gz.snp.gz Alineage/

/data2/yafei/004_Vmap3/VCF/filt_AABB/snp/Blineage/Blineage.vcf.gz
/data2/yafei/004_Vmap3/VCF/filt_AABB/snp/Alineage/Alineage.vcf.gz
/data2/yafei/004_Vmap3/VCF/finalDD_vcf/filt_DD/snp/Dlineage.vcf.gz

#IBS统计--正在进行
nohup java -jar -Xmx100g /data2/yafei/004_Vmap3/JarCode/C41_getIBS_distance2.jar --file1 /data2/yafei/004_Vmap3/VCF/finalDD_vcf/filt_DD/snp/Dlineage.vcf.gz --out Dlineage.ibs.txt > logD.txt &
java -jar -Xmx100g /data2/yafei/004_Vmap3/JarCode/C41_getIBS_distance2.jar --file1 /data2/yafei/004_Vmap3/VCF/filt_AABB/snp/Blineage/Blineage.vcf.gz --out Blineage.ibs.txt

#基本统计
java -jar /data2/yafei/004_Vmap3/JarCode/C36_checkQuality.jar --file /data2/yafei/004_Vmap3/VCF/finalDD_vcf/filt_DD/snp/Dlineage.vcf.gz --out Dlineage_siteQCfile.txt --out2 Dlineage_taxaQCfile.txt > nohupD 2>& 1 
java -jar /data2/yafei/004_Vmap3/JarCode/C36_checkQuality.jar --file /data2/yafei/004_Vmap3/VCF/filt_AABB/snp/Blineage/Blineage.vcf.gz --out Blineage_siteQCfile.txt --out2 Blineage_taxaQCfile.txt

/data2/yafei/004_Vmap3/VCF/TaxaGroup/Ploidy
AABB.txt  AABBDD.txt  DD.txt

#过滤maf提取样本
#DD
for chr in {005,006,011,012,017,018,023,024,029,030,035,036,041,042}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/finalDD_vcf/filt_DD/snp/chr${chr}.mac.bialle.vcf.gz.snp.gz --keep /data2/yafei/004_Vmap3/VCF/TaxaGroup/Ploidy/DD.txt --max-missing 0.8 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/MissVcf/chr${chr}.missDD.vcf.gz 
java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C36_checkQuality.jar --file /data2/yafei/004_Vmap3/VCF/MissVcf/chr${chr}.missDD.vcf.gz --out /data2/yafei/004_Vmap3/VCF/MissVcf/DD/chr${chr}_siteQCfile.txt --out2 /data2/yafei/004_Vmap3/VCF/MissVcf/DD/chr${chr}_taxaQCfile.txt > nohupD 2>& 1 
done

#AABBDD
for chr in {005,006,011,012,017,018,023,024,029,030,035,036,041,042,001,002,007,008,013,014,019,020,025,026,031,032,037,038,003,004,009,010,015,016,021,022,027,028,033,034,039,040}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/filt_AABB/snp/chr${chr}.mac.bialle.vcf.gz.snp.gz --keep /data2/yafei/004_Vmap3/VCF/TaxaGroup/Ploidy/AABBDD.txt --max-missing 0.8 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/MissVcf/chr${chr}.missAABBDD.vcf.gz 
java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C36_checkQuality.jar --file /data2/yafei/004_Vmap3/VCF/MissVcf/chr${chr}.missAABBDD.vcf.gz --out /data2/yafei/004_Vmap3/VCF/MissVcf/AABBDD/chr${chr}_siteQCfile.txt --out2 /data2/yafei/004_Vmap3/VCF/MissVcf/AABBDD/chr${chr}_taxaQCfile.txt > nohupABD 2>& 1 
done

#AABB
for chr in {001,002,007,008,013,014,019,020,025,026,031,032,037,038,003,004,009,010,015,016,021,022,027,028,033,034,039,040}
do
#vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/filt_AABB/snp/chr${chr}.mac.bialle.vcf.gz.snp.gz --keep /data2/yafei/004_Vmap3/VCF/TaxaGroup/Ploidy/AABB.txt --max-missing 0.8 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/MissVcf/chr${chr}.missAABB.vcf.gz 
java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C36_checkQuality.jar --file /data2/yafei/004_Vmap3/VCF/MissVcf/chr${chr}.missAABB.vcf.gz  --out /data2/yafei/004_Vmap3/VCF/MissVcf/AABB/chr${chr}_siteQCfile.txt --out2 /data2/yafei/004_Vmap3/VCF/MissVcf/AABB/chr${chr}_taxaQCfile.txt > nohupAB 2>& 1 
done

#AABBDD
for chr in {005,006,011,012,017,018,023,024,029,030,035,036,041,042,001,002,007,008,013,014,019,020,025,026,031,032,037,038,003,004,009,010,015,016,021,022,027,028,033,034,039,040}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/filt_AABB/snp/chr${chr}.mac.bialle.vcf.gz.snp.gz --keep /data2/yafei/004_Vmap3/VCF/TaxaGroup/Ploidy/AABBDD.txt --max-missing 0.8 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/MissVcf/chr${chr}.missAABBDD.vcf.gz 
java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C36_checkQuality.jar --file /data2/yafei/004_Vmap3/VCF/MissVcf/chr${chr}.missAABBDD.vcf.gz --out /data2/yafei/004_Vmap3/VCF/MissVcf/AABBDD/chr${chr}_siteQCfile.txt --out2 /data2/yafei/004_Vmap3/VCF/MissVcf/AABBDD/chr${chr}_taxaQCfile.txt > nohupABD 2>& 1 
done

#生成test的结果文件，chr001和chr003

java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C41_getIBS_distance2.jar --file1 /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001.mac.bialle.vcf.gz.snp.gz  --out /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001.ibs.txt
java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C41_getIBS_distance2.jar --file1 /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr003.mac.bialle.vcf.gz.snp.gz  --out /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr003.ibs.txt
java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C41_getIBS_distance2.jar --file1 /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr005.mac.bialle.vcf.gz.snp.gz  --out /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr005.ibs.txt
java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C36_checkQuality.jar --file /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001.missAABB.vcf.gz --out /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001_AABB_siteQCfile.txt --out2 /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001_AABB_taxaQCfile.txt > nohupAB 2>& 1 
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/filt_AABB/snp/chr001.mac.bialle.vcf.gz.snp.gz --keep /data2/yafei/004_Vmap3/VCF/TaxaGroup/Ploidy/AABBDD.txt --max-missing 0.8 --mac 2 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001.missAABBDD.vcf.gz
java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/C36_checkQuality.jar --file /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001.missAABBDD.vcf.gz --out /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001_AABBDD_siteQCfile.txt --out2 /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr001_AABBDD_taxaQCfile.txt > nohupABD 2>& 1 
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr005.missDD.vcf.gz --mac 2 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/Vmap3Test/02_Ploidy/chr005.mac.missDD.vcf.gz

sed '1d' chr001_AABB_siteQCfile.txt | shuf -n 50000 |sort -k2n,2 |sed '1i Chr\tPos\tHeterozygous\tProportion\tMissingRate\tMaf' > chr001_AABB_siteQC_50000.txt

java -jar -Xmx200g /data2/yafei/004_Vmap3/JarCode/Depth.jar
