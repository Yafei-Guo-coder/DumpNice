#提取样本
for chr in {05,06,11,12,17,18,23,24,29,30,35,36,41,42}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/DD_VCF/chr0${chr}.vcf.gz --max-alleles 4 --keep /data2/yafei/004_Vmap3/VCF/taxa_DD.txt --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/finalDD_vcf/chr0${chr}.vcf.gz &
done b

#过滤mac以及多等位基因-AABB
for chr in {01,02,03,04,07,08,09,10,13,14,15,16,19,20,21,22,25,26,27,28,31,32,33,34,37,38,39,40}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/AABB_VCF/chr0${chr}.vcf.gz --min-alleles 2 --max-alleles 2 --mac 2 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/chr0${chr}.mac.bialle.vcf.gz &
done
#过滤mac以及多等位基因-DD
for chr in {05,06,11,12,17,18,23,24,29,30,35,36,41,42}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/finalDD_vcf/chr0${chr}.vcf.gz --min-alleles 2 --max-alleles 2 --mac 2 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/finalDD_vcf/chr0${chr}.mac.bialle.vcf.gz
done
#统计-AABB
java -jar /data2/yafei/004_Vmap3/VCF/CheckVcf.jar /data2/yafei/004_Vmap3/VCF/DD_VCF /data2/yafei/004_Vmap3/VCF/DD_VCF/checkVcf.txt

#拆分VCF为INdel和SNP-AABB
for chr in {01,02,03,04,07,08,09,10,13,14,15,16,19,20,21,22,25,26,27,28,31,32,33,34,37,38,39,40}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/filt_AABB/chr0${chr}.mac.bialle.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/Snp/chr0${chr}.snp.vcf.gz &
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/filt_AABB/chr0${chr}.mac.bialle.vcf.gz --keep-only-indels --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/Indel/chr0${chr}.indel.vcf.gz &
done
#拆分VCF为INdel和SNP-DD
for chr in {05,06,11,12,17,18,23,24,29,30,35,36,41,42}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/finalDD_vcf/chr0${chr}.mac.bialle.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/Snp/chr0${chr}.snp.vcf.gz &
vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/finalDD_vcf/chr0${chr}.mac.bialle.vcf.gz --keep-only-indels --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/Indel/chr0${chr}.indel.vcf.gz &
done

vcftools --gzvcf /data2/yafei/004_Vmap3/VCF/AABB_VCF/chr002.vcf.gz --min-alleles 2 --max-alleles 2 --mac 2 --recode --recode-INFO-all --stdout | bgzip -c > /data2/yafei/004_Vmap3/VCF/chr002.mac.bialle.vcf.gz &
