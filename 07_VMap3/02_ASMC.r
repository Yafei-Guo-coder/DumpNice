#vcftools --gzvcf Alineage001_Wild_emmer.vcf.gz --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS'


for i in `cat population1.txt`
do
vcftools --gzvcf Alineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Alineage001_${i}.frq &
vcftools --gzvcf Blineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Blineage001_${i}.frq &
vcftools --gzvcf Dlineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Dlineage001_${i}.frq &
done

for i in `cat population2.txt`
do
vcftools --gzvcf Alineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Alineage001_${i}.frq &
vcftools --gzvcf Blineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Blineage001_${i}.frq &
done

for i in `cat population3.txt`
do
vcftools --gzvcf Dlineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Dlineage001_${i}.frq &
done