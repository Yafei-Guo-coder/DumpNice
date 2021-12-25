Working Directory:
xuebo@204:/data2/xuebo/Projects/Speciation/Variant_Age/Chr
xuebo@204:/data2/xuebo/Projects/Speciation/Variant_Age/Lineage
xuebo@204:/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225/Landrace_225_noAM_220_maf001/LD/chr6.lineage_Landrace_225_noAM_220_maf001_LD.vcf

for i in {1..42}
do
bash 42-21.sh /data2/xuebo/Projects/Speciation/E6/Landrace_locate_225/Landrace_225_noAM_220_maf001/LD/chr${i}.lineage_Landrace_225_noAM_220_maf001_LD.vcf > chr${i}.vcf &
done
vcf-concat chr1.vcf chr2.vcf | bgzip -c > chr1A.vcf.gz

#step1 转格式(vcf文件染色体号不能有字母，重组率文件必须有chr)
#yafei@204:/data1/home/yafei/Test
plink2 --vcf A_position.vcf.recode.vcf --allow-extra-chr --alt1-allele 'force' A_Alt.Anc.vcf.gz 4 3 '#' --export vcf --out A_Anc --autosome-num 42
plink2 --vcf B_position.vcf.recode.vcf --allow-extra-chr --alt1-allele 'force' B_Alt.Anc.vcf.gz 4 3 '#' --export vcf --out B_Anc --autosome-num 42
plink2 --vcf D_position.vcf.recode.vcf --allow-extra-chr --alt1-allele 'force' D_Alt.Anc.vcf.gz 4 3 '#' --export vcf --out D_Anc --autosome-num 42
#xuebo@204:/data2/xuebo/Projects/Speciation/Variant_Age/Chr

for i in {1..7}
do
for j in {"A","B","D"}
do
geva_v1beta --vcf ../Chr/chr${i}${j}.vcf.gz --map recombination_rate_geva2.txt --out chr${i}${j}_rec
done
done

for i in `ls *marker.txt`
do
awk '{print $3}' ${i} | sed '1,2d' |sed '$d' > ${i::-11}.pos &
done

#step2 运行
for i in {1..7}
do
for j in {"A","B","D"}
do
geva_v1beta -i chr${i}${j}_rec.bin -o ${i}${j}_RUN1 --positions chr${i}${j}_rec.pos --Ne 70000 --mut 6.5e-9 --hmm /data1/home/xuebo/software/geva-master/hmm/hmm_initial_probs.txt /data1/home/xuebo/software/geva-master/hmm/hmm_emission_probs.txt
done
done

for i in {1..7}
do
for j in {"A","B","D"}
do
Rscript /data1/home/xuebo/software/geva-master/estimate.R ${i}${j}_RUN1.pairs.txt 70000
done
done

---------------------------------------------本地画图---------------------------------------------
xuebo@204:/data2/xuebo/Projects/Speciation/Variant_Age/rec

setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/06_Allele_Age")

for i in `ls *rec.marker.txt`; do sed '2d' $i | sed '$d' > ${i::-11}.txt; done
for i in `ls *RUN1.sites2.txt`; do grep "J" $i | sed '1i MarkerID\tClock\tN_Concordant\tN_Discordant\tPostMode' > chr${i::-11}.txt; done
