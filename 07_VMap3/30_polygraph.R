# admixtools
# 203或者66: admixtools conda 环境
# export LD_LIBRARY_PATH=/data1/home/yafei/008_Software/anaconda3/envs/admixtools/lib:$LD_LIBRARY_PATH

# Plink GWAS
cat shuf_prec8_*_ABlineage_mlm2.txt | awk '{print $7}' | grep -v p | grep -v NaN | sort -k1,1g > AB.prec8.thresh.txt &
cat shuf_prec8_*_Dlineage_mlm2.txt | awk '{print $7}' | grep -v p | grep -v NaN | sort -k1,1g > D.prec8.thresh.txt &
head -n 977319 AB.prec8.thresh.txt | tail
head -n 138722 D.prec8.thresh.txt | tail
cat prec8_ABlineage_mlm3.txt prec8_Dlineage_mlm3.txt > prec8_mlm3.txt

# Determine ancestry status
for i in {01..42}; do awk 'NR==FNR{a[$1]=$1;b[$1] = $2} NR!=FNR{if($2 in a && b[$2] == $5) print $1"\t"$4-1"\t"$4"\t"b[$2]; else if ($2 in a && b[$2] == $6) print $1"\t"$4-1"\t"$4"\t"b[$2]}' /data1/home/yafei/010_DataSet/Ancestor_state/chr${i}.txt /data4/home/yafei/plink_VCF/chr0${i}.bim > chr${i}.VMap3.bed & done

# Neutral input file
awk '{if($7 > 0.05267) print $3"\t"$4"\t"$2"\t"$7}' prec8_ABlineage_mlm2.txt > prec8_ABlineage_pneutral.txt
awk '{if($7 > 0.054) print $3"\t"$4"\t"$2"\t"$7}' prec8_Dlineage_mlm2.txt > prec8_Dlineage_pneutral.txt
cat prec8_ABlineage_pneutral.txt prec8_Dlineage_pneutral.txt > prec8_pneutral.txt
awk 'NR==FNR{a[$3]=$3;}NR!=FNR{if($2 in a) print $3"\t"$4-1"\t"$4"\t"$2"\t"$5"\t"$6}' prec8_pneutral.txt prec8_mlm3.txt |sort -k1,1n -k2,2n  |sed '1,2d' >  prec8_neutral.txt
for i in {01..42}
do
awk 'NR==FNR{a[$3]=$3;b[$3]=$5} NR!=FNR{if($3 in a && $5 != b[$3] && $5 ~ /(A|T|C|G)/) print $1"\t"$3"\t"$4"\t"$5"\t"$6}' /data1/home/yafei/010_DataSet/Ancestor_state/ancestor/chr${i}.VMap3.bed  prec8_neutral.txt | awk '{ if ( $2 != prev) { print $0; prev = $3 } }' >> prec8.part1_neutral.txt
done
awk '{print $3}' prec8.part1_neutral.txt > pos_neutral.txt
plink --vcf ABlineage_gwas.vcf.gz --freq  --extract pos_neutral.txt --within region3.txt --autosome-num 42 --out ABlineage_neutral
plink --vcf Dlineage_gwas.vcf.gz --freq  --extract pos_neutral.txt --within region3.txt --autosome-num 42 --out Dlineage_neutral

awk '{print $2"\t"$3"\t"$4"\t"$7","$8-$7}' ABlineage_neutral.frq.strat > AB.count_neutral.txt
awk '{print $2"\t"$3"\t"$4"\t"$7","$8-$7}' Dlineage_neutral.frq.strat > D.count_neutral.txt

# R
library(reshape2)
AB <- read.table("AB.count_neutral.txt", header=T,stringsAsFactors=F)
D <- read.table("D.count_neutral.txt", header=T,stringsAsFactors=F)
colnames(AB) <- c("SNP","variable","A1","value")
colnames(D) <- c("SNP","variable","A1","value")

AB_out <- dcast(AB, SNP+A1~variable)
D_out <- dcast(D, SNP+A1~variable)
write.table(AB_out,"AB_out_neutral.txt",quote=F,sep="\t",row.names=F)
write.table(D_out,"D_out_neutral.txt",quote=F,sep="\t",row.names=F)

# Shell
sed 's/-/\t/1' AB_out_neutral.txt | sed '1d' | sort -k1,1n -k2,2n > AB_out2_neutral.txt
sed 's/-/\t/1' D_out_neutral.txt | sed '1d' | sort -k1,1n -k2,2n > D_out2_neutral.txt
cat AB_out2_neutral.txt D_out2_neutral.txt | sort -k1,1n -k2,2n > prec8_allele_neutral.frq

sed -i 's/\t/-/1' prec8_allele_neutral.frq 

awk 'NR==FNR{a[$1]=$1;b[$1]=$0}NR!=FNR{if($3 in a) print $0,b[$3]}' prec8_allele_neutral.frq prec8.part1_neutral.txt | awk '{print $1"\t"$2"\tNA\tNA\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' | sed '1i CHROM\tPOS\tSNPID\tEFFECT\tAF\tAM\tEA\tEU\tIA\tSH\tWA' > neutral.polygraph.txt

# GWAS input file format
awk '{if($7 < 0.05267) print $3"\t"$4"\t"$2"\t"$7}' prec8_ABlineage_mlm2.txt > prec8_ABlineage_psig.txt
awk '{if($7 < 0.054) print $3"\t"$4"\t"$2"\t"$7}' prec8_Dlineage_mlm2.txt > prec8_Dlineage_psig.txt
cat prec8_ABlineage_psig.txt prec8_Dlineage_psig.txt > prec8_psig.txt
awk 'NR==FNR{a[$3]=$3;}NR!=FNR{if($2 in a) print $3"\t"$4-1"\t"$4"\t"$2"\t"$5"\t"$6}' prec8_psig.txt prec8_mlm3.txt | sort -k1,1n -k2,2n  |sed '1,2d' >  prec8_gwas.txt

for i in {01..42}
do
awk 'NR==FNR{a[$3]=$3;b[$3]=$5} NR!=FNR{if($3 in a && $5 != b[$3] && $5 ~ /(A|T|C|G)/) print $1"\t"$3"\t"$4"\t"$5"\t"$6}' /data1/home/yafei/010_DataSet/Ancestor_state/ancestor/chr${i}.VMap3.bed  prec8_gwas.txt | awk '{ if ( $2 != prev) { print $0; prev = $3 } }' >> prec8.part1_gwas.txt
done

awk '{print $3}' prec8.part1_gwas.txt > pos_gwas.txt
plink --vcf ABlineage_gwas.vcf.gz --freq  --extract pos_gwas.txt --within region3.txt --autosome-num 42 --out ABlineage
plink --vcf Dlineage_gwas.vcf.gz --freq  --extract pos_gwas.txt --within region3.txt --autosome-num 42 --out Dlineage

awk '{print $2"\t"$3"\t"$4"\t"$7","$8-$7}' ABlineage.frq.strat > AB.count.txt
awk '{print $2"\t"$3"\t"$4"\t"$7","$8-$7}' Dlineage.frq.strat > D.count.txt

# R
library(reshape2)
AB <- read.table("AB.count.txt", header=T,stringsAsFactors=F)
D <- read.table("D.count.txt", header=T,stringsAsFactors=F)
colnames(AB) <- c("SNP","variable","A1","value")
colnames(D) <- c("SNP","variable","A1","value")
AB_out <- dcast(AB, SNP+A1~variable)
D_out <- dcast(D, SNP+A1~variable)
write.table(AB_out,"AB_out.txt",quote=F,sep="\t",row.names=F)
write.table(D_out,"D_out.txt",quote=F,sep="\t",row.names=F)

# Shell
sed 's/-/\t/1' AB_out.txt | sed '1d' | sort -k1,1n -k2,2n > AB_out2.txt
sed 's/-/\t/1' D_out.txt | sed '1d' | sort -k1,1n -k2,2n > D_out2.txt
cat AB_out2.txt D_out2.txt | sort -k1,1n -k2,2n > prec8_allele.frq

sed -i 's/\t/-/1' prec8_allele.frq 

awk 'NR==FNR{a[$1]=$1;b[$1]=$0}NR!=FNR{if($3 in a) print $0,b[$3]}' prec8_allele.frq prec8.part1_gwas.txt | awk '{print $1"\t"$2"\tNA\tNA\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' | sed '1i CHROM\tPOS\tSNPID\tEFFECT\tAF\tAM\tEA\tEU\tIA\tSH\tWA' > gwas.polygraph.txt

# Graph parameter file

Rscript Plot_Trace.R trace_prec8.txt qfile_prec8.txt prec8 graph.R boxplot_prec8.pdf mcmc_prec8.pdf qb_statistics_prec8.pdf 0.075

aoyue/vmap2/daxing/analysis/008_vmap2_1062_spelt/006_introgression_ancestralSite/004_fdRes/003_byIndividual

cp -r /data4/home/aoyue/vmap2/genotype /mnt/usb2/204data/aoyue/vmap2
cp -r /data1/home/aoyue/biosoftware /mnt/usb2/204data/aoyue
cp -r /data1/home/aoyue/sift /mnt/usb2/204data/aoyue
