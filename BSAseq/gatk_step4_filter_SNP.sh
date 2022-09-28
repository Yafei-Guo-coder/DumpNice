#!/bin/bash
gatk=/mnt/zhou/hangyuan/miniconda3/bin/gatk
reference=/mnt/zhou/hangyuan/BSA-seq/P20210613_AtBSA-raw.bam.results_20210621/01.rawData/TAIR10.fa

#select SNP：提取出SNP
time $gatk SelectVariants -R $reference \
    -V  BSA_no_filter.HC.vcf \
    -select-type SNP \
    -O BSA_combine_SNPs.vcf
#select INDELS：提取出INDELS
time $gatk SelectVariants -R $reference \
    -V  BSA_no_filter.HC.vcf \
    -select-type INDEL  \
    -O BSA_combine_INDELs.vcf

#To filter SNPs：过滤出高质SNPs
time $gatk VariantFiltration  \
    -V BSA_combine_SNPs.vcf \
    -filter "MQ < 40.0" \
    --filter-name "MQ_filter_SNP" \
    -filter "FS > 60.0" \
    --filter-name "FS_filter_SNP" \
    -filter  "QD < 4.0" \
    --filter-name "QD_filter_SNP" \
    -O BSA_SNPs_filter.vcf
 
grep -E '^#|PASS' BSA_SNPs_filter.vcf > BSA_SNPs_filterPASSED.vcf
 
#To filter INDELS：过滤出高质量INDELS
time $gatk VariantFiltration \
    -V BSA_combine_INDELs.vcf \
    -filter "QD < 4.0" \
    --filter-name "QD_filter_INDEL" \
    -filter "FS > 200.0" \
    --filter-name "FS_filter_INDEL" \
    -O BSA_INDELs_filter.vcf

grep -E '^#|PASS' BSA_INDELs_filter.vcf > BSA_INDELs_filterPASSED.vcf

#合并snp过滤结果和indel过滤结果
time $gatk MergeVcfs \
    -I BSA_INDELs_filterPASSED.vcf \
    -I BSA_SNPs_filterPASSED.vcf \
    -O BSA_filter_merged_SNP_INDEL.vcf

#This is the number of sites before filtering:
grep -vc "^#" BSA_combine_SNPs.vcf
grep -vc "^#" BSA_combine_INDELs.vcf
#This is the number of sites retained after filtering:
grep -vc "^#" BSA_SNPs_filter.vcf
grep -vc "^#" BSA_INDELs_filter.vcf
 
time $gatk VariantsToTable \
	 -R $reference \
	 -V BSA_filter_merged_SNP_INDEL.vcf \
	 -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -GF GQ -GF PL \
	 -O BSA.filter.table