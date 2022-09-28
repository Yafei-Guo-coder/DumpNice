# move the gvcf and tbi file into the directory 
gatk=/mnt/zhou/hangyuan/miniconda3/bin/gatk
reference=/mnt/zhou/hangyuan/BSA-seq/P20210613_AtBSA-raw.bam.results_20210621/01.rawData/TAIR10.fa
time $gatk CombineGVCFs \
      -R $reference \
      -V HR.HC.gvcf.gz \
      -V NHR.HC.gvcf.gz \
      -O BSA_with2pools.gvcf.gz
time $gatk GenotypeGVCFs \
     -R  $reference \
     -V "BSA_with2pools.gvcf.gz" \
     -O "BSA_no_filter.HC.vcf" 