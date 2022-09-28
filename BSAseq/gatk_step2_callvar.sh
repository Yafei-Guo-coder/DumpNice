#!/bin/bash
# this is a whole gonome gatk snp calling full process
#Uage: sh gatk_step2.sh sample_1.fastq sample_2.fastq out_dictionary
fastp=/usr/local/bin/fastp
bwa=/usr/bin/bwa
gatk=/mnt/zhou/hangyuan/miniconda3/bin/gatk
samtools=/usr/local/bin/samtools
reference=/mnt/zhou/hangyuan/BSA-seq/P20210613_AtBSA-raw.bam.results_20210621/01.rawData/TAIR10.fa
NTHREADS=30
fq1=$1
fq2=$2
sample=${fq1%%_*}
outdir=$3
outdir=${outdir}/${sample}
 
if [ ! -d "$outdir/cleanfq" ]
then mkdir -p "$outdir/cleanfq"
fi
 
if [ ! -d "$outdir/bwa" ]
then mkdir -p "$outdir/bwa"
fi
 
if [ ! -d "$outdir/gatk" ]
then mkdir -p "$outdir/gatk"
fi
 
time $fastp -i "$fq1" \
            -I "$fq2" \
            -o "$outdir/cleanfq/${sample}_1.fastq" \
            -O "$outdir/cleanfq/${sample}_2.fastq" \
            -h "$outdir/cleanfq/${sample}.html" \
            -j "$outdir/cleanfq/${sample}.json" \
                 && echo '** fq QC done'
time $bwa mem -M -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA" \
          -t $NTHREADS \
          $reference \
          "$outdir/cleanfq/${sample}_1.fastq"  "$outdir/cleanfq/${sample}_2.fastq" \
          | $samtools view -bS - > "$outdir/bwa/${sample}.bam"\
          && echo '** BWA MEM done ** '
time $samtools sort \
          -@ $NTHREADS \
          -O bam \
          -o "$outdir/bwa/${sample}.sorted.bam" \
          "$outdir/bwa/${sample}.bam"\
           && echo "** sorted raw bamfiles done"
time $samtools index "$outdir/bwa/${sample}.sorted.bam" && echo "** ${sample}.sorted.bam index done "
 
# 去除PCR重复
time $gatk MarkDuplicates \
           -I "$outdir/bwa/${sample}.sorted.bam" \
           -M "$outdir/bwa/${sample}.markup_metrics.txt" \
           -O "$outdir/bwa/${sample}.sorted.markup.bam" \
           1> "$outdir/bwa/${sample}.log.mark" 2>&1 \
           && echo "** ${sample}.sorted.bam MarkDuplicates done **"
 
# samtools index 去除PCR标记的bam 文件
time $samtools index "$outdir/bwa/${sample}.sorted.markup.bam" \
              && echo " ** ${sample}.sorted.markup.bam index done **"
 
#gatk开始：必选 -I -O -R，代表输入、输出、参考
#接下来可以按照字母顺序依次写出来，这样比较清晰
#-bamout：将一整套经过gatk程序重新组装的单倍体基因型（haplotypes）输出到文件
#-stand-call-conf :低于这个数字的变异位点被忽略，可以设成标准30（默认是10）
time $gatk --java-options "-Xmx100G -Djava.io.tmpdir=./tmp" HaplotypeCaller \
    -R $reference \
    -I "$outdir/bwa/${sample}.sorted.markup.bam" \
    -O "$outdir/gatk/${sample}.HC.gvcf.gz" \
    --tmp-dir "$outdir/gatk"\
    --emit-ref-confidence GVCF \
    -stand-call-conf 10 \
    && echo "** GVCF ${sample}.HC.gvcf.gz done ***"