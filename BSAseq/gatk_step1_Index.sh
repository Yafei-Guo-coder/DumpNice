#!/bin/bash 
# building sequence alginment dictionary, samtools faidx and gatk creatSequenceDictionary
#Usage: sh gatk_step1.sh /path/your_genome.fasta 
bwa=/usr/bin/bwa                   # set where to find software 
gatk=/mnt/zhou/hangyuan/miniconda3/bin/gatk
samtools=/usr/local/bin/samtools
 
#bwa index
reference=$1
time $bwa index "$reference" && echo "** bwa index done! ** "
#samtools index
time  $samtools faidx $reference && echo "** samtools faidx done! ** "
 
#注意：使用GATK之前，需要先建立参考基因组索引文件.dict和.fai
#.dict中包含了基因组中contigs的名字，也就是一个字典；
#.fai也就是fasta index file，索引文件，可以快速找出参考基因组的碱基，由samtools faidx构建
#构建.dict文件（原来要使用picard的CreateSequenceDictionary模块，但是现在gatk整合了此模块，可以直接使用）
# gatk createSequenceDictionary
time $gatk --java-options "-Xmx100G -Djava.io.tmpdir=./tmp" CreateSequenceDictionary \
    -R "$reference" \
    -O "$reference.dict" \
    && echo "** gatk createSequenceDictionary done! **"