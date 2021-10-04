#######################################################################  bwa  ##########################################################################

zcat file.fq.gz | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > file_sorted.fastq

for ID in {"CRR061683","CRR061684"}
do
ref=/data1/home/xinyue/ref/abd_iwgscV1.fa.gz
in1=/data1/home/xinyue/data/cleandata/${ID}_f1_paired.fq.gz
in2=/data1/home/xinyue/data/cleandata/${ID}_r2_paired.fq.gz
in="$in1 $in2"
out=/data1/home/xinyue/data/bamdata
echo -e "#!/bin/bash \n
bwa mem -t 32 -R '@RG\tID:$ID\tSM:$ID\tLB:SL\tPL:IL' $ref $in | samtools view -S -b - > $out/$ID.bam " > NGS_bwa.$ID.sh
sh NGS_bwa.$ID.sh
done
## sort
samtools sort -@ 2 -O bam -o $out/$ID.sorted.bam $out/$ID.bam
## MarkDuplicates and index
java -jar /data1/home/xuebo/software/picard.jar MarkDuplicates INPUT=$out/$ID.sorted.bam OUTPUT=$out/$ID.rm.bam METRICS_FILE=$out/$ID.metrics.txt
## index
samtools index $out/$ID.rm.bam

