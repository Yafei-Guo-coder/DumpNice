#-------------------------------------------------------------------shell--------------------------------------------------------------------------------------
#Linux Shell中使用awk完成两个文件的关联Join
#外关联
awk 'NR==FNR{a[$2]=$0;}NR!=FNR{print $0,a[$2]}' b.txt a.txt
#内关联
#方法1
awk -F',' 'NR==FNR{a[$1]=$2;}NR!=FNR && a[$1] {print $0,a[$1]}' b.txt a.txt
#方法2
awk -F',' 'NR==FNR{a[$1]=$2;}NR!=FNR && $1 in a {print $0,a[$1]}' b.txt a.txt

awk 'FNR==NR{a[$1];next}($1 in a){next} {print}' a b 


#-------------------------------------------------------------------shell--------------------------------------------------------------------------------------
#提取vcf的特定区域
#yafei@203:/data2/yafei/Project3/make_tree
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
bedtools intersect -a /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr${chr}.withBarley.vcf.gz -b merge_A.bed -header > chr${chr}.withBarley.vcf &
  done

for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
bedtools intersect -a /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr${chr}.withBarley.vcf.gz -b merge_B.bed -header > chr${chr}.withBarley.vcf &
  done

for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
bedtools intersect -a /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr${chr}.withBarley.vcf.gz -b merge_D.bed -header > chr${chr}.withBarley.vcf &
  done

for i in `ls *vcf`
do 
bgzip -c ${i} > ${i}.gz &
  done

awk 'ORS=NR%2?" ":"\n"{print}' 
  
#批量杀死程序
ps aux|grep Volca|tail -n 20 | awk '{print $2}' > id
for i in `cat id`; do kill -9 $i; done


## split by chromosome
for chr in {0..42}
do
samtools view -b $out/$ID.rm.bam $chr > $out/$ID.chr$chr.bam
done
index bam file
for chr in {0..42}
do
samtools index $out/$ID.chr$chr.bam
done

samtools view test.bam | awk '{print $5"\t"$9}'| sed '1i mapping-qulity\tmate-length' > test
awk '{print $5"\t"$9}' standard.bam| sed '1i mapping-qulity\tmate-length' > standard

cat 687_f1_test.fq | paste - - - - | sort -k1,1 -t " "  > 687_f1_test.sorte

grep -v -f 687_uniqid 687_f1_test.sorted | tr "\t" "\n" > 687_f1_test.sorted.fq

nohup zcat CRR061687_r2.filtered.fq.gz |paste - - - - |sort -k1,1 -S 500G | tr '\t' '\n' |gzip > CRR061687_r2.filtered_sorted.fq.gz &
  
nohup seqkit sort -n CRR061687_r2.filtered.fq.gz | gzip -c > CRR061687_r2.filtered_sorted.fq.gz &
  
bcftools merge chr0${chr}.vcf.gz chr0${chr}.vcf.gz chr0${chr}.vcf.gz chr${chr}.vcf.gz -o chr${chr}.all.vcf

bcftools filter 1000Genomes.vcf.gz --regions 9:4700000-4800000 > 4700000-4800000.vcf

#E6:xuebo@204:/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225
#gff & bed 取交集
bedtools intersect -a /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr${chr}.withBarley.vcf.gz -b merge_D.bed -header > chr${chr}.withBarley.vcf &

