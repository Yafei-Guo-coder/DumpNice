#awk
awk 'NR==FNR{a[$2]=$0;}NR!=FNR{print $0,a[$2]}' b.txt a.txt
#内关联
#方法1
awk -F',' 'NR==FNR{a[$1]=$2;}NR!=FNR && a[$1] {print $0,a[$1]}' b.txt a.txt
#方法2
awk -F',' 'NR==FNR{a[$1]=$2;}NR!=FNR && $1 in a {print $0,a[$1]}' b.txt a.txt
awk 'FNR==NR{a[$1];next}($1 in a){next} {print}' a b 
awk 'ORS=NR%2?" ":"\n"{print}' 
  
#shell
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
zcat CRR061687_r2.filtered.fq.gz |paste - - - - |sort -k1,1 -S 500G | tr '\t' '\n' |gzip > CRR061687_r2.filtered_sorted.fq.gz
seqkit sort -n CRR061687_r2.filtered.fq.gz | gzip -c > CRR061687_r2.filtered_sorted.fq.gz
  
#bedtools
bedtools intersect -a /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/chr${chr}.withBarley.vcf.gz -b merge_D.bed -header > chr${chr}.withBarley.vcf
bedtools intersect -a /data1/home/yafei/009_GWAS/gene/gene_v1.1_Lulab.gff3 -b test -wb | awk 'split($9, array, ";") {print $1"\t"$4"\t"$5"\t"array[1]"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' 
#bcftools
bcftools merge chr0${chr}.vcf.gz chr0${chr}.vcf.gz chr0${chr}.vcf.gz chr${chr}.vcf.gz -o chr${chr}.all.vcf
bcftools filter 1000Genomes.vcf.gz --regions 9:4700000-4800000 > 4700000-4800000.vcf

#vcftools
vcftools --gzvcf file.vcf.gz --positions specific_position.txt --recode --recode-INFO-all --out specific_position.vcf


#getVcf.sh
cat $1 |while read gene chr from to taxa
do
    #echo $chr $from $to
    #if echo $2 |grep -q '.*.vcf.gz$';then
    #    vcftools --gzvcf $2 --chr $chr --from-bp $from --to-bp $to  --recode --recode-INFO-all --out $gene.$chr.$from-$to 
    #elif echo $2 |grep -q '.*.vcf$';then
        vcftools --gzvcf /data2/yafei/003_Project3/Vmap1.1/E6/VCF/chr${chr}.E6all.vcf.gz --keep ${taxa}_ancestor.txt --maf 0.0001 --chr $chr --from-bp $from --to-bp $to  --recode --recode-INFO-all --out 5k.$gene.$chr.$from-$to
    #fi
done

vcf-concat AB_noMiss_0.05.vcf.gz D_noMiss_0.05.vcf.gz | bgzip -c > noSort_noMiss_0.05.vcf.gz 


#bedtools
#提取vcf的特定区域
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
  
#删除行首空格
sed 's/^[ \t]*//g'
#删除行尾空格
sed 's/[ \t]*$//g'
#删除空行
sed '/^$/d' 
#替换多个空格为一个逗号
sed 's/\s\+/,/g'
#计算bam的depth
for i in `cat depth.txt`
do
 tail -n1 *${i} | awk '{print $4}' | sed '/^$/d' |awk 'BEGIN{count=0} {count = count+$1} END{print count/NR}'
done

for i in `cat fqlist.txt`
do
bwa mem -t 20 -R '@RG\tID:Aegilops\tPL:illumina\tSM:Aegilops' /data1/publicData/wheat/reference/v1.0/D/bwaLib/d_iwgscV1.fa.gz ${i}_1.fq.gz ${i}_2.fq.gz | samtools view -S -b -> ${i}.bam
samtools sort -n -m 4G -@ 20 -o ${i}.namesort.bam -O bam ${i}.bam && samtools fixmate -@ 20 -m ${i}.namesort.bam ${i}.fixmate.bam && samtools sort -m 4G -@ 20 -o ${i}.fixmate.pos.bam -O bam ${i}.fixmate.bam && rm -f ${i}.namesort.bam && samtools markdup -@ 20 -r ${i}.fixmate.pos.bam ${i}.rmdup.bam && rm -f ${i}.fixmate.bam && rm -f ${i}.fixmate.pos.bam
done

Triticum

#国内github镜像s
https://hub.fastgit.org/

datamash

for i in A B D
do
for j in `cat names`
do
grep -v -f neg_WA_North2_smooth_${i}.top5.txt ${j}_smooth_${i}.top5.txt > ${j}_${i}.go.gene.txt
done
done

tar -zxvf wheatGO-v1.1.tar.gz 

conda create -n R4
conda activate R4
conda install -c bioconda bioconductor-clusterprofiler
conda update R

#计算bam文件的depth
#bams.txt格式为 /data3/wgs/bam/ABD/ABD_0165.bam，一行一个bam文件
for i in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
for j in `cat JIC_bam_D.txt`
do
var=${j##*/} 
mosdepth -c ${i} -n -t 32 JIC_D/out_${i}_${var::-4} $j
done
done
java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage/Alineage_withBarley.vcf.gz --out Alineage_withBarley.all.ibs.txt > logA.txt

#合并test_L1.bam和test_L2.bam文件
samtools merge -h test.sam test_L1_L2.bam test_L1.sorted.bam test_L2.sroted.bam
cat $1 |while read num file
do
	bedtools sample -n ${num} -i ${file} -header | bgzip -c > ${file::-13}.shuf.vcf.gz &
done

samtools mpileup -A -B -q 30 -Q 20 -f /data1/home/xinyue/ref/byChr/chr032.fa.gz /data3/wgs/bam/ABD/ABD_0372.bam -r 32:89999900-89999996

vcf-concat *shuf.vcf.gz |awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1n -k2,2n"}' > Alineage.10000.vcf &
vcftools --vcf file1.snp.vcf --diff file2.snp.vcf --diff-site --out Diff.site


cd Out/vcf/VCF/
sed '/^#/d' chr001.vcf | awk '{print $1"\t"$2"\t"$10"\t"$11"\t"$12}'  > ../../../pos.txt
cd ../../../
zcat trueSet.txt.gz | datamash transpose| awk '{$2=$2$2;if($2==$4 && $4==$5 && $5==$6) next;else print $0}' > test.out
awk 'NR==FNR{a[$2]=$2;b[$2]==$0}NR!=FNR{if($1+1 in a) print $0}' pos.txt test.out
awk 'NR==FNR{a[$2]=$2;b[$2]==$0}NR!=FNR{if($1+1 in a) print "ok";else print $0}' pos.txt test.out
samtools mpileup -A -B -q 30 -Q 20 -f test_1M.fa.gz sample2.rmdup.bam -r 1:489865-489890

make install DESTDIR=/data1/home/yafei/008_Software/gsl-2.7
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data1/home/yafei/008_Software/gsl-2.7/usr/local

#PCA
plink --vcf ABlineage.maf0.05.5k.vcf --pca header tabs -out ABlineage.maf0.05.5k --double-id --autosome-num 42
#MDS
plink --vcf ABlineage.maf0.05.60k.vcf.gz --mds-plot 10 eigendecomp --cluster --double-id --autosome-num 42 --out ABlineage.maf0.05.60k


