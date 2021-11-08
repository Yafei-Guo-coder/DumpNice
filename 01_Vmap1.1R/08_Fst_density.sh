#工作路径：yafei@203:/data2/yafei/003_project3/Project3/FST_group/VmapData/FST_selection
#计算受选择位点Fst的密度分布
for i in `ls *txt`; do  wc -l $i; done | awk '{print $2"\t"$1"\t"$1*0.01"\t"$1*0.05}' | awk -F"[.|\t]" '{print $1"\t"$3"\t"$4"\t"$6}' | awk '{print "tail -n "$4,$1".txt > Top5%_Fst/"$1".top5.bed"}'
for i in `ls *bed`
do
bedtools intersect -b ../gene_v1.1_Lulab.gff3 -a $i -wo |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$9"\t"$13}' | awk -F";" '{print $1}' | awk -F"ID=" '{print $1""$2}'| sed '1i Chr\tRegion1\tRegion2\tXp-clr\tStart\tStop\tID' >${i::-3}gene.txt
done
for i in `ls *gene.txt`
do
awk 'NR==FNR{a[$2]=$1;b[$2]=$2}NR!=FNR{if($7 in b) print $0"\t"a[$7]}' ../all_cloned_gene.txt $i |sed '1i Chr\tRegion1\tRegion2\tXp-clr\tStart\tStop\tID\tName' > ${i::-8}xp-clr.txt
done
for i in `ls *txt`
do
sed '1d' $i > Top5%_Fst/${i::-3}all.xp.txt
done
for i in `ls *all.xp.txt`
do
bedtools intersect -b ../clone.gene.txt -a $i -wo | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$8"\t"$9"\t"$10}' |sed '1i Chr\tRegion1\tRegion2\tXp-clr\tStart\tStop\tID' >${i::-10}all.xp-clr.txt
done
#克隆受选择基因fst的分布
for i in `ls *txt`; do awk '{print $1"\t"$5"\t"$6"\t"7"\t"$8"\t"$2"\t"$3"\t"$4}' $i | sed '1d' >${i::-15}format1.txt ; done
bash mv.sh
for i in `cat name_prefix.txt`
do
bedtools intersect -b ${i}.all.fst -a ${i}.format1.txt  -wo |sed '1i Chr\tStart\tStop\tID\tName\tRegion1\tRegion2\tXp-clr\tCHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$14}' > ${i}.format2.txt
done
#生物胁迫，非生物胁迫，背景基因fst的分布(所有基因)
#input file: ../abioticgene.txt ../backgroudgene.txt ../bioticgene.txt
grep -w -f ../bioticgene.txt ../gene_v1.1_Lulab.gff3 |awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | awk -F"ID=" '{print $1$2}' > bioticgene.bed
for i in `cat name_prefix.txt`
do
bedtools intersect -b ${i}.all.fst -a bioticgene.bed -wo | awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$2"\t"$3"\t"$4}' | sed '1i CHROM\tBIN_START\tBIN_END\tMEAN_FST\tGene_start\tGene_end\tGene_id'> ${i}.bioticgene.txt
don
grep -w -f ../abioticgene.txt ../gene_v1.1_Lulab.gff3 |awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | awk -F"ID=" '{print $1$2}' > abioticgene.bed
for i in `cat name_prefix.txt`
do
bedtools intersect -b ${i}.all.fst -a abioticgene.bed -wo | awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$2"\t"$3"\t"$4}' | sed '1i CHROM\tBIN_START\tBIN_END\tMEAN_FST\tGene_start\tGene_end\tGene_id'> ${i}.abioticgene.txt
done
grep -w -f ../backgroudgene.txt ../gene_v1.1_Lulab.gff3 |awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | awk -F"ID=" '{print $1$2}' > backgroudgene.bed
for i in `cat name_prefix.txt`
do
bedtools intersect -b ${i}.all.fst -a backgroudgene.bed -wo | awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$2"\t"$3"\t"$4}' | sed '1i CHROM\tBIN_START\tBIN_END\tMEAN_FST\tGene_start\tGene_end\tGene_id'> ${i}.backgroudgene.txt
done
#发育，开花，抗病基因fst的分布(已克隆基因)
#input file: ../developgene.txt ../flowergene.txt ../resistgene.txt
grep -w -f ../developgene.txt ../gene_v1.1_Lulab.gff3 |awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | awk -F"ID=" '{print $1$2}' > developgene.bed
for i in `cat name_prefix.txt`
do
bedtools intersect -b ${i}.all.fst -a developgene.bed -wo | awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$2"\t"$3"\t"$4}' | sed '1i CHROM\tBIN_START\tBIN_END\tMEAN_FST\tGene_start\tGene_end\tGene_id'> ${i}.developgene.txt
done
grep -w -f ../flowergene.txt ../gene_v1.1_Lulab.gff3 |awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | awk -F"ID=" '{print $1$2}' > flowergene.bed
for i in `cat name_prefix.txt`
do
bedtools intersect -b ${i}.all.fst -a flowergene.bed -wo | awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$2"\t"$3"\t"$4}' | sed '1i CHROM\tBIN_START\tBIN_END\tMEAN_FST\tGene_start\tGene_end\tGene_id'> ${i}.flowergene.txt
done
grep -w -f ../resistgene.txt ../gene_v1.1_Lulab.gff3 |awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | awk -F"ID=" '{print $1$2}' > resistgene.bed
for i in `cat name_prefix.txt`
do
bedtools intersect -b ${i}.all.fst -a resistgene.bed -wo | awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$2"\t"$3"\t"$4}' | sed '1i CHROM\tBIN_START\tBIN_END\tMEAN_FST\tGene_start\tGene_end\tGene_id'> ${i}.resistgene.txt
done

for i in `cat name_prefix.txt`
do
bedtools intersect -b ${i}.all.fst -a 266gene.bed -wo | awk '{print $5"\t"$6"\t"$7"\t"$10"\t"$2"\t"$3"\t"$4}' | sed '1i CHROM\tBIN_START\tBIN_END\tMEAN_FST\tGene_start\tGene_end\tGene_id'> fst_266clone/${i}.266gene.txt
done

for i in `cat ../name_prefix.txt`
do
awk 'NR==FNR{a[$1]=$1;b[$1]=$2;c[$1]=$3}NR!=FNR{if($7 in a) print $0"\t"b[$7]"\t"c[$7]}' AllGeneAnnoMap.txt $i.266gene.txt |sort -k4,4n | sed '1i CHROM\tBIN_START\tBIN_END\tMEAN_FST\tGene_start\tGene_end\tGene_id\tAnno\tName'> $i.266gene.out.txt
#awk 'NR==FNR{a[$2]=$2;b[$2]=$1}NR!=FNR{if($7 in a) print $0"\t"b[$7]}' ../../all_cloned_gene.txt $i.txt | 
#rm $i.txt
done

#按照MEAN_FST进行排序。
#转移到本地: /Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density
#通过08_Fst_density.r进行统计画图
