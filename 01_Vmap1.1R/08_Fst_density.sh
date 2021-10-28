#工作路径：xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%_Fst
#EU_South*
#North2_South*
#Strang_WA*
#Tibet_South*
#WA_EU*
#WA_South*
#neg_North1_North2*
#neg_WA_North1*
#neg_WA_North2*

#GenWin(smooth)的结果是这样的"WindowStart" "WindowStop" "SNPcount" "MeanY" "Wstat"，没有染色体，要加上染色体的信息
#smooth结果添加染色体号并进行排序，并且合并成A，B，D lineage.

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

#克隆受选择基因fst的分布
