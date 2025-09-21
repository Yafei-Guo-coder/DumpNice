#愈伤芽分化项目
3	3527280	3528440	.	+	.	ID=gene:AT3G11260;Name=WOX5
4	11555039	11560771	.	+	.	ID=gene:AT4G21750;Name=ATML1
1	28041146	28043022	.	+	.	ID=gene:AT1G74650;Name=MYB31
5	26634275	26635885	.	+	.	ID=gene:AT5G66700;Name=HB53

#snpeff
grep AT3G11260 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.genes.txt | datamash transpose
grep AT4G21750 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.genes.txt | datamash transpose
grep AT1G74650 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.genes.txt | datamash transpose
grep AT5G66700 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.genes.txt | datamash transpose

#accession信息
sed 's/[^[:print:]\t]//g' accession.txt > output.txt

#vcf文件下载自网站 https://tools.1001genomes.org/vcfsubset/#select_strains
awk -F"\t" '{if($5 != "." ) print $0}' AT3G11260.vcf > AT3G11260.filter.vcf 
awk -F"\t" '{if($5 != "." ) print $0}' AT4G21750.vcf > AT4G21750.filter.vcf 
awk -F"\t" '{if($5 != "." ) print $0}' AT1G74650.vcf > AT1G74650.filter.vcf 
awk -F"\t" '{if($5 != "." ) print $0}' AT5G66700.vcf > AT5G66700.filter.vcf 

验证思路：
	1. 先画一个单倍型分布吧，看能不能找到相对比较关键的在自然群体中的位点，然后在画个地图上的地理分布，看一下跟经纬度有没有什么关联。预期结果，如果跟群体结构是一致的，那只能说明是正常分化，无法说明是环境适应。
	2. 是否需要做一个跟环境（光照，温度，水分）的关联分析，看这些基因是否参与响应了环境变化。
	3. 基因的单倍型演化过程？跟什么过程相关？比如说地理分布？拟南芥的传播？响应环境变化？
比较基因组上可以做一些什么分析？参考康的光合作用那块分析？

#画单倍型图谱
#准备输入文件
for i in AT1G74650 AT3G11260 AT4G21750 AT5G66700
do
#plink --vcf ${i}.filter.vcf --recode vcf-iid --out ${i}-geno
grep -v "##" <(zcat ${i}-geno.vcf.gz) | awk '{$1=null;$3=null;$6=null;$7=null;$8=null;$9=null;print $0}'  | sed 's/ //' |sed 's/  / /1' | sed 's/    //' | sed 's/0\/0/0/g' | sed 's/0\/1/1/g' | sed 's/1\/1/2/g' | sed 's/\.\/\./-1/g' | sed 's/2\/2/6/g' | sed 's/0\/2/3/g' | sed 's/1\/2/4/g' > ${i}-vcf-geno.txt
done





