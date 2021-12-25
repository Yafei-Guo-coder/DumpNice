#按MAF>0.05和缺失率<0.1过滤
plink --vcf test.imputed.vcf --maf 0.05 --geno 0.1 --recode vcf-iid --out test.filter --allow-extra-chr --double-id --autosome-num 42
--keep-only-indels

#对VCF进行LD筛选
plink --vcf test.filter.vcf --indep-pairwise 50 10 0.2 --out test.filterLD --allow-extra-chr --double-id --autosome-num 42

#提取筛选结果
#plink --vcf test.filter.vcf --make-bed --extract test.filterLD.prune.in --out  test.filter.prune.in --double-id --autosome-num 42
plink --vcf test.filter.vcf --extract test.filterLD.prune.in --recode vcf-iid --out  test.filter.prune.in --double-id --autosome-num 42

#转换成structure/admixture格式
plink --bfile test.filter.prune.in --recode structure --out test.filter.prune.in  #生成. recode.strct_in为structure输入格式
plink --bfile test.filter.prune.in --recode 12 --out test.filter.prune.in --autosome-num 42 #生成.ped为admixture输入格式
plink --bfile test.filter.prune.in --recode vcf --out test.filter.prune.in --autosome-num 42

#admixture
admixture --cv test.filter.ped 1 >>log.txt
admixture --cv test.filter.ped 2 >>log.txt
.......
admixture --cv test.filter.ped 13 >>log.txt
wait
grep "CV error" log.txt > k_1to13
