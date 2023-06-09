#数据筛选及格式转换

#需要安装的软件：plink等
1.按MAF>0.05和缺失率<0.1过滤
plink --vcf test.imputed.vcf --maf 0.05 --geno 0.1 --recode vcf-iid --out test.filter --allow-extra-chr --autosome-num 42 
#（非数字染色体号ChrUn/Sy用此参数, 建议尽量把染色体号转成数字，另外需要对vcf中的标记ID进行编号）
2.对标记进行LD筛选
plink --vcf test.filter.vcf --indep-pairwise 50 10 0.2 --out test.filterLD --allow-extra-chr --autosome-num 42（.in文件里是入选的标记id）
3.提取筛选结果
plink --vcf test.filterLD.vcf --make-bed --extract test.filter.in --out  test.filter.prune.in
4.转换成structure/admixture格式
plink --bfile test.filter.prune.in --recode structure --out test.filter.prune.in  #生成. recode.strct_in为structure输入格式
plink --bfile test.filter.prune.in --recode 12 --out test.filter.prune.in  #生成.ped为admixture输入格式

#群体结构
#需要安装的软件：structure/admixture等
#这里我选了比较简单admixture来做 k值范围1到13
admixture --cv test.filter.ped 1 >>log.txt
admixture --cv test.filter.ped 2 >>log.txt
.......
admixture --cv test.filter.ped 13 >>log.txt
wait
grep "CV error" log.txt >k_1to13

#取CV error最小时的k值=10,  其中test.filter.prune.in.10.Q结果文件作为关联分析的输入源文件（去掉最后一列 添加表头和ID）

#亲缘关系/PCA
#需要安装的软件: tassel等(PCA分析可以用R包，已经做了群体结构这里就没做PCA分析)
run_pipeline.pl -importGuess test_impute.vcf -KinshipPlugin  -method Centered_IBS -endPlugin -export test_kinship.txt -exportType SqrMatrix

#关联分析
#需要安装的软件: tassel/GAPIT/FaSt-LMM等（ GAPIT很强大，要装很多R包，自动做图可视化。FaSt-LMM以后打算尝试一下）
#输入文件的格式需要手动修改一下，也比较简单 (如图)
1.vcf转hapmap格式
run_pipeline.pl -fork1 -vcf test.imputed.vcf  -export test -exportType Hapmap -runfork1
2.SNP位点排序
run_pipeline.pl -SortGenotypeFilePlugin -inputFile test.hmp.txt  -outputFile test_sort -fileType Hapmap
#得到test.hmp.txt
3. GLM模型
run_pipeline.pl -fork1 -h test_sort.hmp.txt -fork2 -r test.trait.txt -fork3 -q test.best_k10.txt -excludeLastTrait -combine4 -input1 -input2 -input3 -intersect -glm -export test_glm -runfork1 -runfork2 -runfork3
4.MLM模型
run_pipeline.pl -fork1 -h test_sort.hmp.txt -fork2 -r test.trait.txt -fork3 -q test.best_k10.txt -excludeLastTrait -fork4 -k test_kinship.txt -combine5  -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export test_mlm -runfork1 -runfork2 -runfork3
