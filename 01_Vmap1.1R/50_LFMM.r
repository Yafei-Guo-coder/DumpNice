#LFMM
#除了用bayenv检测与环境相互关联的信号，也使用LFMM进行检测
#Working directory: yafei@204:/data1/home/yafei/003_Project3/Structure/LFMM
#原始文件是149个在xp-clr中使用的群体, (剔除了TW095,因为没有环境变量数据）的SNP文件(A_noMiss.vcf.gz  B_noMiss.vcf.gz  D_noMiss.vcf.gz)
A:2593384 B:1629490 D:1483980
#准备输入文件
#genotype.lfmm(行是样本，列是位点，无行名列名，以空格分开)
vcftools --gzvcf A_noMiss.vcf.gz --012 --out A.matrix
awk '{$1=null;print $0}' A.matrix.012 | sed 's/\t/ /g'| sed 's/^ //' > A.genotype.lfmm
vcftools --gzvcf B_noMiss.vcf.gz --012 --out B.matrix
awk '{$1=null;print $0}' B.matrix.012 | sed 's/\t/ /g'| sed 's/^ //' > B.genotype.lfmm
vcftools --gzvcf D_noMiss.vcf.gz --012 --out D.matrix
awk '{$1=null;print $0}' D.matrix.012 | sed 's/\t/ /g'| sed 's/^ //' > D.genotype.lfmm

#gradient.env(行是样本，列是环境变量的值)
awk 'NR==FNR{b[$1]=$1; a[$1]=$0}NR!=FNR{if($1 in b) print a[$1]}' 225env.txt 149.taxa.txt | awk '{print $2}' > elevation.txt

#估计群体结构
先对总体进行LD筛选
plink --vcf D_noMiss.vcf.gz --indep-pairwise 50 10 0.2 --out D.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf D_noMiss.vcf.gz --extract D.filterLD.prune.in --recode vcf-iid --out D.prune.in --double-id --autosome-num 42
vcftools --vcf A.prune.in.vcf --012 --out A.LD.matrix
awk '{$1=null;print $0}' A.LD.matrix.012 | sed 's/\t/ /g'| sed 's/^ //' > A.LD.genotype.lfmm
vcftools --vcf B.prune.in.vcf --012 --out B.LD.matrix
awk '{$1=null;print $0}' B.LD.matrix.012 | sed 's/\t/ /g'| sed 's/^ //' > B.LD.genotype.lfmm
vcftools --vcf D.prune.in.vcf --012 --out D.LD.matrix
awk '{$1=null;print $0}' D.LD.matrix.012 | sed 's/\t/ /g'| sed 's/^ //' > D.LD.genotype.lfmm

library(LEA)
genotype = lfmm2geno("A.LD.genotype.lfmm")
obj.snmf = snmf(genotype, K = 1:15, entropy = T, ploidy = 2, project="new")
pdf("A.structure.pdf")
plot(obj.snmf)
barplot(t(Q(obj.snmf, K = 3)), col = 1:3)
barplot(t(Q(obj.snmf, K = 4)), col = 1:4)
barplot(t(Q(obj.snmf, K = 5)), col = 1:5)
barplot(t(Q(obj.snmf, K = 6)), col = 1:6)
barplot(t(Q(obj.snmf, K = 7)), col = 1:7)
barplot(t(Q(obj.snmf, K = 8)), col = 1:8)
barplot(t(Q(obj.snmf, K = 9)), col = 1:9)
barplot(t(Q(obj.snmf, K = 10)), col = 1:10)
barplot(t(Q(obj.snmf, K = 11)), col = 1:11)
barplot(t(Q(obj.snmf, K = 12)), col = 1:12)
barplot(t(Q(obj.snmf, K = 13)), col = 1:13)
dev.off()
genotype = lfmm2geno("D.LD.genotype.lfmm")
obj.snmf = snmf(genotype, K = 1:15, entropy = T, ploidy = 2, project="new")
pdf("D.structure.pdf")
plot(obj.snmf)
barplot(t(Q(obj.snmf, K = 3)), col = 1:3)
barplot(t(Q(obj.snmf, K = 4)), col = 1:4)
barplot(t(Q(obj.snmf, K = 5)), col = 1:5)
barplot(t(Q(obj.snmf, K = 6)), col = 1:6)
barplot(t(Q(obj.snmf, K = 7)), col = 1:7)
barplot(t(Q(obj.snmf, K = 8)), col = 1:8)
barplot(t(Q(obj.snmf, K = 9)), col = 1:9)
barplot(t(Q(obj.snmf, K = 10)), col = 1:10)
barplot(t(Q(obj.snmf, K = 11)), col = 1:11)
barplot(t(Q(obj.snmf, K = 12)), col = 1:12)
barplot(t(Q(obj.snmf, K = 13)), col = 1:13)
dev.off()
genotype = lfmm2geno("B.LD.genotype.lfmm")
obj.snmf = snmf(genotype, K = 1:15, entropy = T, ploidy = 2, project="new")
pdf("B.structure.pdf")
plot(obj.snmf)
barplot(t(Q(obj.snmf, K = 3)), col = 1:3)
barplot(t(Q(obj.snmf, K = 4)), col = 1:4)
barplot(t(Q(obj.snmf, K = 5)), col = 1:5)
barplot(t(Q(obj.snmf, K = 6)), col = 1:6)
barplot(t(Q(obj.snmf, K = 7)), col = 1:7)
barplot(t(Q(obj.snmf, K = 8)), col = 1:8)
barplot(t(Q(obj.snmf, K = 9)), col = 1:9)
barplot(t(Q(obj.snmf, K = 10)), col = 1:10)
barplot(t(Q(obj.snmf, K = 11)), col = 1:11)
barplot(t(Q(obj.snmf, K = 12)), col = 1:12)
barplot(t(Q(obj.snmf, K = 13)), col = 1:13)
dev.off()

#204跑D，k=6(conda R4)。203跑B，k=5(conda r4)。66跑A，k=7(conda r4)。
for i in `cat names`; do sed "s@env/@$i@" test.r |sed "1i setwd(\"/data1/home/yafei/003_Project3/LFMM/A_out/$i\")"> $i/$i.r; done

library(LEA)
obj.lfmm = lfmm("D.genotype.lfmm", "env/bio1.env", K = 6, rep = 5, project="new")
#Record z-scores from the 5 runs in the zs matrix
zs = z.scores(obj.lfmm, K = 6)
#Combine z-scores using the median
zs.median = apply(zs, MARGIN = 1, median)
lambda = median(zs.median^2)/qchisq(0.5, df = 1)
pdf("adj.p.values.pdf")
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values, col = "red")
# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(zs.median^2/.55, df = 1, lower = FALSE)
#histogram of p-values
hist(adj.p.values, col = "green")
dev.off()
## FDR control: Benjamini-Hochberg at level q
## L = number of loci
L = 500
#fdr level q
q = 0.01
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]
write.table(candidates.bh,"site1.txt",row.names=F, quote = F)
## FDR control: Storey's q-values
library(qvalue)
plot(qvalue(adj.p.values))
candidates.qv = which(qvalue(adj.p.values, fdr = .1)$signif)
write.table(candidates.qv,"site2.txt",row.names=F, quote = F)

for i in `cat names`; do sed "s/qq.pdf/$i.qq.pdf/" qq.sh | sed "s/qvalue.pdf/$i.qvalue.pdf/" | sed "s/pos.txt/$i.pos.txt/"> $i/qq.r; done
for i in `cat names` ; do echo "cd /data2/yafei/003_Project3/LFMM/D_Out/"$i;echo "nohup Rscript qq.r &"; done > run.qq.sh
nohup bash run.qq.sh > run.qq.log &
  
for i in `ls *smooth_D.top5.bed`
do
bedtools intersect -a ../D.select.pos.bed -b  $i -wa |sort | uniq > bayenv_xpclr/$i
done
