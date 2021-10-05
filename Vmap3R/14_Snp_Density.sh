#统计VCF文件的snp的密度
#hg19.bed
chr1 248956422
chr2 242193529
chr3 198295559
....

bedtools makewindows -g hg19.bed -w 1000000 > windows.bed
bedtools coverage -a windows.bed -b test.vcf -counts > coverage.txt

for i in {1..42}
do
bedtools coverage -a ${i}.bed -b ../chr00${i}.vcf.gz -counts > chr00${i}.coverage.txt
done

#如果VCF文件过大，gz文件大于20G会报错，解决办法如下：

for i in {1,3,4,5,7,8,9}
#{1,3,7,9}
do
bedmap --echo --count --delim '\t' ${i}.bed <( gunzip -c ../chr00${i}.vcf.gz | vcf2bed --max-mem=100G --sort-tmpdir="~/large/dir" ) > chr00${i}.coverage.txt
done

for i in {10,11,13,14,15,16,17,19,20,21,22,23,25,26,27,28,29,31,32,33,34,35,37,38,39,40,41,42}
#{10,13,14,15,16,22,23,25,26,27,28,29,31,32,33,34,35,37,38,39,40,41,42}
do
bedmap --echo --count --delim '\t' ${i}.bed <( gunzip -c ../chr0${i}.vcf.gz | vcf2bed --max-mem=100G --sort-tmpdir="~/large/dir" ) > chr0${i}.coverage.txt
done


