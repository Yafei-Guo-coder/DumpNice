#############################################volcano Finder#######################################
# VolcanoFinder is a program to perform genome-wide scans of adaptive introgression by using the
# spatial distribution of within-species polymorphism and between species substitution in the genome.
# input: 1.Allele frequency file 2.site frequency spectrum(sfs)3.ancestral file
##工作路径: /data1/home/xuebo/Projects/Speciation/volcanofinder

###### installation ######
wget -c http://degiorgiogroup.fau.edu/vf.html
tar -xzvf volcanofinder_v1.0.tar.gz
cd volcanofinder_v1.0
make
vim ~/.bashrc
source ~/.bashrc

###### prepare allele frequency file ######
# 1.calculate allele numbers of ref and al
# ex:chr1	821443	821444	T	1137	C	3
vcftools --gzvcf /data1/home/xuebo/Projects/Speciation/E3/chr36.all.vcf.gz  --counts --out chr36_part1
cat chr36_part1.frq.count |tail -n +2 |awk '{print $1"\t"$2-1"\t"$2"\t"substr($5,1,1)"\t"substr($5,3,length($5))"\t"substr($6,1,1)"\t"substr($6,3,length($6))}'>chr36_part1.temp
# 2.extract and non-polarized sites to obtain Freq file
# ex:position x n folded
# 460000 9 100 0
# 460010 100 100 0
# 460210 30 78 1
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' A_pos_ances.txt >  A_pos_ances.bed 
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' B_pos_ances.txt >  B_pos_ances.bed 
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' D_pos_ances.txt >  D_pos_ances.bed 
mkdir -p Dsplit_ances
for chr in `cut -f 1 D_pos_ances.bed | sort | uniq`; do
                grep -w $chr D_pos_ances.bed > Dsplit_ances/$chr.pos_ances.bed
done
***intersect,取交集，并且找出谁是ancestral，格式是pos\tderived\tderived+ancestral
bedtools intersect -a ../anc/Dsplit_ances/36.pos_ances.bed  -b chr36_part1.temp -wo |awk '{if($4==$8)print $2"\t"$11"\t"$11+$9"\t""0";else print $2"\t"$9"\t"$9+$11"\t""0" }' >polarized.freq
***找出-a里面的特有的pos
bedtools subtract -b ../anc/Dsplit_ances/36.pos_ances.bed  -a chr36_part1.temp |awk '{print $2"\t"$7"\t"$5+$7"\t""1"}' >>polarized.freq
# remove undistinguished sites
awk '{if($2==0)print $0 }' polarized.freq >x
awk '{if($2==$3&&$4==1)print $0 }' polarized.freq >>x
sort -k1,1n x>xx
awk 'FNR==NR{a[$0]++;next}(!($0 in a))' xx polarized.freq >FreqFile
sort -k1,1n FreqFile |sed '1i position\tx\tn\tfolded' > FreqFilesort
# 3.sum to sfs Spect file
# 604 is alleles number
for var in {1..604};do awk -v j=$var '{if($2==j)print $2}' FreqFile| wc -l >> tempp ; done


# 104072 is snp number  zcat /data1/home/xuebo/Projects/Speciation/E3/chr36.all.vcf.gz | grep -v "#" | wc -l
awk '{print NR"\t"$1/104072}' tempp >SpectFile

###### run VolcanoFinder by blocks #######
# g:block length(ex:1000bp) G:step length(ex:1000snp)
#ex:VolcanoFinder -big 1000 FreqFile SpectFile -1 1 1 xbsj 1 50
#VolcanoFinder -big g FreqFile SpectFile D P MODEL OutFile BLOCK NBLOCK
#VolcanoFinder -bi G FreqFile SpectFile D P MODEL OutFile BLOCK NBLOCK
VolcanoFinder -big 100 FreqFilesort SpectFile -1 1 1 xbsj 1 50

# merge results
#VolcanoFinder -m xbsj 50
VolcanoFinder -m OutFile NBLOCK