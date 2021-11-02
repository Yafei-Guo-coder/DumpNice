#!/bin/bash
for i in {1..20}
do
mkdir /data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/bootstrap${i}
cd /data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/input/out_dir_msmc_${i}
nohup msmc2 -i 20 -t 20 -p 1*2+15*1+1*2 -o /data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/bootstrap${i}/within1_durum_bootstrap${i} -I 0,1,2,3 bootstrap_multihetsep.chr1.txt bootstrap_multihetsep.chr2.txt bootstrap_multihetsep.chr7.txt bootstrap_multihetsep.chr8.txt bootstrap_multihetsep.chr13.txt bootstrap_multihetsep.chr14.txt bootstrap_multihetsep.chr19.txt bootstrap_multihetsep.chr20.txt bootstrap_multihetsep.chr25.txt bootstrap_multihetsep.chr26.txt bootstrap_multihetsep.chr31.txt bootstrap_multihetsep.chr32.txt bootstrap_multihetsep.chr37.txt bootstrap_multihetsep.chr38.txt bootstrap_multihetsep.chr3.txt bootstrap_multihetsep.chr4.txt bootstrap_multihetsep.chr9.txt bootstrap_multihetsep.chr10.txt bootstrap_multihetsep.chr15.txt bootstrap_multihetsep.chr16.txt bootstrap_multihetsep.chr21.txt bootstrap_multihetsep.chr22.txt bootstrap_multihetsep.chr27.txt bootstrap_multihetsep.chr28.txt bootstrap_multihetsep.chr33.txt bootstrap_multihetsep.chr34.txt bootstrap_multihetsep.chr39.txt bootstrap_multihetsep.chr40.txt > nohupout1_${i} 2>& 1 &
  nohup msmc2 -i 20 -t 20 -p 1*2+15*1+1*2 -o /data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/bootstrap${i}/within2_landraceAB_bootstrap${i} -I 4,5,6,7 bootstrap_multihetsep.chr1.txt bootstrap_multihetsep.chr2.txt bootstrap_multihetsep.chr7.txt bootstrap_multihetsep.chr8.txt bootstrap_multihetsep.chr13.txt bootstrap_multihetsep.chr14.txt bootstrap_multihetsep.chr19.txt bootstrap_multihetsep.chr20.txt bootstrap_multihetsep.chr25.txt bootstrap_multihetsep.chr26.txt bootstrap_multihetsep.chr31.txt bootstrap_multihetsep.chr32.txt bootstrap_multihetsep.chr37.txt bootstrap_multihetsep.chr38.txt bootstrap_multihetsep.chr3.txt bootstrap_multihetsep.chr4.txt bootstrap_multihetsep.chr9.txt bootstrap_multihetsep.chr10.txt bootstrap_multihetsep.chr15.txt bootstrap_multihetsep.chr16.txt bootstrap_multihetsep.chr21.txt bootstrap_multihetsep.chr22.txt bootstrap_multihetsep.chr27.txt bootstrap_multihetsep.chr28.txt bootstrap_multihetsep.chr33.txt bootstrap_multihetsep.chr34.txt bootstrap_multihetsep.chr39.txt bootstrap_multihetsep.chr40.txt > nohupout2_${i} 2>& 1 &
  nohup msmc2 -i 20 -t 20 -p 1*2+15*1+1*2 -o /data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/bootstrap${i}/across_durum_landraceAB_bootstrap${i} -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 bootstrap_multihetsep.chr1.txt bootstrap_multihetsep.chr2.txt bootstrap_multihetsep.chr7.txt bootstrap_multihetsep.chr8.txt bootstrap_multihetsep.chr13.txt bootstrap_multihetsep.chr14.txt bootstrap_multihetsep.chr19.txt bootstrap_multihetsep.chr20.txt bootstrap_multihetsep.chr25.txt bootstrap_multihetsep.chr26.txt bootstrap_multihetsep.chr31.txt bootstrap_multihetsep.chr32.txt bootstrap_multihetsep.chr37.txt bootstrap_multihetsep.chr38.txt bootstrap_multihetsep.chr3.txt bootstrap_multihetsep.chr4.txt bootstrap_multihetsep.chr9.txt bootstrap_multihetsep.chr10.txt bootstrap_multihetsep.chr15.txt bootstrap_multihetsep.chr16.txt bootstrap_multihetsep.chr21.txt bootstrap_multihetsep.chr22.txt bootstrap_multihetsep.chr27.txt bootstrap_multihetsep.chr28.txt bootstrap_multihetsep.chr33.txt bootstrap_multihetsep.chr34.txt bootstrap_multihetsep.chr39.txt bootstrap_multihetsep.chr40.txt > nohupout3_${i} 2>& 1 &
  cd ../../
  monitor msmc2 21 30s 
done

[xuebo@lulab3 input]$ cat getchr1.sh
#!/bin/bash
/data1/home/xuebo/software/msmc-tools/generate_multihetsep.py \
--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr1.unmasked.bed.gz \
--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr1.unmasked.bed.gz \
--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr1.unmasked.bed.gz \
--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr1.unmasked.bed.gz \
--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr1.unmasked.bed.gz \
--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr1.unmasked.bed.gz \
--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr1.unmasked.bed.gz \
--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr1.unmasked.bed.gz \
--negative_mask=/data1/home/xuebo/Projects/Evo/smcpp/maskfile/mask1/chr1.masked.bed.gz \
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/B116.chr1.recode.vcf.gz \
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/B113.chr1.recode.vcf.gz \
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/B122.chr1.recode.vcf.gz \
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/XI_S4.chr1.recode.vcf.gz \
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/TW135.chr1.recode.vcf.gz \
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/TW160.chr1.recode.vcf.gz \
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/XI_33.chr1.recode.vcf.gz \
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/TW136.chr1.recode.vcf.gz \
> chr1.txt

(base) [xuebo@lulab3 input]$ head chr1.txt
1	1298138	77	AAGGAAAAGGAAAAAA
1	1298156	18	AAGGAAGGGGAAAAAA
1	1298251	95	AAGGAAAGGGAAAAAA
1	1298340	89	CCGGCCCGGGCCCCCC
1	1298348	8	TTCCTTTTCCTTTTTT
1	1298492	144	TTCCTTTCCCTTTTTT
1	1298493	1	TTCCTTTCCCTTTTTT
1	1298603	110	GGGGGGGAGGGGGGGG
1	1298660	57	AACCAAAACCAAAAAA
1	1298663	3	CCCCCCCTCCCCCCCC


(base) [xuebo@lulab3 input]$ cat getlist.sh
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
cat ../durum_landraceAB_8.txt | while read line
do
echo "/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/$line.chr$i.recode.vcf.gz" >> chr$i.list
done	
done

(base) [xuebo@lulab3 input]$ cat chr1.list
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/B116.chr1.recode.vcf.gz
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/B113.chr1.recode.vcf.gz
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/B122.chr1.recode.vcf.gz
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/XI_S4.chr1.recode.vcf.gz
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/TW135.chr1.recode.vcf.gz
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/TW160.chr1.recode.vcf.gz
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/XI_33.chr1.recode.vcf.gz
/data1/home/xuebo/Projects/Speciation/msmc/split_CCR/durum_landraceAB/vcf/TW136.chr1.recode.vcf.gz

(base) [xuebo@lulab3 input]$ cat gettxt.sh
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
num=$(wc -l chr$i.list | awk '{print $1}')
echo -e "#!/bin/bash" > getchr$i.sh
echo "/data1/home/xuebo/software/msmc-tools/generate_multihetsep.py \\" >> getchr$i.sh
for j in $(seq 1 $num)
do
echo "--mask=/data1/home/xuebo/Projects/Speciation/E3/chr/chr$i.unmasked.bed.gz \\" >> getchr$i.sh
done
echo "--negative_mask=/data1/home/xuebo/Projects/Evo/smcpp/maskfile/mask1/chr$i.masked.bed.gz \\" >> getchr$i.sh
cat chr$i.list | while read line
do
echo "$line \\" >> getchr$i.sh
done
echo " > chr$i.txt" >> getchr$i.sh 
echo "nohup sh getchr$i.sh &" >> run.sh
done

(base) [xuebo@lulab3 input]$ cat run.sh 
nohup sh getchr1.sh &
  nohup sh getchr2.sh &
  nohup sh getchr7.sh &
  nohup sh getchr8.sh &
  
  
  (base) [xuebo@lulab3 input]$ cd out_dir_msmc_1
(base) [xuebo@lulab3 out_dir_msmc_1]$ ls
bootstrap_multihetsep.chr1.txt   bootstrap_multihetsep.chr18.txt  bootstrap_multihetsep.chr26.txt  bootstrap_multihetsep.chr34.txt  bootstrap_multihetsep.chr42.txt
bootstrap_multihetsep.chr10.txt  bootstrap_multihetsep.chr19.txt  bootstrap_multihetsep.chr27.txt  bootstrap_multihetsep.chr35.txt  bootstrap_multihetsep.chr5.txt
bootstrap_multihetsep.chr11.txt  bootstrap_multihetsep.chr2.txt   bootstrap_multihetsep.chr28.txt  bootstrap_multihetsep.chr36.txt  bootstrap_multihetsep.chr6.txt
bootstrap_multihetsep.chr12.txt  bootstrap_multihetsep.chr20.txt  bootstrap_multihetsep.chr29.txt  bootstrap_multihetsep.chr37.txt  bootstrap_multihetsep.chr7.txt
bootstrap_multihetsep.chr13.txt  bootstrap_multihetsep.chr21.txt  bootstrap_multihetsep.chr3.txt   bootstrap_multihetsep.chr38.txt  bootstrap_multihetsep.chr8.txt
bootstrap_multihetsep.chr14.txt  bootstrap_multihetsep.chr22.txt  bootstrap_multihetsep.chr30.txt  bootstrap_multihetsep.chr39.txt  bootstrap_multihetsep.chr9.txt
bootstrap_multihetsep.chr15.txt  bootstrap_multihetsep.chr23.txt  bootstrap_multihetsep.chr31.txt  bootstrap_multihetsep.chr4.txt   nohupout1_1
bootstrap_multihetsep.chr16.txt  bootstrap_multihetsep.chr24.txt  bootstrap_multihetsep.chr32.txt  bootstrap_multihetsep.chr40.txt  nohupout2_1
bootstrap_multihetsep.chr17.txt  bootstrap_multihetsep.chr25.txt  bootstrap_multihetsep.chr33.txt  bootstrap_multihetsep.chr41.txt  nohupout3_1
(base) [xuebo@lulab3 out_dir_msmc_1]$ head bootstrap_multihetsep.chr1.txt
1	26076	312	GGAAGGGGGGGGAAGG
1	26634	558	CTCCCCCCCCCCCCCC
1	26642	8	CTCCCCCCCCCCCCCC
1	26661	19	GGGGAAGGGGGGGGGG
1	62304	198	TTCCCCTTTTTCCCTT
1	63994	321	CCCCCCCTCCCCCCCC
1	64007	13	GGAAAAGAGGGAAAGG
1	65093	61	TTCCCCTTTTTCCCTT
1	65102	9	GGAAGGGGGGGAAAGG
1	65109	7	TTCCCCTTTTTTCCTT


(base) [xuebo@lulab3 vcf]$ cat cp1.sh 
#!/bin/bash
cat ../durum_4.txt | while read line
do
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
cp /data1/home/xuebo/Projects/Speciation/msmc/Durum/vcf/$line.chr$i.recode.vcf.gz ./
  cp /data1/home/xuebo/Projects/Speciation/msmc/Durum/vcf/$line.chr$i.recode.vcf.gz.tbi ./
  done
done