reference
for i in {1..42}
do
WGS --model file --type ByChr --chr $i --file unmasked.num.txt --out unmasked.chr$i &
  done


unmask
#!/bin/bash
for chr in {1..42}
do
WGS --model vcf --type siteOverlap --file /data1/home/xuebo/Projects/Speciation/E3/chr${chr}.all.vcf.gz --file2 /data1/home/xuebo/Projects/Speciation/smcpp/reference/unmasked.chr${chr}.txt --out /data1/home/xuebo/Projects/Speciation/smcpp/E3unmasked/chr${chr}.vcf &
  done

exonremove
#!/bin/bash
for chr in {1..42}
do
WGS --model vcf --type unintersect --file /data1/home/xuebo/Projects/Speciation/smcpp/E3unmasked/chr${chr}.vcf.gz --file2 /data1/home/xuebo/Projects/Speciation/gff/chr${chr}.exon --out /data1/home/xuebo/Projects/Speciation/smcpp/E3removeExon/chr${chr}.vcf &
  done

for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
zcat chr${i}.vcf.gz | wc -l >> getAlineagesnp.txt
done

#!/bin/bash
for i in {1..42}
do
java -Xss64g -Xmx500g -jar /data1/home/xuebo/software/beagle.28Sep18.793.jar nthreads=20 gt=/data1/home/xuebo/Projects/Speciation/smcpp/E3removeExon/"chr"$i".vcf.gz" out="chr"$i".beagle" phase-segment=1 burnin=3 iterations=6 phase-states=100 window=1 overlap=0.001  &
  monitor java 10  60s
done




#!/bin/bash
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
sed -i '15 r /data1/home/xuebo/Projects/Speciation/smcpp/sample15/group/vcf_contifheader.txt' chr$i.vcf
bgzip -c chr$i.vcf >  chr$i.vcf.gz
tabix chr$i.vcf.gz
done


#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/sample15/domemmer/V1
cd /data1/home/xuebo/Projects/Speciation/smcpp/sample15/domemmer/V1
sed 's/\t/_/g' /data1/home/xuebo/Projects/Speciation/smcpp/sample15/group/dom_emmer_15.txt |  tr "\n" "," | sed s'/.$//' > group_domemmer_header.txt &
  mkdir data
name=$(sed -n 1p group_domemmer_header.txt)
echo $name
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
smc++ vcf2smc --mask /data1/home/xuebo/Projects/Evo/smcpp/maskfile/mask1/chr$i.masked.bed.gz /data1/home/xuebo/Projects/Speciation/smcpp/sample15/domemmer/vcf/chr$i.vcf.gz \
/data1/home/xuebo/Projects/Speciation/smcpp/sample15/domemmer/V1/data/chr$i.smc.gz $i domemmer:$name  &
  done

get.mask
#!/bin/bash
for i in {1..10}
do
cd ${i}x
for chr in {1..2}
do
WGS --model vcf --type unintersect --file chr${chr}.unmask.vcf --file2 /data1/home/xuebo/Projects/Speciation/gff/chr${chr}.exon --out chr${chr}.mask.vcf &
  done
cd ..
done

get.unmask
#!/bin/bash
for i in {1..10}
do
mkdir ${i}x
cd ${i}x
for chr in {1..2}
do
WGS --model vcf --type siteOverlap --file /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/depth/vcf/final_vcf/${i}x/chr${chr}.final_snp.vcf.gz \
--file2 /data1/home/xuebo/Projects/Speciation/smcpp/reference/unmasked.chr${chr}.txt --out chr${chr}.unmask.vcf &
  done
cd ..
done



get.diplod
#!/bin/bash
for i in {1..10}
do
mkdir ${i}x
cd ${i}x
mkdir vcf
for chr in {1..2}
do
WGS --model vcf --type GenerateDiploid --file /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/depth/smcpp/maskfile/${i}x/chr${chr}.mask.vcf.gz --file2 ../Landrace_15.txt --out ./vcf/chr${chr}.vcf &
  done
cd ..
done

cat getsmcformat.sh 
#!/bin/bash
for i in {1..10}
do
cd ${i}x
mkdir data
name=$(sed -n 1p ../group_LandraceAB_header.txt)
echo $name
for chr in {1..2}
do
smc++ vcf2smc --mask /data1/home/xuebo/Projects/Evo/smcpp/maskfile/mask1/chr${chr}.masked.bed.gz ./vcf/chr${chr}.vcf.gz ./data/chr${chr}.smc.gz ${chr} LandraceAB:$name
done
cd ..
done
cat estimate.sh 
#!/bin/bash
for i in {1..10}
do
cd ${i}x
mkdir analysis
OMP_NUM_THREADS=1 smc++ estimate 6.5e-9 ./data/chr*.smc.gz --polarization-error 0.5 -o ./analysis --cores 30 --spline cubic --timepoints 1e2 1e5 --regularization-penalty 6.0 --knots 30 &
  cd ..
done



head unmasked.num.txt
1	44123
1	44124
1	57087
1	57088
1	57089
1	57090

head chr23.masked.txt 
23	1	757
23	763	2403
23	2406	5235
23	5242	5939
23	5982	6692

head unmasked.chr9.txt
9	1
9	2
9	3
9	4
9	5
9	6
9	7

head chr.A
chr1.fa chr2.fa chr7.fa chr8.fa chr13.fa chr14.fa chr19.fa chr20.fa chr25.fa chr26.fa chr31.fa chr32.fa chr37.fa chr38.fa 

cat ByChr.sh
for i in {1..42}
do
WGS --model file --type ByChr --chr $i --file unmasked.num.txt --out unmasked.chr$i &
  done


for i in 21 30
do
zcat /data1/home/yaozhou/Projects/EVO/data/merge/lineage/chr/chr$i.masked.bed.gz > chr$i.bed
cat chr$i.bed /data1/home/yaozhou/data/ref/wheat/gff/chr$i.masked.txt /data2/yaozhou/Projects/reference/chr$i.masked.txt > chr$i.masked.txt
sort -k2n -k3n chr$i.masked.txt > chr$i.sorted.txt
bedtools merge -i chr$i.sorted.txt -d 100 > chr$i.masked.bed &
  done


#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/sample_all/domemmer/vcf
cd /data1/home/xuebo/Projects/Speciation/smcpp/sample_all/domemmer/vcf
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
WGS --model vcf --type GenerateDiploid --file /data1/home/xuebo/Projects/Speciation/smcpp/E3removeExon/chr$i.vcf.gz --file2 /data1/home/xuebo/Projects/Speciation/smcpp/sample_all/group/dom_emmer_100.txt --out chr$i.vcf &
  done


