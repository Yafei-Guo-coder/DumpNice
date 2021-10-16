#!/bin/bash
sleep 4800s
for i in {1..42}
do
  WGS --model vcf --type toXPCLR --group1 ../groupTaxa/EA.txt --rec ../../slidewindow_recomrate_updown_20M.txt --chr $i --file /data2/xuebo/Projects/Speciation/E6/chr${i}.Land.vcf.gz --out groupEAChr${i} &
done

#!/bin/bash
for i in {1..42}
do
	XPCLR -xpclr ../groupEA/groupEAChr${i}.geno ../groupWA/groupWAChr${i}.geno ../groupWA/groupWAChr${i}.snp EA_WA_10kchr${i} -w1 0.005 500 10000 $i -p1 0.95 &
done
