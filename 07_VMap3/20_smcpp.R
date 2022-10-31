#split fold
-rw-rw-r--. 1 xuebo xuebo  467 Jul 24 15:43 getdiploidvcf.sh
-rw-rw-r--. 1 xuebo xuebo 1.2K Jul 24 16:09 getsmcformat.sh
-rw-rw-r--. 1 xuebo xuebo  650 Jul 24 16:23 getcp.sh
-rw-rw-r--. 1 xuebo xuebo  617 Jul 24 16:29 split_1.sh
-rw-rw-r--. 1 xuebo xuebo  617 Jul 24 16:31 split_2.sh
-rw-rw-r--. 1 xuebo xuebo  617 Jul 24 16:32 split_3.sh
-rw-rw-r--. 1 xuebo xuebo  626 Jul 24 16:36 split_4.sh
-rw-rw-r--. 1 xuebo xuebo  837 Jul 25 11:54 split_bootstrap20.sh
-rw-rw-r--. 1 xuebo xuebo  837 Jul 25 11:58 split_bootstrap20_2.sh
-rw-rw-r--. 1 xuebo xuebo  837 Jul 25 11:59 split_bootstrap20_3.sh
-rw-rw-r--. 1 xuebo xuebo  846 Jul 25 12:01 split_bootstrap20_4.sh
-rw-rw-r--. 1 xuebo xuebo  830 Jul 25 17:13 split_bootstrap20_5.sh
-rw-rw-r--. 1 xuebo xuebo  255 Jul 26 10:58 plot.sh

#Ne fold
#ls -lh *sh |sort -k7
-rw-rw-r--. 1 xuebo xuebo 513 Jul 12 13:19 getdiploidvcf.sh
-rw-rw-r--. 1 xuebo xuebo 844 Jul 12 13:41 getsmcformat.sh
-rw-rw-r--. 1 xuebo xuebo 471 Jul 12 13:45 estimate.sh
-rw-rw-r--. 1 xuebo xuebo 575 Jul 12 14:08 estimate_V_65_8.sh
-rw-rw-r--. 1 xuebo xuebo 579 Jul 12 14:11 estimate_V_65_10.sh
-rw-rw-r--. 1 xuebo xuebo 569 Jul 12 14:13 estimate_V_1_7.sh
-rw-rw-r--. 1 xuebo xuebo 570 Jul 12 14:16 estimate_V_1_8.sh
-rw-rw-r--. 1 xuebo xuebo 570 Jul 12 14:17 estimate_V_1_9.sh
-rw-rw-r--. 1 xuebo xuebo 573 Jul 12 14:19 estimate_V_1_10.sh
-rw-rw-r--. 1 xuebo xuebo 574 Jul 12 14:21 estimate_V_1_11.sh
-rw-rw-r--. 1 xuebo xuebo 603 Jul 19 20:44 estimate_V2_e_8.sh
-rw-rw-r--. 1 xuebo xuebo 582 Jul 19 20:49 estimate_V2_e_9.sh
-rw-rw-r--. 1 xuebo xuebo 585 Jul 19 22:29 estimate_V2_6.5e-9.sh

#getsmcformat.sh 
#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V1
cd /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V1
sed 's/\t/_/g' /data1/home/xuebo/Projects/Speciation/smcpp/sample15_2/group/wild_emmer_15.txt |  tr "\n" "," | sed s'/.$//' > group_wildemmer_header.txt &
  mkdir data
name=$(sed -n 1p group_wildemmer_header.txt)
echo $name
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
smc++ vcf2smc --mask /data1/home/xuebo/Projects/Evo/smcpp/maskfile/mask1/chr$i.masked.bed.gz /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/vcf/chr$i.vcf.gz \
/data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V1/data/chr$i.smc.gz $i wildemmer:$name  &
  done




#getdiploidvcf.sh
#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/vcf
cd /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/vcf
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
WGS --model vcf --type GenerateDiploid --file /data1/home/xuebo/Projects/Speciation/smcpp/E3removeExon/chr$i.vcf.gz --file2 /data1/home/xuebo/Projects/Speciation/smcpp/sample15_2/group/wild_emmer_15.txt --out chr$i.vcf &
  done


(base) [xuebo@lulab2 wildemmer]$ cat estimate_V2_e-9.sh
#!/bin/bash
for i in {"7e-9","8e-9","9e-9"}
do
mkdir /data2/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V2_${i}
mkdir /data2/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V2_${i}/analysis
OMP_NUM_THREADS=1 smc++ estimate ${i} /data2/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V1/data/chr*.smc.gz --polarization-error 0.5 \
-o /data2/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V2_${i}/analysis --cores 30 --spline cubic --timepoints 1e2 1e5 --regularization-penalty 6.0 --knots 30 &
  done



#split.sh
#!/bin/bash
for i in {"1e-9","2e-9","3e-9"}
do
mkdri /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/split/split_${i}
OMP_NUM_THREADS=1 smc++ split -o /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/split/split_${i} /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V2_${i}/analysis/model.final.json \
/data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/domemmer/V2_${i}/analysis/model.final.json /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/split/V1/data/chr*.gz --cores 32 &
  done

#plot.sh
#!/bin/bash
for i in {"1e-9","2e-9","3e-9","4e-9","5e-9","6e-9","7e-9","8e-9","9e-9","6.5e-9","1e-8","2e-8","3e-8"}
do
for j in {1..20}
do
smc++ plot ./split_${i}/split_wildemmer_domemmer_${j}.pdf ./split_${i}/split${j}/model.final.json -c
done
done

#plot.sh 
#!/bin/bash
for i in {"1e-9","2e-9","3e-9","4e-9","5e-9","6e-9","7e-9","8e-9","9e-9","6.5e-9","1e-8","2e-8","3e-8"}
do
for j in {1..20}
do
smc++ plot ./split_${i}/split_wildemmer_domemmer_${j}.pdf ./split_${i}/split${j}/model.final.json -c
done
done
#plot_dist.R 
for(name in c("1e-9","2e-9","3e-9","4e-9","5e-9","6e-9","7e-9","8e-9","9e-9","6.5e-9","1e-8","2e-8","3e-8")){
  setwd(paste("/data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/split/split_",name,sep = ""))
  minx_matrix = data.frame(matrix(NA, ncol = 3, nrow = 20))
  colnames(minx_matrix) = c("Name","Bootstrap","Split_time")
  for(i in 1:20){
    dat1 = read.csv(paste("split_wildemmer_domemmer_",i,".csv",sep=""),head=T)
    minx = dat1[length(dat1$x),]$x
    minx_matrix[i,1] = paste("mu_",name,sep = "")
    minx_matrix[i,2] = paste("bootstrap",i,sep = "")
    minx_matrix[i,3] = minx
  }
  write.table(minx_matrix,"dis_matrix_split_wildemmer_domemmer.txt",col.names = F,row.names = F,quote=F,sep="\t")
}


#getcp.sh 
#!/bin/bash
mkdir  /data1/home/xuebo/Projects/Speciation/smcpp/split_sample15_speciation/split_wildemmer_domemmer/V1/split
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
cp /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/wildemmer/V1/data/chr$i.smc.gz  /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/split/V1/data/chr$i.smcpop1.gz
cp /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/domemmer/V1/data/chr$i.smc.gz  /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/split/V1/data/chr$i.smcpop2.gz
done


#-------------------->
cat getdiploidvcf.sh
#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/vcf
cd /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/vcf
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
WGS --model vcf --type GenerateDiploid --file /data1/home/xuebo/Projects/Speciation/smcpp/E3removeExon/chr$i.vcf.gz --file2 /data1/home/xuebo/Projects/Speciation/smcpp/sample15_2/group/Durum_15.txt --out chr$i.vcf &
  done
(base) [xuebo@lulab2 Durum]$ cat getsmcformat.sh
#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/V1
cd /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/V1
sed 's/\t/_/g' /data1/home/xuebo/Projects/Speciation/smcpp/sample15_2/group/Durum_15.txt |  tr "\n" "," | sed s'/.$//' > group_Durum_header.txt &
  mkdir data
name=$(sed -n 1p group_Durum_header.txt)
echo $name
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
smc++ vcf2smc --mask /data1/home/xuebo/Projects/Evo/smcpp/maskfile/mask1/chr$i.masked.bed.gz /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/vcf/chr$i.vcf.gz \
/data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/V1/data/chr$i.smc.gz $i Durum:$name  &
  done
(base) [xuebo@lulab2 Durum]$ cat estimate_V1_e_8.sh
#!/bin/bash
for i in {"1e-8","2e-8","3e-8"}
do
mkdir /data2/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/V1_${i}
mkdir /data2/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/V1_${i}/analysis
OMP_NUM_THREADS=1 smc++ estimate ${i} /data2/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/V1/data/chr*.smc.gz --polarization-error 0.5 \
-o /data2/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/SplitAB/Durum/V1_${i}/analysis --cores 30 --spline cubic --timepoints 1e2 1e5 --regularization-penalty 6.0 --knots 30 &
done

### getdiploidvcf.sh
#!/bin/bash
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/vcf
cd /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/vcf
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
WGS --model vcf --type GenerateDiploid --file /data1/home/xuebo/Projects/Speciation/smcpp/E3removeExon/chr$i.vcf.gz --file2 /data1/home/xuebo/Projects/Speciation/smcpp/sample15_2/group/landrace_15.txt --out chr$i.vcf &
  
  ### getsmcformat.sh
  #!/bin/bash
  mkdir /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/V1
cd /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/V1
sed 's/\t/_/g' /data1/home/xuebo/Projects/Speciation/smcpp/sample15_2/group/landrace_15.txt |  tr "\n" "," | sed s'/.$//' > group_landrace_AB_header.txt &
  mkdir data
name=$(sed -n 1p group_landrace_AB_header.txt)
echo $name
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
smc++ vcf2smc --mask /data1/home/xuebo/Projects/Evo/smcpp/maskfile/mask1/chr$i.masked.bed.gz /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/vcf/chr$i.vcf.gz \
/data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/V1/data/chr$i.smc.gz $i landrace_AB:$name  &
  done

### estimate.sh
#!/bin/bash
for i in {"1e-8","2e-8","3e-8"}
do
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/V2_${i}
mkdir /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/V2_${i}/analysis
OMP_NUM_THREADS=1 smc++ estimate ${i} /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/V1/data/chr*.smc.gz --polarization-error 0.5 \
-o /data1/home/xuebo/Projects/Speciation/smcpp/NP_edit/mutation_rate/M_65_9/landrace_AB/V2_${i}/analysis --cores 30 --spline cubic --timepoints 1e2 1e5 --regularization-penalty 6.0 --knots 30 &
  done


