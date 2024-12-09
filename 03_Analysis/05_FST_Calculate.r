#------------------------------------------------------------------------计算FST.sh--------------------------------------------------------------------------------------
#FST_group: /data2/yafei/Project3/FST_group/subgroups/ 共27个groups。
#工作目录：/data2/yafei/Project3/FST_group/VmapData/
#AABB
for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38,3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
vcftools --vcf /data2/yafei/Project3/V+W_E2/E2_all/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/AABB_AABBDD_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W_E2/E2_maf/AABB_AABBDD/chr${chr}.mafAABB_AABBDD.vcf.gz 
done
vcf-concat *mafAABB_AABBDD.vcf.gz  > mafAABB_AABBDD.vcf 
bgzip -c mafAABB_AABBDD.vcf > mafAABB_AABBDD.vcf.gz

#pop1:	A_Wild_einkorn	
#pop2:	A_Domesticated_einkorn	
#pop3:	A_Urartu
#pop4:	B_Speltoides
#pop5:	AABB_Wild_emmer
#pop6:	AABB_Domesticated_emmer
#pop7:	AABB_Georgian_wheat
#pop8:	AABB_Ispahanicum
#pop9:	AABB_Rivet_wheat
#pop10:	AABB_Polish_wheat
#pop11:	AABB_Persian_wheat
#pop12:	AABB_Khorasan_wheat
#pop13:	AABB_Durum
#pop14:	AABBDD_Spelt
#pop15:	AABBDD_Macha
#pop16:	AABBDD_Club_wheat
#pop17:	AABBDD_Indian_dwarf_wheat
#pop18:	AABBDD_Yunan_wheat
#pop19:	AABBDD_Xinjiang_wheat
#pop20:	AABBDD_Tibetan_semi_wild
#pop21:	AABBDD_Vavilovii
#pop22:	AABBDD_Landrace_new
#pop23:	AABBDD_Cultivar
#pop24:	AABBDD_Synthetic
#pop25:	DD_Strangulata
#pop26:	DD_Meyeri
#pop27:	DD_Anathera

#A: /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz
#B: /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz
#D: /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz
#AABB: /data2/yafei/Project3/V+W_E2/E2_maf/AABB_AABBDD/mafAABB_AABBDD.vcf.gz
#AABBDD: /data2/yafei/Project3/V+W_E2/E2_maf/AABBDD/AABBDD.vcf.gz

#分lineage计算fst: /data2/yafei/Project3/FST_group/lineage/
#Alineage
#pop1:
for i in {2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/pop1.txt --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/A/"p_1_"$i 
done

#pop2:
for i in {3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/pop2.txt --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/A/"p_2_"$i &
  done
#pop3:
for i in {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/pop3.txt --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/A/"p_3_"$i &
  done

#pop5-24:
thread_num=20
tempfifo="my_temp_fifo"
mkfifo ${tempfifo}	
exec 6<>${tempfifo}
rm -f ${tempfifo}
for ((i=1;i<=${thread_num};i++))
  do
{
  echo 
}
done >&6 

for i in {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
do
for((j=i+1;j<=24;j++))
  do
{
  read -u6
  {
    vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$j".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/A/"p_"$i"_"$j 
    echo "" >&6
  } & 
} 
done 
done
wait

exec 6>&-
  
  
  #Blineage
  #pop5-24:
  
  thread_num=15
tempfifo="my_temp_fifo"
mkfifo ${tempfifo}	
exec 6<>${tempfifo}
rm -f ${tempfifo}
for ((i=1;i<=${thread_num};i++))
  do
{
  echo 
}
done >&6 

for i in {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
do
for((j=i+1;j<=24;j++))
  do
{
  read -u6
  {
    vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$j".txt" --fst-window-size 1000000 --out /data2/yafei/Project3/FST_group/VmapData/B/"p_"$i"_"$j > Bnohup$i 2>& 1
    echo "" >&6
  } & 
} 
done 
done

wait

exec 6>&-
  
  #Dlineage
  #pop14-27:
  thread_num=60
tempfifo="my_temp_fifo"
mkfifo ${tempfifo}	
exec 6<>${tempfifo}
rm -f ${tempfifo}
for ((i=1;i<=${thread_num};i++))
  do
{
  echo 
}
done >&6 
for i in {14,15,16,17,18,19,20,21,22,23,24,25,26}
do
for((j=i+1;j<=27;j++))
  do
{
  read -u6
  {
    vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$j".txt" --out /data2/yafei/Project3/FST_group/VmapData/D/site_fst/"p_"$i"_"$j >> Dnohup$i 2>& 1
    echo "" >&6
  } & 
} 
done 
done

wait
exec 6>&-
  
  #grep "out p" nohup.out | awk '{print $2}'
  #grep "weight" nohup.out | awk '{print $7}'
  
  #A/B/D.fst 修改输出文件名
  for i in `ls *windowed.weir.fst`
do
awk 'BEGIN {count = 0} {if($5<10){sum+=$5; count++}} END {print FILENAME,"Average = ", sum/count}' $i >> A_weight.fst
done

#求中位数
for i in `ls *windowed.weir.fst`
do
echo $i
awk '{print $5}' $i |sed '1d' |sort |awk -f medium.awk
done |xargs -n 2 > D_medium.fst

#AB/ABD.fst
for i in `cat ABfiles`
do
echo $i 
cat /data2/yafei/Project3/FST_group/VmapData/A/$i  /data2/yafei/Project3/FST_group/VmapData/B/$i | sed '/WEIGHTED/d'|awk 'BEGIN {count = 0} {if($5<10){sum+=$5; count++}} END {print "Average = ", sum/count}' >> AB_weight.fst 
done 

#画热图 修改文件名，数据集名，输出文件名
#yafei@203:/data2/yafei/Project3/FST_group/VmapData/heatmap

library("corrplot")
data <- read.table("D_weight.fst_ibs",header=T,stringsAsFactors=F,sep="\t")
rownames(data) <- data$Dlineage
data<- data[,-1]
data<- as.matrix(data)
pdf("fst_ibs_D.pdf",width=30,height=30)
corrplot(data,method = "color",tl.col="black",tl.srt = 45,addrect=4,addCoef.col = "grey",number.cex=4,number.digits=3,tl.cex=2.5,cl.cex=3,cl.lim = c(0, 1))
dev.off()
