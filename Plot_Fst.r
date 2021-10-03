#画每组比较的boxplot图
library(stringr)
path <- "/data2/yafei/Project3/FST_group/VmapData/A/boxplot"
name <- c("5_10", "5_11", "5_12", "5_13", "5_14", "5_15", "5_16", "5_17", "5_18", "5_19", "5_20", "5_21", "5_22", "5_23", "5_24", "5_6", "5_7", "5_8", "5_9", "6_10", "6_11", "6_12", "6_13", "6_14", "6_15", "6_16", "6_17", "6_18", "6_19", "6_20", "6_21", "6_22", "6_23", "6_24", "6_7", "6_8", "6_9")
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 
length(data)
for(i in 1:length(data)){
  data[[i]]$Name <- name[i]
  data[[i]]$Before <-  str_split_fixed(data[[i]]$Name, "_", 2)[1]
  data[[i]]$After <- sapply(strsplit(as.character(data[[i]]$Name),'_'), "[", 2)
}

all <- data.frame()
for(i in 1:length(data)){
  all <- rbind(all, data[[i]])
}

library(ggplot2)
pdf(file = "FST_box_weight.pdf",width=200,height=100) #结果保存
ggplot(all, aes(After, WEIGHTED_FST, fill=factor(Before))) + geom_boxplot() +facet_grid(.~After, scales = "free_x")+
  scale_x_discrete()
dev.off()
pdf(file = "FST_box_mean.pdf",width=200,height=100) #结果保存
ggplot(all, aes(After, MEAN_FST, fill=factor(Before))) + geom_boxplot() +facet_grid(.~After, scales = "free_x")+
  scale_x_discrete()
dev.off()

#计算位点FST
#Alineage
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

for i in {14,15,16,17,18,19,20,21,22,23}
do
for((j=i+1;j<=24;j++))
  do
{
  read -u6
  {
    vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$j".txt"  --out /data2/yafei/Project3/FST_group/VmapData/A/site_fst/"p_"$i"_"$j >> Anohup$i 2>& 1
    echo "" >&6
  } & 
} 
done 
done
wait

exec 6>&-
  
  #Blineage
  #pop5-24:
  
  for i in {14,15,16,17,18,19,20,21,22,23}
do
for((j=i+1;j<=24;j++))
  do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --weir-fst-pop /data2/yafei/Project3/FST_group/subgroups/"pop"$j".txt" --out /data2/yafei/Project3/FST_group/VmapData/A/site_fst/"p_"$i"_"$j >> Anohup$i 2>& 1 &
  done 
done

#提取D site FST top 5%的位点
for i in `ls *.weir.fst`
do
sed '/-nan/d' $i  | awk  '{if ($3 < 0 ) {$3=0;print $0} else{print $0}}'  | sort -k3 -gr  > $i.txt &
  done
#计算每个文件的行数，并计算%5是多少行，重定向到*pos文件中
for i in `ls *txt`
do
wc -l $i
done

#提取fst大于0.2的位点
for i in `ls *txt`
do
awk '{if($3 > 0.2) {print $0} else {exit}}' $i  > $i.pos2 &
  done

#粘到Excel里

cat *pos2 | awk '{print $1"\t"$2}' | sort | uniq |awk '{print $1"\t"$2-1"\t"$2}' |sort -k1,1n -k2,2n > D_highfst.bed

awk '{if($1 == 10) print $2"\t"$3-1"\t"$3}' count.0.2thresh |sort -k1,1n -k2,2n > my.bed

#提取
yafei@203:/data2/yafei/Project3/FST_group/VmapData/D/site_fst/D_highfst.bed


for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
bedtools intersect -a /data1/home/yafei/Project3/Vmap1.1/chr${chr}.all.vcf.gz -b /data2/yafei/Project3/FST_group/VmapData/D/site_fst/my/myD.bed -header > chr${chr}.D.vcf &
  done