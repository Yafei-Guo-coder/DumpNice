#-----------------------------------------------------Landrace & Wild/Domesticated emmer / Free IBS--------------------------------------------------------------------------------------
#工作目录：yafei@203:/data2/yafei/Project3/IBS
#画Landrace & Wild emmer的IBS
#Landrace包括5种：Landrace/Club Wheat/Indian Dwarf/Yunan Wheat/Vovilovi 
#修改 输入文件名
#读取group名的文件

path <- "/data2/yafei/Project3/group/subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)
A <- read.table("/data2/yafei/Project3/IBS/spelt/Alineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
A$Lineage <- NA
A$Fill <- NA
nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_Durum","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_Khorasan_wheat","AABBDD_Landrace","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_Polish_wheat","AABB_Rivet_wheat","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
#A:487  485-487
#B:406  404-406
#D:306  304-306
1.SS_taxa.txt
2.sub_Anathera.txt
3.sub_Club_wheat.txt
4.sub_Cultivar.txt
sub_Domesticated_einkorn.txt
sub_Domesticated_emmer.txt
sub_Durum.txt
sub_Georgian_wheat.txt
sub_Indian_dwarf_wheat.txt
sub_Ispahanicum.txt
sub_Khorasan_wheat.txt
sub_Landrace.txt
sub_Landrace_new.txt
sub_Macha.txt
sub_Meyeri.txt
sub_Persian_wheat.txt
sub_Polish_wheat.txt
sub_Rivet_wheat.txt
sub_Spelt.txt
sub_Speltoides.txt
sub_Strangulata.txt
sub_Synthetic.txt
sub_Tibetan_semi_wild.txt
sub_Urartu.txt
sub_Vavilovii.txt
sub_Wild_einkorn.txt
sub_Wild_emmer.txt
sub_Xinjiang_wheat.txt
sub_Yunan_wheat.txt
#A
for(i in 1:length(data)){
  A[which(A$Dxy %in% data[[i]][,1]),486] <- nameABD[i]
  if(i %in% c(3,4,9,12,13,14,19,22,23,25,28,29)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "pink"
  }else if(i %in% c(5,24,26)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "lightblue"
  }else if(i %in% c(6,7,8,10,11,16,17,18,27)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "orange"
  }else if(i %in% c(1,20)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "red"
  }else if(i %in% c(2,15,21)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "yellow"
  }
}
rownames(A) <- A$Dxy

#B
for(i in 1:length(data)){
  A[which(A$Dxy %in% data[[i]][,1]),405] <- nameABD[i]
  if(i %in% c(3,4,9,12,13,14,19,22,23,25,28,29)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "pink"
  }else if(i %in% c(5,24,26)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "lightblue"
  }else if(i %in% c(6,7,8,10,11,16,17,18,27)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "orange"
  }else if(i %in% c(1,20)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "red"
  }else if(i %in% c(2,15,21)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "yellow"
  }
}
rownames(A) <- A$Dxy

#D
for(i in 1:length(data)){
  A[which(A$Dxy %in% data[[i]][,1]),305] <- nameABD[i]
  if(i %in% c(3,4,9,12,13,14,19,22,23,25,28,29)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "pink"
  }else if(i %in% c(5,24,26)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "lightblue"
  }else if(i %in% c(6,7,8,10,11,16,17,18,27)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "orange"
  }else if(i %in% c(1,20)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "red"
  }else if(i %in% c(2,15,21)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "yellow"
  }
}
rownames(A) <- A$Dxy

#Land <- A[which(A$Lineage == "AABBDD_Landrace" | A$Lineage == "AABBDD_Club_wheat" | A$Lineage == "AABBDD_Indian_dwarf_wheat" | A$Lineage == "AABBDD_Yunan_wheat" | A$Lineage == "AABBDD_Vavilovii"),]
Land <- A[which(A$Lineage == "AABBDD_Landrace"),]
wildname <- c("B023","B024","B025","B026","B027","B028","B029","B030","B031","B032","B033","B034","B035","B036","B037","B038","B039","B040","B041","B042","B044","B045","B046","B047","B048","B049","B050","B051","B052","TWS04","XI_S12","XI_S13","XI_S14","XI_S15","XI_S16","XI_S17","XI_S18","XI_S19","XI_S20","XI_S21","XI_W1","XI_W2","XI_W3","XI_W4","XI_W5","XI_W6","XI_W7","XI_W8","XI_W9","XI_W10")
Wild <- Land[,wildname]
Wild$mean <- apply(Wild,1,mean)
dim(Wild)

Domesname <- c("B063","B064","B065","B066","B067","B068","B069","B070","B071","B072","B073","B074","B075","B076","B077","B078","B079","B080","B081","B082","B083","B084","B085","B086","B087","B088","B089","B090","B091")
Domes <- Land[,Domesname]
Domes$mean <- apply(Domes,1,mean)
dim(Domes)

Freename <- c("B001","B002","B003","B004","B005","B006","B007","B008","B009","B010","B011","B012","B013","B014","B015","B016","B017","B018","B019","B020","B021","B022","B053","B054","B055","B056","B057","B058","B059","B060","B061","B062","B093","B094","B096","B097","B098","B099","B100","B101","B102","B113","B114","B115","B116","B117","B118","B119","B120","B121","B122","B123","B124","B125","XI_S1","XI_S2","XI_S3","XI_S4","XI_S5")
Free <- Land[,Freename]
Free$mean <- apply(Free,1,mean)
dim(Free)

Land <- A[which(A$Lineage == "AABBDD_Landrace"),]
Strangname <- c("D001","D010","D014","D015","D016","D022","D023","D026","D028","D031","D032","D033","D034","D035","D036","D037","D038","D039","D041","D042","D043","D044","D046","D047","D048","D049","D050","D051","D052","D054","D055","D056","XI_A2","XI_A3","XI_A4","XI_A5")
Strang <- Land[,Strangname]
Strang$mean <- apply(Strang,1,mean)
dim(Strang)

#对Landrace做admixture structure
nohup vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --keep /data2/yafei/Project3/IBS/145Landrace --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/IBS/A_Landrace.vcf.gz  > outA &
  nohup vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz --keep /data2/yafei/Project3/IBS/145Landrace --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/IBS/B_Landrace.vcf.gz > outB &
  nohup vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz --keep /data2/yafei/Project3/IBS/145Landrace --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/IBS/D_Landrace.vcf.gz > outC &
  
  plink --vcf D_Landrace.vcf.gz --indep-pairwise 50 10 0.2 --out D.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf D_Landrace.vcf.gz --make-bed --extract D.filterLD.prune.in --out  D.filter.prune.in --double-id --autosome-num 42
plink --bfile D.filter.prune.in --recode 12 --out D.filter.prune.in --autosome-num 42

#A_wild_land
nohup bash A_wild_land.sh > A_wild_land.log &
  #A_dome_land
  nohup bash A_dome_land.sh > A_dome_land.log &
  #A_free_land
  nohup bash A_free_land.sh > A_free_land.log &
  #B_wild_land
  nohup bash B_wild_land.sh > B_wild_land.log &
  #B_dome_land
  nohup bash B_dome_land.sh > B_dome_land.log &
  #B_free_land
  nohup bash B_free_land.sh > B_free_land.log &
  #D_strang_land
  nohup bash D_strang_land.sh > D_strang_land.log &
  #AB_wild_land
  nohup vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/ABlineage.vcf.gz --keep /data2/yafei/Project3/IBS/wild_land --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/IBS/AB_wild_land.vcf.gz  > AB_wild_land.log &
  #AB_dome_land
  nohup vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/ABlineage.vcf.gz --keep /data2/yafei/Project3/IBS/dome_land --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/IBS/AB_dome_land.vcf.gz  > AB_dome_land.log &
  #AB_free_land
  nohup vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/ABlineage.vcf.gz --keep /data2/yafei/Project3/IBS/free_land --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/IBS/AB_free_land.vcf.gz  > AB_free_land.log &
  
  #A_wild_land
  plink --vcf A_wild_land.vcf.gz --indep-pairwise 50 10 0.2 --out A_wild_land.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf A_wild_land.vcf.gz --make-bed --extract A_wild_land.filterLD.prune.in --out  A_wild_land.prune.in --double-id --autosome-num 42
plink --bfile A_wild_land.prune.in --recode 12 --out A_wild_land.prune.in  --autosome-num 42
#A_dome_land
plink --vcf A_dome_land.vcf.gz --indep-pairwise 50 10 0.2 --out A_dome_land.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf A_dome_land.vcf.gz --make-bed --extract A_dome_land.filterLD.prune.in --out  A_dome_land.prune.in --double-id --autosome-num 42
plink --bfile A_dome_land.prune.in --recode 12 --out A_dome_land.prune.in  --autosome-num 42
#A_free_land
plink --vcf A_free_land.vcf.gz --indep-pairwise 50 10 0.2 --out A_free_land.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf A_free_land.vcf.gz --make-bed --extract A_free_land.filterLD.prune.in --out  A_free_land.prune.in --double-id --autosome-num 42
plink --bfile A_free_land.prune.in --recode 12 --out A_free_land.prune.in  --autosome-num 42
#B_wild_land
plink --vcf B_wild_land.vcf.gz --indep-pairwise 50 10 0.2 --out B_wild_land.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf B_wild_land.vcf.gz --make-bed --extract B_wild_land.filterLD.prune.in --out  B_wild_land.prune.in --double-id --autosome-num 42
plink --bfile B_wild_land.prune.in --recode 12 --out B_wild_land.prune.in  --autosome-num 42
#B_dome_land
plink --vcf B_dome_land.vcf.gz --indep-pairwise 50 10 0.2 --out B_dome_land.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf B_dome_land.vcf.gz --make-bed --extract B_dome_land.filterLD.prune.in --out  B_dome_land.prune.in --double-id --autosome-num 42
plink --bfile B_dome_land.prune.in --recode 12 --out B_dome_land.prune.in  --autosome-num 42
#B_free_land
plink --vcf B_free_land.vcf.gz --indep-pairwise 50 10 0.2 --out B_free_land.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf B_free_land.vcf.gz --make-bed --extract B_free_land.filterLD.prune.in --out  B_free_land.prune.in --double-id --autosome-num 42
plink --bfile B_free_land.prune.in --recode 12 --out B_free_land.prune.in  --autosome-num 42
#D_strang_land
plink --vcf D_strang_land.vcf.gz --indep-pairwise 50 10 0.2 --out D_strang_land.filterLD --allow-extra-chr --double-id --autosome-num 42
plink --vcf D_strang_land.vcf.gz --make-bed --extract D_strang_land.filterLD.prune.in --out  D_strang_land.prune.in --double-id --autosome-num 42
plink --bfile D_strang_land.prune.in --recode 12 --out D_strang_land.prune.in  --autosome-num 42

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
for i in {1..10}  
do
{
  read -u6
  {
    admixture --cv A_wild_land.prune.in.ped $i >> logA_wild_land.txt
    admixture --cv A_dome_land.prune.in.ped $i >> logA_dome_land.txt
    admixture --cv A_free_land.prune.in.ped $i >> logA_free_land.txt
    admixture --cv B_wild_land.prune.in.ped $i >> logB_wild_land.txt
    admixture --cv B_dome_land.prune.in.ped $i >> logB_dome_land.txt
    admixture --cv B_free_land.prune.in.ped $i >> logB_free_land.txt
    admixture --cv D_strang_land.prune.in.ped $i >> logD_strang_land.txt
  } & 
} 
done 
wait
exec 6>&-
  
  grep "CV error" log.txt 

for i in `cat Lulab4T_59.txt`
do
start_time=`date +%s`
samtools sort -n -m 8G -o /data4/home/yafei/vmap3/bam/${i}.namesort.bam -O bam -@ 4 /data4/home/yafei/vmap3/bam/${i}.bam && samtools fixmate -@ 4 -m /data4/home/yafei/vmap3/bam/${i}.namesort.bam /data4/home/yafei/vmap3/bam/${i}.fixmate.bam && samtools sort -m 8G -o /data4/home/yafei/vmap3/bam/${i}.fixmate.pos.bam -O bam -@ 4 /data4/home/yafei/vmap3/bam/${i}.fixmate.bam && rm -f /data4/home/yafei/vmap3/bam/${i}.namesort.bam && samtools markdup -@ 4 -r /data4/home/yafei/vmap3/bam/${i}.fixmate.pos.bam /data4/home/yafei/vmap3/bam/${i}.rmdup.bam && rm -f /data4/home/yafei/vmap3/bam/${i}.fixmate.bam && rm -f /data4/home/yafei/vmap3/bam/${i}.fixmate.pos.bam && samtools index /data4/home/yafei/vmap3/bam/${i}.rmdup.bam
stop_time=`date +%s`
echo "TIME:`expr $stop_time - $start_time`"
done

#计算IBS
#工作目录：203:/data2/yafei/Project3/CS_Vmap/CS/
nohup java -jar /data1/home/yafei/Project3/Vmap1.1/C41_getIBS_distance2.jar --file1 /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --out Alineage.all.ibs.txt > logA.txt &
  nohup java -jar /data1/home/yafei/Project3/Vmap1.1/C41_getIBS_distance2.jar --file1 /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz --out Blineage.all.ibs.txt > logB.txt &
  nohup java -jar /data1/home/yafei/Project3/Vmap1.1/C41_getIBS_distance2.jar --file1 /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz --out Dlineage.all.ibs.txt > logD.txt &
  
  #计算spelt和domesticated emmer的IBS（A）
  
  setwd("/data2/yafei/Project3/IBS/spelt")
path <- "/data2/yafei/Project3/group/subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)
A <- read.table("/data2/yafei/Project3/IBS/Dlineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
A$Lineage <- NA
A$Fill <- NA
nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_Durum","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_Khorasan_wheat","AABBDD_Landrace","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_Polish_wheat","AABB_Rivet_wheat","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
#A:487  485-487
#B:406  404-406
#D:306  304-306
#A
for(i in 1:length(data)){
  A[which(A$Dxy %in% data[[i]][,1]),486] <- nameABD[i]
  if(i %in% c(3,4,9,12,13,14,19,22,23,25,28,29)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "pink"
  }else if(i %in% c(5,24,26)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "lightblue"
  }else if(i %in% c(6,7,8,10,11,16,17,18,27)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "orange"
  }else if(i %in% c(1,20)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "red"
  }else if(i %in% c(2,15,21)){
    A[which(A$Dxy %in% data[[i]][,1]),487] <- "yellow"
  }
}
rownames(A) <- A$Dxy

#Land <- A[which(A$Lineage == "AABBDD_Landrace" | A$Lineage == "AABBDD_Club_wheat" | A$Lineage == "AABBDD_Indian_dwarf_wheat" | A$Lineage == "AABBDD_Yunan_wheat" | A$Lineage == "AABBDD_Vavilovii"),]
Domesticated_emmer <- rownames(A[which(A$Lineage == "AABB_Domesticated_emmer"),])
Land <- rownames(A[which(A$Lineage == "AABBDD_Landrace" | A$Lineage == "AABBDD_Club_wheat" | A$Lineage == "AABBDD_Indian_dwarf_wheat" | A$Lineage == "AABBDD_Yunan_wheat" | A$Lineage == "AABBDD_Vavilovii"),])
Wild_emmer <- rownames(A[which(A$Lineage == "AABB_Wild_emmer"),])
Free <- rownames(A[which(A$Lineage == "AABB_Durum" | A$Lineage == "AABB_Khorasan_wheat" | A$Lineage == "AABB_Polish_wheat" |A$Lineage == "AABB_Rivet_wheat"),])
Spelt <- rownames(A[which(A$Lineage == "AABBDD_Spelt"),])

Xinjiang <- rownames(A[which(A$Lineage == "AABBDD_Xinjiang_wheat"),])
Macha <- rownames(A[which(A$Lineage == "AABBDD_Macha"),])
Persian <- rownames(A[which(A$Lineage == "AABB_Persian_wheat"),])

S <- A[c(Domesticated_emmer,Land),Spelt]
S$mean <- apply(S,1,mean)
dim(S)

X <- A[c(Free, Land), Xinjiang]
X$mean <- apply(X,1,mean)

M <- A[c(Free,Wild_emmer, Land), Macha]
M <- A[c(Domesticated_emmer,Free,Wild_emmer, Land), Macha]

M$mean <- apply(M,1,mean)

P <- A[c(Free, Land), Persian]
P$mean <- apply(P,1,mean)

#B
for(i in 1:length(data)){
  A[which(A$Dxy %in% data[[i]][,1]),405] <- nameABD[i]
  if(i %in% c(3,4,9,12,13,14,19,22,23,25,28,29)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "pink"
  }else if(i %in% c(5,24,26)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "lightblue"
  }else if(i %in% c(6,7,8,10,11,16,17,18,27)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "orange"
  }else if(i %in% c(1,20)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "red"
  }else if(i %in% c(2,15,21)){
    A[which(A$Dxy %in% data[[i]][,1]),406] <- "yellow"
  }
}
rownames(A) <- A$Dxy
#Land <- A[which(A$Lineage == "AABBDD_Landrace" | A$Lineage == "AABBDD_Club_wheat" | A$Lineage == "AABBDD_Indian_dwarf_wheat" | A$Lineage == "AABBDD_Yunan_wheat" | A$Lineage == "AABBDD_Vavilovii"),]
Domesticated_emmer <- rownames(A[which(A$Lineage == "AABB_Domesticated_emmer"),])
Land <- rownames(A[which(A$Lineage == "AABBDD_Landrace" | A$Lineage == "AABBDD_Club_wheat" | A$Lineage == "AABBDD_Indian_dwarf_wheat" | A$Lineage == "AABBDD_Yunan_wheat" | A$Lineage == "AABBDD_Vavilovii"),])

Spelt <- c("TW011","TW012","TW013","TW014","TW015","TW016","TW017","TW018","TW019","TW020","TW021","TW022","TW023","TW024","TWA01","TWA02","TWA03","TWA04","TWA05","TWA06","TWA07","TWA08","TWA09","TWA10")

S <- A[c(Domesticated_emmer,Land),Spelt]
S$mean <- apply(S,1,mean)
dim(S)
#D
for(i in 1:length(data)){
  A[which(A$Dxy %in% data[[i]][,1]),305] <- nameABD[i]
  if(i %in% c(3,4,9,12,13,14,19,22,23,25,28,29)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "pink"
  }else if(i %in% c(5,24,26)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "lightblue"
  }else if(i %in% c(6,7,8,10,11,16,17,18,27)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "orange"
  }else if(i %in% c(1,20)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "red"
  }else if(i %in% c(2,15,21)){
    A[which(A$Dxy %in% data[[i]][,1]),306] <- "yellow"
  }
}

rownames(A) <- A$Dxy
Land <- rownames(A[which(A$Lineage == "AABBDD_Landrace" | A$Lineage == "AABBDD_Club_wheat" | A$Lineage == "AABBDD_Indian_dwarf_wheat" | A$Lineage == "AABBDD_Yunan_wheat" | A$Lineage == "AABBDD_Vavilovii"),])

Spelt <- c("TW011","TW012","TW013","TW014","TW015","TW016","TW017","TW018","TW019","TW020","TW021","TW022","TW023","TW024","TWA01","TWA02","TWA03","TWA04","TWA05","TWA06","TWA07","TWA08","TWA09","TWA10")

S <- A[c(Land),Spelt]
S$mean <- apply(S,1,mean)
dim(S)


