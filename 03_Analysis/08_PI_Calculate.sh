#-------------------------------------------------------------计算 pi .sh--------------------------------------------------------------------------------------

#分lineage计算pi: 
#Alineage
for i in {1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do 
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --site-pi --out /data2/yafei/Project3/Pi/A/"pop"$i"1M.group" 
done

#Blineage
for i in {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --site-pi --out /data2/yafei/Project3/Pi/B/"pop"$i"1M.group"
done

#Dlineage
for i in {14,15,16,17,18,19,20,21,22,23,24,25,26,27}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --site-pi --out /data2/yafei/Project3/Pi/D/"pop"$i"1M.group"
done

#分染色体计算Windows pi

#Alineage:{1,2,7,8,13,14,19,20,25,26,31,32,37,38}
#Blineage:{3,4,9,10,15,16,21,22,27,28,33,34,39,40}
#Dlineage:{5,6,11,12,17,18,23,24,29,30,35,36,41,42}

for i in `ls *group.sites.pi`
do
for j in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
awk '{if($1=='"$j"') print $0}'  ${i} >> chr/${j}_${i}
done
done
#/data2/yafei/Project3/Pi/A/chr
#修改文件表头
for i in `ls *group.sites.pi`
do
sed -i '1i CHROM\tPOS\tPI' ${i}
done

#由于gcc版本的原因，转移到204上做
#/data1/home/yafei/Project3/Pi/A/chr
#计算windows pi
#Alineage name
#Blineage name
#Dlineage name

for i in `cat D.name`
do
for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/home/yafei/Project3/Pi/D/chr/${chr}_${i} --file2 /data1/home/yafei/Project3/Pi/SNPbed_1M/chr${chr}_1M.bed --out /data1/home/yafei/Project3/Pi/D/1M_window/${chr}_${i} &
  done
done

#把chr Windows pi合并成 群体 pi

for i in `cat D.name`
do
cat D/1M_window/ *_${i} >> D/by_chr/${i}
done

cat A.name B.name |sort |uniq -c |grep "2 " |awk '{print $2}' > ABfiles
cat A.name B.name D.name|sort |uniq -c |grep "3 " |awk '{print $2}' > ABDfiles

#ABlineage
for i in `cat ABfiles`
do
cat /data1/home/yafei/Project3/Pi/A/by_chr/$i  /data1/home/yafei/Project3/Pi/B/by_chr/$i  > /data1/home/yafei/Project3/Pi/AB/$i
done
#ABDlineage
for i in `cat ABDfiles`
do
cat /data1/home/yafei/Project3/Pi/A/by_chr/$i  /data1/home/yafei/Project3/Pi/B/by_chr/$i  /data1/home/yafei/Project3/Pi/D/by_chr/$i > /data1/home/yafei/Project3/Pi/ABD/$i
done

#重新计算新疆的PI值
#目录xuebo@204:/data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf  pop19.txt

nohup vcf-concat chr01.vcf chr02.vcf chr07.vcf chr08.vcf chr013.vcf chr014.vcf chr019.vcf chr020.vcf chr025.vcf chr026.vcf chr031.vcf chr032.vcf chr037.vcf chr038.vcf > Alineage.vcf &
  nohup vcf-concat chr03.vcf chr04.vcf chr09.vcf chr010.vcf chr015.vcf chr016.vcf chr021.vcf chr022.vcf chr027.vcf chr028.vcf chr033.vcf chr034.vcf chr039.vcf chr040.vcf > Blineage.vcf &
  nohup vcf-concat chr05.vcf chr06.vcf chr011.vcf chr012.vcf chr017.vcf chr018.vcf chr023.vcf chr024.vcf chr029.vcf chr030.vcf chr035.vcf chr036.vcf chr041.vcf chr042.vcf > Dlineage.vcf &
  bgzip -c ${i} > ${i}.gz
nohup tabix -p vcf Alineage.vcf.gz &
  nohup tabix -p vcf Blineage.vcf.gz &
  nohup tabix -p vcf Dlineage.vcf.gz &
  #分染色体计算PI
  for chr in {1..42}
do
nohup vcftools --vcf /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/chr0${chr}.vcf --keep /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/pop19.txt --site-pi --out /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/xinjiang/${chr}_pop191M.group.sites.pi > nohupVmap 2>& 1 &
  done
#转移到yafei204对应的位置:/data1/home/yafei/Project3/Pi/xinjiang

for chr in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/home/yafei/Project3/Pi/xinjiang/${chr}_pop191M.group.sites.pi.sites.pi --file2 /data1/home/yafei/Project3/Pi/SNPbed_1M/chr${chr}_1M.bed --out Awindow/${chr}_pop191M.group.sites.pi &
  done


for chr in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/home/yafei/Project3/Pi/xinjiang/${chr}_pop191M.group.sites.pi.sites.pi --file2 /data1/home/yafei/Project3/Pi/SNPbed_1M/chr${chr}_1M.bed --out Bwindow/${chr}_pop191M.group.sites.pi &
  done

for chr in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/home/yafei/Project3/Pi/xinjiang/${chr}_pop191M.group.sites.pi.sites.pi --file2 /data1/home/yafei/Project3/Pi/SNPbed_1M/chr${chr}_1M.bed --out Dwindow/${chr}_pop191M.group.sites.pi &
  done

#把Windows pi合并成染色体

for i in `cat D.name`
do
cat D/1M_window/ *_${i} >> D/by_chr/${i}
done

#ggplot2画图
#改名字
library(ggplot2)
nameA <- c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic")
nameB <- c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic")
nameD <- c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic","DD_Strangulata", "DD_Meyeri", "DD_Anathera")
nameAB <- c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic")
nameABD <- c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","AABBDD_Synthetic")

#修改 路径和name以及输出文件名
path <- "/data1/home/yafei/Project3/Pi/D/by_chr"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 

length(data)

#A: 4, 13
#B: 2, 11
#D: 11
#AB: 10
#ABD:
for(i in 1:length(data)){
  data[[i]]$Name <- nameD[i]
  if (i<11) {
    data[[i]]$Fill <- "pink"
  } else if (i<20) {
    data[[i]]$Fill <- "lightblue"
  } else {
    data[[i]]$Fill <- "orange"
  }
}
all <- data.frame()
for(i in 1:length(data)){
  all <- rbind(all, data[[i]])
}
head(all)
tail(all)

#绘图：修改横坐标顺序和输出文件名
library(ggplot2)
P2<- ggplot(all, aes(x=Name, y=V4,fill=Fill),na.rm=TRUE) + 
  #geom_violin() + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  scale_fill_manual(values = c("#6000b4", "#1f97ff", "#ffce6e"))+
  #Alineage
  #scale_x_discrete(limits= c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Blineage
  #scale_x_discrete(limits= c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Dlineage
  scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera"))+
  #ABlineage
  #scale_x_discrete(limits= c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #ABDlineage
  #scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  theme_bw()+ #背景变为白色
  #B:0.01 D:0.004 AB:0.01
  ylab("Pi")+xlab("Lineages")+scale_y_continuous(limits=c(0,0.004)) + 
  theme(axis.text.x=element_text(angle=80,hjust = 1,colour="black",family="Times",size=200),
        axis.text.y=element_text(family="Times",size=300,face="plain"),
        axis.title.y=element_text(family="Times",size = 300,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())

pdf(file = "D_Pi.pdf",width=200,height=100) #结果保存
print(P2)
dev.off()


#计算平均Pi值
#工作目录：yafei@204:/data1/home/yafei/Project3/Pi/AB/by_chr
for i in `ls`; do awk '{print $4*1000000000}' $i|awk '{sum+=$1} END {printf("%.6f\n",sum/4942000000000)}'; done
#计算标准误.R
path <- "/data1/home/yafei/Project3/Pi/D/by_chr"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
all <- vector()
std <- function(x) sd(x,na.rm=T)/sqrt(length(x))
for(i in 1:length(data)){
  all[i] <-  std(data[[i]][,4])
}
