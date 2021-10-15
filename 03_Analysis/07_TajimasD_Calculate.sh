#-------------------------------------------------------------计算 Tajima's D .sh--------------------------------------------------------------------------------------
#调色板：display.brewer.pal(n = 8, name = 'Dark2')
#分lineage计算Tajima's D: /data2/yafei/Project3/FST_group/lineage/
#工作目录：203:/data2/yafei/Project3/TajimasD/VmapData
#Alineage
for i in {1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --TajimaD 1000000 --out /data2/yafei/Project3/TajimasD/VmapData/A/"pop"$i"_1M.group" &
  done

#Blineage
for i in {4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --TajimaD 1000000 --out /data2/yafei/Project3/TajimasD/VmapData/B/"pop"$i"_1M.group" &
  done

#Dlineage
for i in {14,15,16,17,18,19,20,21,22,23,24,25,26,27}
do
vcftools --gzvcf /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz --keep /data2/yafei/Project3/FST_group/subgroups/"pop"$i".txt" --TajimaD 1000000 --out /data2/yafei/Project3/TajimasD/VmapData/D/"pop"$i"_1M.group" &
  done

#ABlineage

for i in `cat ABfiles`
do
cat /data2/yafei/Project3/TajimasD/VmapData/A/$i  /data2/yafei/Project3/TajimasD/VmapData/B/$i | sed '/CHROM/d' | sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > /data2/yafei/Project3/TajimasD/VmapData/AB/$i
done

#ABDlineage
for i in `cat ABDfiles`
do
cat /data2/yafei/Project3/TajimasD/VmapData/A/$i  /data2/yafei/Project3/TajimasD/VmapData/B/$i  /data2/yafei/Project3/TajimasD/VmapData/D/$i| sed '/CHROM/d' | sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > /data2/yafei/Project3/TajimasD/VmapData/ABD/$i
done

#重新计算新疆的TajimaD
for chr in {1..42}
do
nohup vcftools --gzvcf /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/chr0${chr}.vcf --keep /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/pop19.txt --TajimaD 1000000 --out /data2/xuebo/Projects/Speciation/hapScan/scan_AABBDD/vcf/xinjiang_Tajima/${chr}_pop19_1M.group &
  done

nohup cat 1_pop19_1M.group.Tajima.D 2_pop19_1M.group.Tajima.D 7_pop19_1M.group.Tajima.D 8_pop19_1M.group.Tajima.D 13_pop19_1M.group.Tajima.D 14_pop19_1M.group.Tajima.D 19_pop19_1M.group.Tajima.D 20_pop19_1M.group.Tajima.D 25_pop19_1M.group.Tajima.D 26_pop19_1M.group.Tajima.D 31_pop19_1M.group.Tajima.D 32_pop19_1M.group.Tajima.D 37_pop19_1M.group.Tajima.D 38_pop19_1M.group.Tajima.D | sed '/CHROM/d' |sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > A/pop19_1M.group &
  nohup cat 3_pop19_1M.group.Tajima.D 4_pop19_1M.group.Tajima.D 9_pop19_1M.group.Tajima.D 10_pop19_1M.group.Tajima.D 15_pop19_1M.group.Tajima.D 16_pop19_1M.group.Tajima.D 21_pop19_1M.group.Tajima.D 22_pop19_1M.group.Tajima.D 27_pop19_1M.group.Tajima.D 28_pop19_1M.group.Tajima.D 33_pop19_1M.group.Tajima.D 34_pop19_1M.group.Tajima.D 39_pop19_1M.group.Tajima.D 40_pop19_1M.group.Tajima.D | sed '/CHROM/d' |sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > B/pop19_1M.group &
  nohup cat 5_pop19_1M.group.Tajima.D 6_pop19_1M.group.Tajima.D 11_pop19_1M.group.Tajima.D 12_pop19_1M.group.Tajima.D 17_pop19_1M.group.Tajima.D 18_pop19_1M.group.Tajima.D 23_pop19_1M.group.Tajima.D 24_pop19_1M.group.Tajima.D 29_pop19_1M.group.Tajima.D 30_pop19_1M.group.Tajima.D 35_pop19_1M.group.Tajima.D 36_pop19_1M.group.Tajima.D 41_pop19_1M.group.Tajima.D 42_pop19_1M.group.Tajima.D | sed '/CHROM/d' |sed '1i CHROM\tBIN_START\tN_SNPS\tTajimaD'  > D/pop19_1M.group &
  #转移到yafei203对应位置
  
  #改名字
  library(ggplot2)
nameA <- c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar")
nameB <- c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar")
nameD <- c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera")
nameAB <- c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar")
nameABD <- c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera")

#修改 路径和name以及输出文件名
path <- "/data2/yafei/Project3/TajimasD/VmapData/ABD"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 
length(data)

#A: 4, 13
#B: 2, 11
#D: 11
#AB: 10
#ABD:
for(i in 1:length(data)){
  data[[i]]$Name <- nameABD[i]
  if (i<20) {
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
#绘图：修改输出文件名
library(ggplot2)
P2<- ggplot(all, aes(x=Name, y=TajimaD,fill=Fill),na.rm=TRUE) + 
  #geom_violin() + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5, alpha = 0.8)+ #绘制箱线图
  scale_fill_manual(values = c("#6000b4", "#1f97ff", "#ffce6e"))+
  #Alineage
  #scale_x_discrete(limits= c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Blineage
  #scale_x_discrete(limits= c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Dlineage
  #scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera"))+
  #ABlineage
  #scale_x_discrete(limits= c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #ABDlineage
  scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  theme_bw()+ #背景变为白色
  ylab("TajimaD")+xlab("Lineages")+
  theme(axis.text.x=element_text(angle=80,hjust = 1,colour="black",family="Times",size=200),
        axis.text.y=element_text(family="Times",size=300,face="plain"),
        axis.title.y=element_text(family="Times",size = 300,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")

pdf(file = "ABD_TajimaD.pdf",width=200,height=100) #结果保存
print(P2)
dev.off()


#计算平均TajimaD值以及标准误
#工作目录：yafei@203:/data2/yafei/Project3/TajimasD/VmapData/A
path <- "/data2/yafei/Project3/TajimasD/VmapData/D"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 
SE <- vector()
Mean <- vector()
std <- function(x) sd(x,na.rm=T)/sqrt(length(x))
for(i in 1:length(data)){
  Mean[i] <- mean(data[[i]][,4],na.rm=T)
  SE[i] <-  std(data[[i]][,4])
}
