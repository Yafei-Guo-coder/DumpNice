#-------------------------------------------------------------------CS IBS--------------------------------------------------------------------------------------
#提取CS的pos和posAllele文件做hapscanner。工作目录：203:/data1/home/yafei/Project3/Vmap1.1
for chr in {1..42}
do
zcat chr${chr}.all.vcf.gz |awk '{if(NF>100) print $1"\t"$2}' |sed '1d' > CS_pos/chr${chr}_pos.txt &
  zcat chr${chr}.all.vcf.gz |awk '{if(NF>100) print $1"\t"$2"\t"$4"\t"$5}' |sed '1c Chr\tPos\tRef\tAlt' > CS_posAllele/chr${chr}_posAllele.txt  &
  done

#生成参数文件 工作目录：204:/data2/xuebo/Projects/Speciation/More_accessions/hapScan/
nohup java -jar Yafei_Guo.jar &
  
  #进行CS的hapscanner扫描
  nohup bash bashCS.sh &
  #生成的结果文件 工作目录：203:/data2/yafei/Project3/CS_Vmap/CS/
  #合并lineage文件
  nohup vcf-concat chr01.vcf chr02.vcf chr07.vcf chr08.vcf chr013.vcf chr014.vcf chr019.vcf chr020.vcf chr025.vcf chr026.vcf chr031.vcf chr032.vcf chr037.vcf chr038.vcf | bgzip -c > AlineageCS.vcf.gz &
  nohup vcf-concat chr03.vcf chr04.vcf chr09.vcf chr010.vcf chr015.vcf chr016.vcf chr021.vcf chr022.vcf chr027.vcf chr028.vcf chr033.vcf chr034.vcf chr039.vcf chr040.vcf | bgzip -c > BlineageCS.vcf.gz &
  nohup vcf-concat chr05.vcf chr06.vcf chr011.vcf chr012.vcf chr017.vcf chr018.vcf chr023.vcf chr024.vcf chr029.vcf chr030.vcf chr035.vcf chr036.vcf chr041.vcf chr042.vcf | bgzip -c > DlineageCs.vcf.gz &
  nohup bcftools merge /data1/home/yafei/Project3/Vmap1.1/Alineage.vcf.gz /data2/yafei/Project3/CS_Vmap/CS/AlineageCS.vcf.gz -o Alineage_withCS.vcf &
  nohup bcftools merge /data1/home/yafei/Project3/Vmap1.1/Blineage.vcf.gz /data2/yafei/Project3/CS_Vmap/CS/BlineageCS.vcf.gz -o Blineage_withCS.vcf &
  nohup bcftools merge /data1/home/yafei/Project3/Vmap1.1/Dlineage.vcf.gz /data2/yafei/Project3/CS_Vmap/CS/DlineageCS.vcf.gz -o Dlineage_withCS.vcf &
  
  #提取barely的vcf
  这个文件夹是分染色体的: /data2/xuebo/Projects/Speciation/tree/withBarley_segregate
这个是合起来的: /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage

#计算IBS:barely
nohup java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage/Alineage_withBarley.vcf.gz --out Alineage_withBarley.all.ibs.txt > logA.txt &
  nohup java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage/Blineage_withBarley.vcf.gz --out Blineage_withBarley.all.ibs.txt > logB.txt &
  nohup java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 /data2/xuebo/Projects/Speciation/tree/withBarley_segregate/lineage/Dlineage_withBarley.vcf.gz --out Dlineage_withBarley.all.ibs.txt > logD.txt &
  
  #计算IBS:CS
  nohup java -jar -Xmx100g /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 Alineage_withCS.vcf.gz --out Alineage_withCS.all.ibs.txt > logA.txt &
  nohup java -jar /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 Blineage_withCS.vcf.gz --out Blineage_withCS.all.ibs.txt > logB.txt &
  nohup java -jar -Xmx100g /data2/xuebo/Projects/Speciation/javaCode/C41_getIBS_distance2.jar --file1 Dlineage_withCS.vcf.gz --out Dlineage_withCS.all.ibs.txt > logD.txt &
  
  #IBS画图 分不同的lineage
  #工作目录：yafei@203:/data2/yafei/Project3/CS_Vmap
  #读取group名的文件
  path <- "/data2/yafei/Project3/group/subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)

#修改 输入文件名
A <- read.table("Alineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
A$Lineage <- NA
A$Fill <- NA
nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_Durum","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_Khorasan_wheat","AABBDD_Landrace","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_Polish_wheat","AABB_Rivet_wheat","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
#A:487  485-487
#B:406  404-406
#D:306  304-306
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
sub <- A[,485:487]

library(ggplot2)
P2<- ggplot(sub, aes(x=Lineage, y=CS, fill=Fill),na.rm=TRUE) + 
  #geom_violin() + 
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  scale_fill_manual(values = c("#6000b4", "#1f97ff", "#ffce6e"))+
  #Alineage
  #scale_x_discrete(limits= c("A_Wild_einkorn","A_Domesticated_einkorn","A_Urartu","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Blineage
  scale_x_discrete(limits= c("B_Speltoides","AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #Dlineage
  #scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar","DD_Strangulata", "DD_Meyeri", "DD_Anathera"))+
  #ABlineage
  #scale_x_discrete(limits= c("AABB_Wild_emmer","AABB_Domesticated_emmer","AABB_Georgian_wheat","AABB_Ispahanicum","AABB_Rivet_wheat","AABB_Polish_wheat","AABB_Persian_wheat","AABB_Khorasan_wheat","AABB_Durum","AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  #ABDlineage
  #scale_x_discrete(limits= c("AABBDD_Spelt","AABBDD_Macha","AABBDD_Club_wheat","AABBDD_Indian_dwarf_wheat","AABBDD_Yunan_wheat","AABBDD_Xinjiang_wheat","AABBDD_Tibetan_semi_wild","AABBDD_Vavilovii","AABBDD_Landrace","AABBDD_Cultivar"))+
  theme_bw()+ #背景变为白色
  #A:0.4 B:0.4 D:0.7
  ylab("IBS_CS")+xlab("Lineages")+scale_y_continuous(limits=c(0,0.4)) + 
  theme(axis.text.x=element_text(angle=80, hjust = 1, colour="black", family="Times", size=200),
        axis.text.y=element_text(family="Times",size=300,face="plain"),
        axis.title.y=element_text(family="Times",size = 300,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())

pdf(file = "B_IBS_CS.pdf",width=200,height=100) #结果保存
print(P2)
dev.off()

#画热图文件的产生
library("corrplot")
library("reshape")
path <- "/data2/yafei/Project3/group/subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)
D <- read.table("Blineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
rownames(D) <- D$Dxy
D<- D[,-1]
D<- as.matrix(D)
all <- as.data.frame(melt(D))
#colnames(all) <- c("name1","name2","value")
all$X1 <- as.character(all$X1)
all$X2 <- as.character(all$X2)
nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_Durum","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_Khorasan_wheat","AABBDD_Landrace","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_Polish_wheat","AABB_Rivet_wheat","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
for(i in 1:length(data)){
  all[which(all$X1 %in% data[[i]][,1]),1] <- nameABD[i]
}
for(i in 1:length(data)){
  all[which(all$X2 %in% data[[i]][,1]),2] <- nameABD[i]
}
colnames(all) <- c("id1","id2","value")
all$variable <- "var"
sub <- all[-which(all$value == 0),]

casted <- cast(sub,id1+id2~variable,mean)
ibs <- cast(casted, id1~id2)
ibs[is.na(ibs)] <- 0
names <- ibs$id1
write.table(ibs,"AB_ibs.txt",quote=F) 

#热图  
library("corrplot")
ibs <- read.table("D_ibs.txt",header=T,stringsAsFactors=F)
names <- ibs$id1
ibs <- ibs[,-1]
ibs <- as.matrix(ibs)
rownames(ibs) <- names
colnames(ibs) <- names

#A:0.5 B:0.5 D:0.7
pdf("IBS_D_CS_heat.pdf",width=100,height=100)
#corrplot(data,method = "color",col = col3(10),tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=2,tl.cex=2,cl.cex=2,cl.lim = c(0, 1))
corrplot(ibs,method = "color",tl.col="black",tl.srt = 45, addrect=4,addCoef.col = "grey", type = "lower",number.cex=6,number.digits=3,tl.cex=10,cl.cex=12,cl.lim = c(0, 0.7))
dev.off()

