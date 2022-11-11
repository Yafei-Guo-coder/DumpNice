#yafei204: /data4/home/yafei/plink_VCF/
#过滤: maf0.01, R2: 0.9 (/data4/home/yafei/plink_VCF/maf001r209)
for i in {001..042}
do
plink --bfile /data4/home/yafei/plink_VCF/chr${i} --maf 0.01 --make-bed --out chr${i}.filter --allow-extra-chr --autosome-num 42 &
done
for i in {001..042}
do
plink --bfile chr${i}.filter --indep-pairwise 50 10 0.9 --out chr${i}.filterLD --allow-extra-chr --autosome-num 42 &
done
for i in {001..042}
do
plink --bfile chr${i}.filter --make-bed --extract chr${i}.filterLD.prune.in --out chr${i}.filter.prune.in --autosome-num 42 &
done
#提取群体样本(/data2/yafei/004_Vmap3/Group/type_8/plink_group/)
#(domemmer.txt,landrace.txt,other_tetraploid.txt,wildemmer.txt,cultivar.txt,freethresh.txt,other_hexaploid.txt,strangulata.txt)
for i in {027..042}
do
plink --bfile /data4/home/yafei/plink_VCF/maf001r209/bfile/chr${i}.filter.prune.in --keep /data2/yafei/004_Vmap3/Group/type_8/plink_group/landrace.txt --recode vcf-iid --out group/landrace.chr${i} --autosome-num 42 &
plink --bfile /data4/home/yafei/plink_VCF/maf001r209/bfile/chr${i}.filter.prune.in --keep /data2/yafei/004_Vmap3/Group/type_8/plink_group/cultivar.txt --recode vcf-iid --out group/cultivar.chr${i} --autosome-num 42 &
done
for i in {001,002,003,004,007,008,009,010,013,014,015,016,019,020,021,022,025,026,027,028,031,032,033,034,037,038,039,040}
do
plink --bfile /data4/home/yafei/plink_VCF/maf001r209/bfile/chr${i}.filter.prune.in --keep /data2/yafei/004_Vmap3/Group/type_8/plink_group/domemmer.txt --recode vcf-iid --out group/domemmer.chr${i} --autosome-num 42 &
plink --bfile /data4/home/yafei/plink_VCF/maf001r209/bfile/chr${i}.filter.prune.in --keep /data2/yafei/004_Vmap3/Group/type_8/plink_group/wildemmer.txt --recode vcf-iid --out group/wildemmer.chr${i} --autosome-num 42 &
plink --bfile /data4/home/yafei/plink_VCF/maf001r209/bfile/chr${i}.filter.prune.in --keep /data2/yafei/004_Vmap3/Group/type_8/plink_group/freethresh.txt --recode vcf-iid --out group/freethresh.chr${i} --autosome-num 42 &
done
for i in {005,006,011,012,017,018,023,024,029,030,035,036,041,042}
do
plink --bfile /data4/home/yafei/plink_VCF/maf001r209/bfile/chr${i}.filter.prune.in --keep /data2/yafei/004_Vmap3/Group/type_8/plink_group/strangulata.txt --recode vcf-iid --out group/strangulata.chr${i} --autosome-num 42 &
done
vcf-concat domemmer.chr001.vcf.gz domemmer.chr002.vcf.gz domemmer.chr003.vcf.gz domemmer.chr004.vcf.gz domemmer.chr007.vcf.gz domemmer.chr008.vcf.gz domemmer.chr009.vcf.gz domemmer.chr010.vcf.gz domemmer.chr013.vcf.gz domemmer.chr014.vcf.gz domemmer.chr015.vcf.gz domemmer.chr016.vcf.gz domemmer.chr019.vcf.gz domemmer.chr020.vcf.gz domemmer.chr021.vcf.gz domemmer.chr022.vcf.gz domemmer.chr025.vcf.gz domemmer.chr026.vcf.gz domemmer.chr027.vcf.gz domemmer.chr028.vcf.gz domemmer.chr031.vcf.gz domemmer.chr032.vcf.gz domemmer.chr033.vcf.gz domemmer.chr034.vcf.gz domemmer.chr037.vcf.gz domemmer.chr038.vcf.gz domemmer.chr039.vcf.gz domemmer.chr040.vcf.gz| bgzip -c > lineage/dommer.vcf.gz
vcf-concat landrace.chr001.vcf.gz landrace.chr002.vcf.gz landrace.chr003.vcf.gz landrace.chr004.vcf.gz landrace.chr005.vcf.gz landrace.chr006.vcf.gz landrace.chr007.vcf.gz landrace.chr008.vcf.gz landrace.chr009.vcf.gz landrace.chr010.vcf.gz landrace.chr011.vcf.gz landrace.chr012.vcf.gz landrace.chr013.vcf.gz landrace.chr014.vcf.gz landrace.chr015.vcf.gz landrace.chr016.vcf.gz landrace.chr017.vcf.gz landrace.chr018.vcf.gz landrace.chr019.vcf.gz landrace.chr020.vcf.gz landrace.chr021.vcf.gz landrace.chr022.vcf.gz landrace.chr023.vcf.gz landrace.chr024.vcf.gz landrace.chr025.vcf.gz landrace.chr026.vcf.gz landrace.chr027.vcf.gz landrace.chr028.vcf.gz landrace.chr029.vcf.gz landrace.chr030.vcf.gz landrace.chr031.vcf.gz landrace.chr032.vcf.gz landrace.chr033.vcf.gz landrace.chr034.vcf.gz landrace.chr035.vcf.gz landrace.chr036.vcf.gz landrace.chr037.vcf.gz landrace.chr038.vcf.gz landrace.chr039.vcf.gz landrace.chr040.vcf.gz landrace.chr041.vcf.gz landrace.chr042.vcf.gz | bgzip -c > landrace.vcf.gz
vcf-concat cultivar.chr001.vcf.gz cultivar.chr002.vcf.gz cultivar.chr003.vcf.gz cultivar.chr004.vcf.gz cultivar.chr005.vcf.gz cultivar.chr006.vcf.gz cultivar.chr007.vcf.gz cultivar.chr008.vcf.gz cultivar.chr009.vcf.gz cultivar.chr010.vcf.gz cultivar.chr011.vcf.gz cultivar.chr012.vcf.gz cultivar.chr013.vcf.gz cultivar.chr014.vcf.gz cultivar.chr015.vcf.gz cultivar.chr016.vcf.gz cultivar.chr017.vcf.gz cultivar.chr018.vcf.gz cultivar.chr019.vcf.gz cultivar.chr020.vcf.gz cultivar.chr021.vcf.gz cultivar.chr022.vcf.gz cultivar.chr023.vcf.gz cultivar.chr024.vcf.gz cultivar.chr025.vcf.gz cultivar.chr026.vcf.gz cultivar.chr027.vcf.gz cultivar.chr028.vcf.gz cultivar.chr029.vcf.gz cultivar.chr030.vcf.gz cultivar.chr031.vcf.gz cultivar.chr032.vcf.gz cultivar.chr033.vcf.gz cultivar.chr034.vcf.gz cultivar.chr035.vcf.gz cultivar.chr036.vcf.gz cultivar.chr037.vcf.gz cultivar.chr038.vcf.gz cultivar.chr039.vcf.gz cultivar.chr040.vcf.gz cultivar.chr041.vcf.gz cultivar.chr042.vcf.gz | bgzip -c > cultivar.vcf.gz
vcf-concat freethresh.chr001.vcf.gz freethresh.chr002.vcf.gz freethresh.chr003.vcf.gz freethresh.chr004.vcf.gz freethresh.chr007.vcf.gz freethresh.chr008.vcf.gz freethresh.chr009.vcf.gz freethresh.chr010.vcf.gz freethresh.chr013.vcf.gz freethresh.chr014.vcf.gz freethresh.chr015.vcf.gz freethresh.chr016.vcf.gz freethresh.chr019.vcf.gz freethresh.chr020.vcf.gz freethresh.chr021.vcf.gz freethresh.chr022.vcf.gz freethresh.chr025.vcf.gz freethresh.chr026.vcf.gz freethresh.chr027.vcf.gz freethresh.chr028.vcf.gz freethresh.chr031.vcf.gz freethresh.chr032.vcf.gz freethresh.chr033.vcf.gz freethresh.chr034.vcf.gz freethresh.chr037.vcf.gz freethresh.chr038.vcf.gz freethresh.chr039.vcf.gz freethresh.chr040.vcf.gz | bgzip -c > freethresh.vcf.gz
vcf-concat wildemmer.chr001.vcf.gz wildemmer.chr002.vcf.gz wildemmer.chr003.vcf.gz wildemmer.chr004.vcf.gz wildemmer.chr007.vcf.gz wildemmer.chr008.vcf.gz wildemmer.chr009.vcf.gz wildemmer.chr010.vcf.gz wildemmer.chr013.vcf.gz wildemmer.chr014.vcf.gz wildemmer.chr015.vcf.gz wildemmer.chr016.vcf.gz wildemmer.chr019.vcf.gz wildemmer.chr020.vcf.gz wildemmer.chr021.vcf.gz wildemmer.chr022.vcf.gz wildemmer.chr025.vcf.gz wildemmer.chr026.vcf.gz wildemmer.chr027.vcf.gz wildemmer.chr028.vcf.gz wildemmer.chr031.vcf.gz wildemmer.chr032.vcf.gz wildemmer.chr033.vcf.gz wildemmer.chr034.vcf.gz wildemmer.chr037.vcf.gz wildemmer.chr038.vcf.gz wildemmer.chr039.vcf.gz wildemmer.chr040.vcf.gz | bgzip -c > wildemmer.vcf.gz
vcf-concat strangulata.chr005.vcf.gz strangulata.chr006.vcf.gz strangulata.chr011.vcf.gz strangulata.chr012.vcf.gz strangulata.chr017.vcf.gz strangulata.chr018.vcf.gz strangulata.chr023.vcf.gz strangulata.chr024.vcf.gz strangulata.chr029.vcf.gz strangulata.chr030.vcf.gz strangulata.chr035.vcf.gz strangulata.chr036.vcf.gz strangulata.chr041.vcf.gz strangulata.chr042.vcf.gz | bgzip -c > strangulata.vcf.gz

#sitePI
for i in {"cultivar","landrace","domemmer","freethresh","wildemmer"}
do
vcftools --gzvcf ${i}.A.vcf.gz --site-pi --out ${i}.A &
vcftools --gzvcf ${i}.B.vcf.gz --site-pi --out ${i}.B &
done
for i in {"cultivar","landrace","strangulata"}
do
vcftools --gzvcf ${i}.D.vcf.gz --site-pi --out ${i}.D &
done
#admixture
for i in {"cultivar","landrace","strangulata.D","domemmer","freethresh","wildemmer"}
do
  #zcat ${i}.vcf.gz | awk '{print $1"\t"$2}' | grep -v "#" | shuf -n 500000 | sort -k1,1n -k2,2n > ${i}.pos & 
  #vcftools --gzvcf ${i}.vcf.gz --positions ${i}.pos --recode --recode-INFO-all --out ${i}.05M.vcf &
  plink --vcf ${i}.05M.vcf.recode.vcf --recode 12 --out ${i}.admix --autosome-num 42 &
done   
for i in {"cultivar","landrace","strangulata.D","domemmer","freethresh","wildemmer"}
do
  admixture --cv ${i}.admix.ped 2 >> ${i}.log.txt &
  admixture --cv ${i}.admix.ped 3 >> ${i}.log.txt & 
  admixture --cv ${i}.admix.ped 4 >> ${i}.log.txt &
  admixture --cv ${i}.admix.ped 5 >> ${i}.log.txt &
  admixture --cv ${i}.admix.ped 6 >> ${i}.log.txt &
  admixture --cv ${i}.admix.ped 7 >> ${i}.log.txt &
done
  
#画图
  wildemmer k=7; domemmer k=5; freethresh k=4;
  #统一设置：注释内容及颜色
  library(pheatmap)
  require(reshape)
  require (rworldmap)
  require(rworldxtra)
  library(pophelper)
  library(ggplot2)
  require(gridExtra)
  library(RColorBrewer)
  setwd("/Users/guoyafei/Documents/02_VmapIII/04_Statistics/02_分群/04_admixture/")
  #load Qmatrix files
  sfiles <- list.files(path="/Users/guoyafei/Documents/02_VmapIII/04_Statistics/02_分群/04_admixture/", full.names=T)
  slist <- readQ(files=sfiles)
  tabulateQ(qlist=readQ(sfiles))
  summariseQ(tabulateQ(qlist=readQ(sfiles)))
  
  pdf("test.pdf",width = 7,height = 5)
  p1 <- plotQ(slist[1:3],returnplot=T,exportplot=F,quiet=T,basesize=11,
              sortind="all",showindlab=T,showyaxis=T,showticks=T,sharedindlab=T,clustercol=brewer.pal(7, "Set2"))
  grid.arrange(p1$plot[[3]],p1$plot[[1]],p1$plot[[2]],nrow=3)
  dev.off()
#fst
#/data2/yafei/004_Vmap3/Group/type_8/subGroup/wildemmer_sub7.txt;domemmer_sub5.txt;freethresh_sub4.txt;
#/data4/home/yafei/plink_VCF/maf001r209/group/lineage/wildemmer.vcf.gz;domemmer.vcf.gz;freethresh.vcf.gz;landrace.vcf.gz;cultivar.vcf.gz;
  for i in {"sub1","sub2","sub3","sub4"}
  do  
  vcftools --gzvcf wildemmer.vcf.gz --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/${pop1} --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/${pop2} --fst-window-size 100000 --out wildemmer.${pop1}_${pop2}
  vcftools --gzvcf domemmer.vcf.gz --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/${pop1} --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/${pop2} --fst-window-size 100000 --out domemmer.${pop1}_${pop2}
  vcftools --gzvcf freethresh.vcf.gz --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/${pop1} --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/${pop2} --fst-window-size 100000 --out freethresh.${pop1}_${pop2}
  done
  
cat wildemmer_fst.txt | while read pop1 pop2
do
vcftools --gzvcf wildemmer.vcf.gz --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/wildemmer_${pop1}.txt --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/wildemmer_${pop2}.txt --fst-window-size 100000 --out wildemmer.${pop1}_${pop2}
done
cat domemmer_fst.txt | while read pop1 pop2
do
vcftools --gzvcf domemmer.vcf.gz --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/domemmer_${pop1}.txt --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/domemmer_${pop2}.txt --fst-window-size 100000 --out domemmer.${pop1}_${pop2}
done
cat freethresh_fst.txt | while read pop1 pop2
do
vcftools --gzvcf freethresh.vcf.gz --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/freethresh_${pop1}.txt --weir-fst-pop /data2/yafei/004_Vmap3/Group/type_8/subGroup/freethresh_${pop2}.txt --fst-window-size 100000 --out freethresh.${pop1}_${pop2}
done


setwd("/Users/guoyafei/Documents/02_VmapIII/11_piratio")
data <- read.table("landrace-cultivar.pi.ratio", header=F,stringsAsFactors = F)
ggplot(data, aes(x=V4, y=log(V3), color=V4)) + 
  geom_boxplot() + 
  #facet_grid(V4~.)+
  #scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  xlab("")+ylab("")+theme_bw()+
  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.position = "null",axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
library(tidyverse)
x_order <- data %>%
  group_by(V4) %>%
  summarize(mean_y=mean(V3))%>%
  ungroup()%>%
  arrange(desc(mean_y))%>%
  select(V4);
data$V4<-factor(data$V4,levels=as.character(x_order$V4),ordered = TRUE)

setwd("/Users/guoyafei/Documents/02_VmapIII/12_H12-H2:H1")
landrace <- read.table("landrace.txt",header=F,stringsAsFactors = F)
cultivar <- read.table("cultivar.txt",header=F,stringsAsFactors = F)
landrace$type <- "landrace"
cultivar$type <- "cultivar"

all <- rbind(landrace,cultivar)
ggplot(all, aes(x=V1, y=V2, color=type)) + 
  geom_bar(stat="identity",width=0.5,position='dodge') + 
  #geom_hline(yintercept = 0.059, color = 'gray', size = 0.5) +
  #facet_grid(variable~.)+
  #scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  xlab("")+ylab("")+theme_bw()+
  
  3.376

  theme(plot.title = element_text(color="red", size=10, face="bold.italic"),legend.text = element_text(size=12),legend.title=element_blank(),axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y = element_text(size = 10))




