#----------------------------------------------------------计算Vmap和WW的AABBDD的分离位点并进行统计画图 .sh & .R--------------------------------------------------------------------------------------
#203: /data2/yafei/Project3/V+W/E2_maf/Vmap
#WW
#AABBDD
for chr in {1..42}
do
vcftools --vcf /data2/yafei/Project3/V+W/E2_all_ChangeName/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/WW_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W/E2_maf/Vmap/chr${chr}.WWmafAABBDD.vcf.gz 
done

vcf-concat *WWmafAABBDD.vcf.gz  > WWmafAABBDD.vcf 
bgzip -c WWmafAABBDD.vcf > WWmafAABBDD.vcf.gz

#Vmap
#AABBDD
for chr in {1..42}
do
vcftools --vcf /data2/yafei/Project3/V+W/E2_all_ChangeName/chr${chr}.all.vcf --keep /data2/yafei/Project3/group/AABBDD_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > /data2/yafei/Project3/V+W/E2_maf/Vmap/chr${chr}.VmapmafAABBDD.vcf.gz 
done
vcf-concat *VmapmafAABBDD.vcf.gz  > VmapmafAABBDD.vcf 
bgzip -c VmapmafAABBDD.vcf > VmapmafAABBDD.vcf.gz

#统计
nohup java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W/E2_maf/Vmap/VmapmafAABBDD.vcf.gz --out /data2/yafei/Project3/V+W/E2_maf/Vmap/Vmap_siteQCfile.txt --out2 /data2/yafei/Project3/V+W/E2_maf/Vmap/Vmap_taxaQCfile.txt > nohupVmap 2>& 1 &
  nohup java -jar /data1/home/yafei/C36_checkQuality.jar --file /data2/yafei/Project3/V+W/E2_maf/Vmap/WWmafAABBDD.vcf.gz --out /data2/yafei/Project3/V+W/E2_maf/Vmap/WW_siteQCfile.txt --out2 /data2/yafei/Project3/V+W/E2_maf/Vmap/WW_taxaQCfile.txt > nohupWW 2>& 1 &
  
  #画图
  library(ggplot2)
#site frequency
WW <- read.table("/data2/yafei/Project3/V+W/E2_maf/Vmap/WW_siteQCfile.txt",header=T,stringsAsFactors=F)
Vmap <- read.table("/data2/yafei/Project3/V+W/E2_maf/Vmap/Vmap_siteQCfile.txt",header=T,stringsAsFactors=F)

WW$Lineage <- "WW"
Vmap$Lineage <- "Vmap"

site <- rbind(WW,Vmap)
pdf("Site_Maf.pdf",width=20,height=8)
ggplot(site, aes(x=Maf))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
pdf("Site_MissRate.pdf",width=20,height=8)
ggplot(site, aes(x=MissingRate))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#indi site
Vmap$ID <- paste(Vmap$Chr, Vmap$Pos,sep="-")
WW$ID <- paste(WW$Chr, WW$Pos,sep="-")
WW_a <- WW[which(WW$ID %in% Vmap$ID),]
Vmap_a <- Vmap[which(Vmap$ID %in% WW$ID),]
WW_a$Index <- c(1:60969)
Vmap_a$Index <- c(1:60969)
pdf("WW_Maf.pdf",width=20,height=8)
ggplot(data = WW_a, mapping = aes(x = Index, y = Maf)) + geom_bar(stat = 'identity')+theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
pdf("Vmap_Maf.pdf",width=20,height=8)
ggplot(data = Vmap_a, mapping = aes(x = Index, y = Maf)) + geom_bar(stat = 'identity')+theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

WW_a$Lineage <- "WW"
Vmap_a$Lineage <- "Vmap"
all <- rbind(WW_a,Vmap_a)
pdf("All_Maf.pdf",width=20,height=8)
ggplot(data = all, mapping = aes(x = Index, y = Maf)) + geom_bar(stat = 'identity')+facet_grid(Lineage ~ .)+theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#taxa
WW <- read.table("/data2/yafei/Project3/V+W/E2_maf/Vmap/WW_taxaQCfile.txt",header=T,stringsAsFactors=F)
Vmap <- read.table("/data2/yafei/Project3/V+W/E2_maf/Vmap/Vmap_taxaQCfile.txt",header=T,stringsAsFactors=F)
WW$Lineage <- "WW"
Vmap$Lineage <- "Vmap"

taxa <- rbind(WW,Vmap)
pdf("Mind_MissRate.pdf",width=20,height=8)
ggplot(taxa, aes(x=MissRate))+ geom_histogram() + labs(x="Individual Missing Rate", y="Count")+facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
