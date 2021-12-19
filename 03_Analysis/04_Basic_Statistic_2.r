#---------------------------------------------------plot A/B/Dlineage Taxa/Site MissRate/Maf.R---------------------------------------------------------------------------------------------------
library(ggplot2)
#ABDlineage site
A <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Alineage_siteQCfile.txt",header=T,stringsAsFactors=F)
B <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Blineage_siteQCfile.txt",header=T,stringsAsFactors=F)
D <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Dlineage_siteQCfile.txt",header=T,stringsAsFactors=F)
A$Lineage <- "A-lineage"
B$Lineage <- "B-lineage"
D$Lineage <- "D-lineage"
#colnames(A) <- c("Lineage","A","B","MissingRate","Maf")
#colnames(B) <- c("Lineage","A","B","MissingRate","Maf")
#colnames(D) <- c("Lineage","A","B","MissingRate","Maf")
site <- rbind(A,B,D)
pdf("LineageSite_Maf.pdf",width=20,height=8)
ggplot(site, aes(x=Maf))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
pdf("LineageSite_MissRate.pdf",width=20,height=8)
ggplot(site, aes(x=MissingRate))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#ABDlineage taxa
A <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Alineage_taxaQCfile.txt",header=T,stringsAsFactors=F)
B <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Blineage_taxaQCfile.txt",header=T,stringsAsFactors=F)
D <- read.table("/data2/yafei/Project3/V+W_E2/E2_all/lineage/Dlineage_taxaQCfile.txt",header=T,stringsAsFactors=F)
A$Taxa <- "A-lineage"
B$Taxa <- "B-lineage"
D$Taxa <- "D-lineage"
taxa <- rbind(A,B,D)
colnames(taxa)[1] <- "Lineage"
pdf("LineageMind_MissRate.pdf",width=20,height=8)
ggplot(taxa, aes(x=MissRate))+ geom_histogram() + labs(x="Individual Missing Rate", y="Count")+facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()

#-----------------------------------------------plot AA/BB/DD/AABB/AABBDD/Landrace Taxa/Site MissRate/Maf.R------------------------------------------------------------------------------------------------------------
#AA
nohup cat chr1_siteQCfile.txt chr2_siteQCfile.txt chr7_siteQCfile.txt chr8_siteQCfile.txt chr13_siteQCfile.txt  chr14_siteQCfile.txt chr19_siteQCfile.txt  chr20_siteQCfile.txt  chr25_siteQCfile.txt chr26_siteQCfile.txt  chr31_siteQCfile.txt  chr32_siteQCfile.txt chr37_siteQCfile.txt chr38_siteQCfile.txt > AA_siteQCfile.txt &
  #AABB
  nohup cat chr10_siteQCfile.txt  chr16_siteQCfile.txt  chr21_siteQCfile.txt  chr27_siteQCfile.txt  chr32_siteQCfile.txt  chr38_siteQCfile.txt  chr4_siteQCfile.txt chr13_siteQCfile.txt  chr19_siteQCfile.txt  chr22_siteQCfile.txt  chr28_siteQCfile.txt  chr33_siteQCfile.txt  chr39_siteQCfile.txt  chr7_siteQCfile.txt chr14_siteQCfile.txt  chr1_siteQCfile.txt   chr25_siteQCfile.txt  chr2_siteQCfile.txt   chr34_siteQCfile.txt  chr3_siteQCfile.txt   chr8_siteQCfile.txt chr15_siteQCfile.txt  chr20_siteQCfile.txt  chr26_siteQCfile.txt  chr31_siteQCfile.txt  chr37_siteQCfile.txt  chr40_siteQCfile.txt  chr9_siteQCfile.txt > AABB_siteQCfile.txt &
  #AABBDD
  nohup cat chr10_siteQCfile.txt  chr16_siteQCfile.txt  chr21_siteQCfile.txt  chr27_siteQCfile.txt  chr32_siteQCfile.txt  chr38_siteQCfile.txt  chr4_siteQCfile.txt chr11_siteQCfile.txt  chr17_siteQCfile.txt  chr22_siteQCfile.txt  chr28_siteQCfile.txt  chr33_siteQCfile.txt  chr39_siteQCfile.txt  chr5_siteQCfile.txt chr12_siteQCfile.txt  chr18_siteQCfile.txt  chr23_siteQCfile.txt  chr29_siteQCfile.txt  chr34_siteQCfile.txt  chr3_siteQCfile.txt   chr6_siteQCfile.txt chr13_siteQCfile.txt  chr19_siteQCfile.txt  chr24_siteQCfile.txt  chr2_siteQCfile.txt   chr35_siteQCfile.txt  chr40_siteQCfile.txt  chr7_siteQCfile.txt chr14_siteQCfile.txt  chr1_siteQCfile.txt   chr25_siteQCfile.txt  chr30_siteQCfile.txt  chr36_siteQCfile.txt  chr41_siteQCfile.txt  chr8_siteQCfile.txt chr15_siteQCfile.txt  chr20_siteQCfile.txt  chr26_siteQCfile.txt  chr31_siteQCfile.txt  chr37_siteQCfile.txt  chr42_siteQCfile.txt  chr9_siteQCfile.txt > AABBDD_siteQCfile.txt &
  #DD
  nohup cat chr11_siteQCfile.txt  chr17_siteQCfile.txt  chr23_siteQCfile.txt  chr29_siteQCfile.txt  chr35_siteQCfile.txt  chr41_siteQCfile.txt  chr5_siteQCfile.txt chr12_siteQCfile.txt  chr18_siteQCfile.txt  chr24_siteQCfile.txt  chr30_siteQCfile.txt  chr36_siteQCfile.txt  chr42_siteQCfile.txt  chr6_siteQCfile.txt  > DD_siteQCfile.txt &
  #Landrace
  nohup cat chr10_siteQCfile.txt  chr16_siteQCfile.txt  chr21_siteQCfile.txt  chr27_siteQCfile.txt  chr32_siteQCfile.txt  chr38_siteQCfile.txt  chr4_siteQCfile.txt chr11_siteQCfile.txt  chr17_siteQCfile.txt  chr22_siteQCfile.txt  chr28_siteQCfile.txt  chr33_siteQCfile.txt  chr39_siteQCfile.txt  chr5_siteQCfile.txt chr12_siteQCfile.txt  chr18_siteQCfile.txt  chr23_siteQCfile.txt  chr29_siteQCfile.txt  chr34_siteQCfile.txt  chr3_siteQCfile.txt   chr6_siteQCfile.txt chr13_siteQCfile.txt  chr19_siteQCfile.txt  chr24_siteQCfile.txt  chr2_siteQCfile.txt   chr35_siteQCfile.txt  chr40_siteQCfile.txt  chr7_siteQCfile.txt chr14_siteQCfile.txt  chr1_siteQCfile.txt   chr25_siteQCfile.txt  chr30_siteQCfile.txt  chr36_siteQCfile.txt  chr41_siteQCfile.txt  chr8_siteQCfile.txt chr15_siteQCfile.txt  chr20_siteQCfile.txt  chr26_siteQCfile.txt  chr31_siteQCfile.txt  chr37_siteQCfile.txt  chr42_siteQCfile.txt  chr9_siteQCfile.txt Landrace_siteQCfile.txt &
  #SS
  nohup cat chr10_siteQCfile.txt  chr16_siteQCfile.txt  chr22_siteQCfile.txt  chr28_siteQCfile.txt  chr34_siteQCfile.txt  chr3_siteQCfile.txt   chr4_siteQCfile.txt chr15_siteQCfile.txt  chr21_siteQCfile.txt  chr27_siteQCfile.txt  chr33_siteQCfile.txt  chr39_siteQCfile.txt  chr40_siteQCfile.txt  chr9_siteQCfile.txt > SS_siteQCfile.txt &
  
  #site_miss & site_maf
AA <- read.table("/data1/home/yafei/Project3/Maf/mafAA/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
AA_maf <- AA[,c(1,5)]
AA_miss <- AA[,c(1,4)]
AA_maf$V1 <- "AA"
AA_miss$V1 <- "AA"

SS <- read.table("/data1/home/yafei/Project3/Maf/mafSS/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
SS_maf <- SS[,c(1,5)]
SS_miss <- SS[,c(1,4)]
SS_maf$V1 <- "SS"
SS_miss$V1 <- "SS"

DD <- read.table("/data1/home/yafei/Project3/Maf/mafDD/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
DD_maf <- DD[,c(1,5)]
DD_miss <- DD[,c(1,4)]
DD_maf$V1 <- "DD"
DD_miss$V1 <- "DD"

AABB <- read.table("/data1/home/yafei/Project3/Maf/mafAABB/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
AABB_maf <- AABB[,c(1,5)]
AABB_miss <- AABB[,c(1,4)]
AABB_maf$V1 <- "AABB"
AABB_miss$V1 <- "AABB"

AABBDD <- read.table("/data1/home/yafei/Project3/Maf/mafAABBDD/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
AABBDD_maf <- AABBDD[,c(1,5)]
AABBDD_miss <- AABBDD[,c(1,4)]
AABBDD_maf$V1 <- "AABBDD"
AABBDD_miss$V1 <- "AABBDD"

Landrace <- read.table("/data1/home/yafei/Project3/Maf/mafLand/100k_siteQCfile.txt",header=F,stringsAsFactors=F)
Landrace_maf <- Landrace[,c(1,5)]
Landrace_miss <- Landrace[,c(1,4)]
Landrace_maf$V1 <- "Landrace"
Landrace_miss$V1 <- "Landrace"

maf <- rbind(AA_maf,SS_maf,DD_maf,AABB_maf,AABBDD_maf)
colnames(maf) <- c("Lineage","Maf")
maf$Maf <- as.numeric(maf$Maf)

miss <- rbind(AA_miss,SS_miss,DD_miss,AABB_miss,AABBDD_miss)
colnames(miss) <- c("Lineage","MissRate")
miss$MissRate <- as.numeric(miss$MissRate)

library(ggplot2)
pdf("site_Maf.pdf",width=20,height=8)
ggplot(maf, aes(x=Maf))+ geom_histogram() + facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=25),panel.border = element_blank())
dev.off()

pdf("site_MissRate.pdf",width=20,height=8)
ggplot(miss, aes(x=MissRate))+ geom_histogram() + labs(x="Site Missing Rate", y="count")+facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=25),panel.border = element_blank())
dev.off()

**mind_miss
AA <- read.table("/data1/home/yafei/Project3/Maf/mafAA/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
SS <- read.table("/data1/home/yafei/Project3/Maf/mafSS/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
DD <- read.table("/data1/home/yafei/Project3/Maf/mafDD/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
AABB <- read.table("/data1/home/yafei/Project3/Maf/mafAABB/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
AABBDD <- read.table("/data1/home/yafei/Project3/Maf/mafAABBDD/all_taxaQCfile.txt",header=F,stringsAsFactors=F)
Landrace <- read.table("/data1/home/yafei/Project3/Maf/mafLand/all_taxaQCfile.txt",header=F,stringsAsFactors=F)

AA$V1 <- "AA"
SS$V1 <- "SS"
DD$V1 <- "DD"
AABB$V1 <- "AABB"
AABBDD$V1 <- "AABBDD"
Landrace$V1 <- "Landrace"

mind<- rbind(AA,SS,DD,AABB,AABBDD,Landrace)
colnames(mind) <- c("Lineage","Het","MissingRate")
mind$MissingRate <- as.numeric(mind$MissingRate)

**应该画密度
pdf("Mind_MissRate.pdf",width=20,height=8)
ggplot(mind, aes(x=MissingRate,y=..density..))+ geom_histogram() + labs(x="Individual Missing Rate", y="Frequency")+facet_grid(. ~ Lineage) +theme(panel.background = element_blank(),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 20),strip.text.x = element_text(size=15),panel.border = element_blank())
dev.off()
