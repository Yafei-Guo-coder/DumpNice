#shell requirement------
#计算RDA：环境变量和遗传变异
#工作目录：204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225/Lineages
#计算目录：203:yafei:/data1/home/yafei/003_Project3/RDA

#全基因组随机抽取了20k snp,把缺失的填充成major allele.
#1. 提取计算迁徙路径的样本的VCF文件
vcftools --gzvcf ABlineage.E6_Landrace_locate.vcf.gz --max-missing 1 --maf 0.05 --recode --stdout | bgzip -c > AB_noMiss_0.05.vcf.gz
vcftools --gzvcf Dlineage.E6_Landrace_locate.vcf.gz --max-missing 1 --maf 0.05 --recode --stdout | bgzip -c > D_noMiss_0.05.vcf.gz
vcf-concat AB_noMiss_0.05.vcf.gz D_noMiss_0.05.vcf.gz | bgzip -c > noSort_noMiss_0.05.vcf.gz 
zcat noSort_noMiss_0.05.vcf.gz | vcf-sort |bgzip -c > Sorted_noMiss_0.05.vcf.gz

#2. 提取基因区的VCF文件
zcat D_noMiss_0.05.vcf.gz | head -n 40 > header.txt
bedtools intersect -b gene_v1.1_Lulab.gff3 -a sorted_noMiss_0.05.vcf.gz -wa | shuf -n 3000 | sort -k1,1n -k2,2n > Select_gene_noHeader.vcf
cat header.txt Select_gene_noHeader.vcf > gene_3000.vcf
sed '/^##/d' gene_3000.vcf | awk '{$1=null;$2=null;$3=null;$4=null;$5=null;$6=null;$7=null;$8=null;$9=null;print $0'} | sed 's/^[ \t]*//g' > all_noMiss_0.05_3000.txt
sed 's@0/0@0@g' all_noMiss_0.05_3000.txt | sed 's@0/1@1@g' |sed 's@1/1@2@g' |sed 's@1/0@1@g' > All_noMiss_0.05_3000.txt

#分亚基因组提取基因区的VCF文件，并随机抽取3000个SNP(实际运行多抽一些，可以设置30-50k)
bedtools intersect -b gene_v1.1_Lulab.gff3 -a sorted_noMiss_0.05.vcf.gz -wa > sorted_noMiss_gene_noHeader.vcf 
awk '{output="chr"$1".noheader.vcf"; print $0 > $output}' sorted_noMiss_gene_noHeader.vcf
cat chr1.noheader.vcf chr2.noheader.vcf chr7.noheader.vcf chr8.noheader.vcf chr13.noheader.vcf chr14.noheader.vcf chr19.noheader.vcf chr20.noheader.vcf chr25.noheader.vcf chr26.noheader.vcf chr31.noheader.vcf chr32.noheader.vcf chr37.noheader.vcf chr38.noheader.vcf | shuf -n 3000 | sort -k1,1n -k2,2n > Select_gene_noHeader_A.vcf
cat chr3.noheader.vcf chr4.noheader.vcf chr9.noheader.vcf chr10.noheader.vcf chr15.noheader.vcf chr16.noheader.vcf chr21.noheader.vcf chr22.noheader.vcf chr27.noheader.vcf chr28.noheader.vcf chr33.noheader.vcf chr34.noheader.vcf chr39.noheader.vcf chr40.noheader.vcf | shuf -n 3000 | sort -k1,1n -k2,2n > Select_gene_noHeader_B.vcf
cat chr5.noheader.vcf chr6.noheader.vcf chr11.noheader.vcf chr12.noheader.vcf chr17.noheader.vcf chr18.noheader.vcf chr23.noheader.vcf chr24.noheader.vcf chr29.noheader.vcf chr30.noheader.vcf chr35.noheader.vcf chr36.noheader.vcf chr41.noheader.vcf chr42.noheader.vcf | shuf -n 3000 | sort -k1,1n -k2,2n > Select_gene_noHeader_D.vcf
cat header.txt Select_gene_noHeader_A.vcf > gene_3000_A.vcf
cat header.txt Select_gene_noHeader_B.vcf > gene_3000_B.vcf
cat header.txt Select_gene_noHeader_D.vcf > gene_3000_D.vcf
sed '/^##/d' gene_3000_A.vcf | awk '{$1=null;$2=null;$3=null;$4=null;$5=null;$6=null;$7=null;$8=null;$9=null;print $0'} | sed 's/^[ \t]*//g' > all_noMiss_0.05_3000_A.txt
sed 's@0/0@0@g' all_noMiss_0.05_3000_A.txt | sed 's@0/1@1@g' |sed 's@1/1@2@g' |sed 's@1/0@1@g' > All_noMiss_0.05_3000_A.txt
sed '/^##/d' gene_3000_B.vcf | awk '{$1=null;$2=null;$3=null;$4=null;$5=null;$6=null;$7=null;$8=null;$9=null;print $0'} | sed 's/^[ \t]*//g' > all_noMiss_0.05_3000_B.txt
sed 's@0/0@0@g' all_noMiss_0.05_3000_B.txt | sed 's@0/1@1@g' |sed 's@1/1@2@g' |sed 's@1/0@1@g' > All_noMiss_0.05_3000_B.txt
sed '/^##/d' gene_3000_D.vcf | awk '{$1=null;$2=null;$3=null;$4=null;$5=null;$6=null;$7=null;$8=null;$9=null;print $0'} | sed 's/^[ \t]*//g' > all_noMiss_0.05_3000_D.txt
sed 's@0/0@0@g' all_noMiss_0.05_3000_D.txt | sed 's@0/1@1@g' |sed 's@1/1@2@g' |sed 's@1/0@1@g' > All_noMiss_0.05_3000_D.txt
awk '{for (i = 1; i <= NF; ++i) {split($i, array, ":"); print array[1]}}' All_noMiss_0.05_3000_A.txt | xargs -n225 > All_noMiss_0.05_3000_A2.txt
awk '{for (i = 1; i <= NF; ++i) {split($i, array, ":"); print array[1]}}' All_noMiss_0.05_3000_B.txt | xargs -n225 > All_noMiss_0.05_3000_B2.txt
awk '{for (i = 1; i <= NF; ++i) {split($i, array, ":"); print array[1]}}' All_noMiss_0.05_3000_D.txt | xargs -n225 > All_noMiss_0.05_3000_D2.txt

#3. 样本分区：EU, WA, SCA, EA_North, EA_SouthWest
#4. RDA analysis(bar plot)
#input files:
#env_table: data frame, row: sample name, col: environment variables
#genotype_table: data frame, row: sample name, col: snp site
library(vegan)
library(RColorBrewer)
library(ggplot2)
phylum <- read.delim('/Users/guoyafei/Documents/01_VMap1.1/02_Environment/01_RDA_plot/All_noMiss_0.05_2000.txt',  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
row.names(phylum) <- c(1:20000)
phylum <- data.frame(t(phylum))
colnames(phylum) <- c(1:20000)
phylum <- phylum[which(rownames(phylum)!="TW095"),]
env <- read.delim('select_bio2.txt', row.names = 1, header=T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
env_all <- data.frame(env[,1:20])
env_all <- env_all[which(rownames(env_all)!="TW095"),]
env_ele <- env_all[,1,drop=F]
env_temp <- env_all[,2:12]
env_prec <- env_all[,13:20]
#直接使用原始数据，不做转化，对于群落物种组成数据来讲（因为通常包含很多 0 值），不是很推荐
#rda_result <- rda(phylum~., env, scale = FALSE)
##tb-RDA
#物种数据 Hellinger 预转化（处理包含很多 0 值的群落物种数据时，推荐使用）
phylum_hel <- decostand(phylum, method = 'hellinger')
#环境数据
rda_tb_all <- rda(phylum_hel~., env_all, scale = FALSE)
rda_tb_temp <- rda(phylum_hel~., env_temp, scale = FALSE)
rda_tb_prec <- rda(phylum_hel~., env_prec, scale = FALSE)
#rda_tb_ele <- rda(phylum_hel~., env_ele, scale = FALSE)
#分类样本提取
taxa <- read.table("select_taxa4.txt",header=T,stringsAsFactors = F,sep="\t")
taxa_EA_N <- taxa[which(taxa$Region=="Cen_A" | taxa$Region=="NE_A" ),1]
taxa_EA_S <- taxa[which(taxa$Region=="SW_A" | taxa$Region=="SA"|taxa$Region =="Tibet"),1]
taxa_WA <- taxa[which(taxa$Region=="WA"),1]
taxa_Tibet <- taxa[which(taxa$Region=="Tibet"),1]
taxa_SCA <- taxa[which(taxa$Region=="NW_A" | taxa$Region=="CA"),1]
#taxa_AF <- taxa[which(taxa$Region=="AF"),1]
taxa_EU <- taxa[which(taxa$Region=="EU"),1]
#画不同区域的对应的不同变量的解释度barplot，重复500次----
alltemp1 <- vector()
allprec1 <- vector()
allele1 <- vector()
alltemp2 <- vector()
allprec2 <- vector()
allele2 <- vector()
#100次重复，计算SE
x <- 1
while (x < 100){
  #选择TAXA_new-----
  taxa_north <- taxa_EA_N[sort(sample(c(1:length(taxa_EA_N)),size=20))]
  taxa_south <- taxa_EA_S[sort(sample(c(1:length(taxa_EA_S)),size=20))]
  taxa_WA_25 <- taxa_WA[sort(sample(c(1:length(taxa_WA)),size=20))]
  taxa_EU_25 <- taxa_EU[sort(sample(c(1:length(taxa_EU)),size=20))]
  taxa_SCA_25 <- taxa_SCA[sort(sample(c(1:length(taxa_SCA)),size=20))]
  #taxa_Tibet_25 <- taxa_SCA[sort(sample(c(1:length(taxa_Tibet)),size=10))]
  #筛选样本变异----
  phylum_EA_N <- phylum_hel[taxa_north,]
  phylum_EA_S <- phylum_hel[taxa_south,]
  phylum_WA <- phylum_hel[taxa_WA_25,]
  phylum_SCA <- phylum_hel[taxa_SCA_25,]
  phylum_EU <- phylum_hel[taxa_EU_25,]
  #phylum_Tibet <- phylum_hel[taxa_Tibet_25,]
  #筛选样本的环境（海拔）变量----
  ele_EA_N <- env_ele[taxa_north,,drop=F]
  ele_EA_S <- env_ele[taxa_south,,drop=F]
  ele_WA <- env_ele[taxa_WA_25,,drop=F]
  ele_SCA <- env_ele[taxa_SCA_25,,drop=F]
  ele_EU <- env_ele[taxa_EU_25,,drop=F]
  #ele_Tibet <- env_ele[taxa_Tibet_25,,drop=F]
  #筛选样本的环境（温度）变量----
  temp_EA_N <- env_temp[taxa_north,,drop=F]
  temp_EA_S <- env_temp[taxa_south,,drop=F]
  temp_WA <- env_temp[taxa_WA_25,,drop=F]
  temp_SCA <- env_temp[taxa_SCA_25,,drop=F]
  temp_EU <- env_temp[taxa_EU_25,,drop=F]
  #temp_Tibet <- env_temp[taxa_Tibet_25,,drop=F]
  #筛选样本的环境（降水）变量----
  prec_EA_N <- env_prec[taxa_north,]
  prec_EA_S <- env_prec[taxa_south,]
  prec_WA <- env_prec[taxa_WA_25,]
  prec_SCA <- env_prec[taxa_SCA_25,]
  prec_EU <- env_prec[taxa_EU_25,]
  #prec_Tibet <- env_prec[taxa_Tibet_25,]
  #筛选样本的环境（all）变量----
  env_EA_N <- env_all[taxa_north,]
  env_EA_S <- env_all[taxa_south,]
  env_WA <- env_all[taxa_WA_25,]
  env_SCA <- env_all[taxa_SCA_25,]
  env_EU <- env_all[taxa_EU_25,]
  #env_Tibet <- env_all[taxa_Tibet_25,]
  #RDA分析(EU_25)----
  #EU
  EU_temp_rda <- rda(phylum_EU~., temp_EU, scale = FALSE)
  RsquareAdj(EU_temp_rda)
  EU_prec_rda <- rda(phylum_EU~., prec_EU, scale = FALSE)
  RsquareAdj(EU_prec_rda)
  #EU_all_rda <- rda(phylum_EU~., env_EU, scale = FALSE)
  #RsquareAdj(EU_all_rda)
  EU_ele_rda <- rda(phylum_EU~., ele_EU, scale = FALSE)
  RsquareAdj(EU_ele_rda)
  #RDA分析(SCA_25)----
  #SCA
  SCA_temp_rda <- rda(phylum_SCA~., temp_SCA, scale = FALSE)
  RsquareAdj(SCA_temp_rda)
  SCA_prec_rda <- rda(phylum_SCA~., prec_SCA, scale = FALSE)
  RsquareAdj(SCA_prec_rda)
  #SCA_all_rda <- rda(phylum_SCA~., env_SCA, scale = FALSE)
  #RsquareAdj(SCA_all_rda)
  SCA_ele_rda <- rda(phylum_SCA~., ele_SCA, scale = FALSE)
  RsquareAdj(SCA_ele_rda)
  #RDA分析(WA_25)----
  #WA
  WA_temp_rda <- rda(phylum_WA~., temp_WA, scale = FALSE)
  RsquareAdj(WA_temp_rda)
  WA_prec_rda <- rda(phylum_WA~., prec_WA, scale = FALSE)
  RsquareAdj(WA_prec_rda)
  #WA_all_rda <- rda(phylum_WA~., env_WA, scale = FALSE)
  #RsquareAdj(WA_all_rda)
  WA_ele_rda <- rda(phylum_WA~., ele_WA, scale = FALSE)
  RsquareAdj(WA_ele_rda)
  #RDA分析(north_25)----
  #north
  north_temp_rda <- rda(phylum_EA_N~., temp_EA_N, scale = FALSE)
  RsquareAdj(north_temp_rda)
  north_prec_rda <- rda(phylum_EA_N~., prec_EA_N, scale = FALSE)
  RsquareAdj(north_prec_rda)
  #north_all_rda <- rda(phylum_EA_N~., env_EA_N, scale = FALSE)
  #RsquareAdj(north_all_rda)
  north_ele_rda <- rda(phylum_EA_N~., ele_EA_N, scale = FALSE)
  RsquareAdj(north_ele_rda)
  #RDA分析(south_25)----
  #South
  South_temp_rda <- rda(phylum_EA_S~., temp_EA_S, scale = FALSE)
  RsquareAdj(South_temp_rda)
  South_prec_rda <- rda(phylum_EA_S~., prec_EA_S, scale = FALSE)
  RsquareAdj(South_prec_rda)
  South_ele_rda <- rda(phylum_EA_S~., ele_EA_S, scale = FALSE)
  RsquareAdj(South_ele_rda)
  #Tibet
  #Tibet_temp_rda <- rda(phylum_Tibet~., temp_Tibet, scale = FALSE)
  #RsquareAdj(Tibet_temp_rda)
  #Tibet_prec_rda <- rda(phylum_Tibet~., prec_Tibet, scale = FALSE)
  #RsquareAdj(Tibet_prec_rda)
  #Tibet_all_rda <- rda(phylum_Tibet~., env_Tibet, scale = FALSE)
  #RsquareAdj(Tibet_all_rda)
  #Tibet_ele_rda <- rda(phylum_Tibet~., ele_Tibet, scale = FALSE)
  #RsquareAdj(Tibet_ele_rda)
  #生成画图输入文件----
  #Rsq
  #tempName1 <- c(as.numeric(RsquareAdj(WA_temp_rda)[1]),as.numeric(RsquareAdj(EU_temp_rda)[1]),as.numeric(RsquareAdj(SCA_temp_rda)[1]),as.numeric(RsquareAdj(north_temp_rda)[1]),as.numeric(RsquareAdj(South_temp_rda)[1]))
  #precName1 <- c(as.numeric(RsquareAdj(WA_prec_rda)[1]),as.numeric(RsquareAdj(EU_prec_rda)[1]),as.numeric(RsquareAdj(SCA_prec_rda)[1]),as.numeric(RsquareAdj(north_prec_rda)[1]),as.numeric(RsquareAdj(South_prec_rda)[1]))
  #eleName1 <- c(as.numeric(RsquareAdj(WA_ele_rda)[1]),as.numeric(RsquareAdj(EU_ele_rda)[1]),as.numeric(RsquareAdj(SCA_ele_rda)[1]),as.numeric(RsquareAdj(north_ele_rda)[1]),as.numeric(RsquareAdj(South_ele_rda)[1]))
  #alltemp1 <- cbind(alltemp1,tempName1)
  #allprec1 <- cbind(allprec1,precName1)
  #allele1 <- cbind(allele1,eleName1)
  #Adjust Rsq
  tempName2 <- c(as.numeric(RsquareAdj(WA_temp_rda)[2]),as.numeric(RsquareAdj(EU_temp_rda)[2]),as.numeric(RsquareAdj(SCA_temp_rda)[2]),as.numeric(RsquareAdj(north_temp_rda)[2]),as.numeric(RsquareAdj(South_temp_rda)[2]))
  precName2 <- c(as.numeric(RsquareAdj(WA_prec_rda)[2]),as.numeric(RsquareAdj(EU_prec_rda)[2]),as.numeric(RsquareAdj(SCA_prec_rda)[2]),as.numeric(RsquareAdj(north_prec_rda)[2]),as.numeric(RsquareAdj(South_prec_rda)[2]))
  eleName2 <- c(as.numeric(RsquareAdj(WA_ele_rda)[2]),as.numeric(RsquareAdj(EU_ele_rda)[2]),as.numeric(RsquareAdj(SCA_ele_rda)[2]),as.numeric(RsquareAdj(north_ele_rda)[2]),as.numeric(RsquareAdj(South_ele_rda)[2]))
  alltemp2 <- cbind(alltemp2,tempName2)
  allprec2 <- cbind(allprec2,precName2)
  allele2 <- cbind(allele2,eleName2)
  x <- x+1
}
#Rsq
#rownames(alltemp1) <- c("WA","EU","SCA","EA_N","EA_S")
#rownames(allprec1) <- c("WA","EU","SCA","EA_N","EA_S")
#Rsq <- cbind(apply(alltemp1,1,mean),apply(allprec1,1,mean),apply(alltemp1,1,sd),apply(allprec1,1,sd))
#colnames(Rsq)<- c("temp_mean","prec_mead","temp_sd","prec_sd")
#write.table(Rsq, "RDA_Rsq.txt", sep="\t")

#Adjust Rsq(后续使用)
rownames(alltemp2) <- c("WA","EU","CA","EA","SA")
rownames(allprec2) <- c("WA","EU","CA","EA","SA")
rownames(allele2) <- c("WA","EU","CA","EA","SA")

AdjRsq <- cbind(apply(alltemp2,1,mean),apply(allprec2,1,mean),apply(allele2,1,mean),apply(alltemp2,1,sd),apply(allprec2,1,sd),apply(allele2,1,sd))
colnames(AdjRsq)<- c("temp_mean","prec_mean","ele_mean","temp_sd","prec_sd","ele_sd")
write.table(AdjRsq, "RDA_AdjRsq_AllRegion.txt", row.names = T,sep="\t",col.names = T,quote=F)
#修改RDA_AdjRsq_new.txt的格式，接下来画图

color2 <- brewer.pal(n = 3, name = "Accent")
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/V3/")
AdjRsq <- read.table("RDA_AdjRsq_AllRegion.txt", header=T,sep="\t")
AdjRsq <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/V2/RDA_AdjRsq_AllRegion_noTibet.txt", header=T,sep="\t")
AdjRsq <- AdjRsq[1:10,]
AdjRsq <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/V2/RDA_AdjRsq_AllRegion.txt", header=T,sep="\t")
AdjRsq <- AdjRsq[1:10,]
AdjRsq <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/V2/RDA_AdjRsq_GeneRegion_noTibet.txt", header=T,sep="\t")
AdjRsq <- AdjRsq[1:10,]
AdjRsq <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/V2/RDA_AdjRsq_GeneRegion.txt", header=T,sep="\t")
AdjRsq <- AdjRsq[1:10,]

AdjRsq$Region = factor(AdjRsq$Region, levels=c('EU','WA','CA','EA','SA'))
AdjRsq$Type = factor(AdjRsq$Type, levels=c("Temperature","Precipitation"))
ggplot(AdjRsq, aes(x=Region, y=Adjusted.R.squared, group=Type,fill=Type)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Adjusted.R.squared-SE, ymax=Adjusted.R.squared+SE),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+
  theme_classic()+
  ylab("Proportion variance explained")+
  scale_fill_manual(values = color2) 

##整体点线图
library(vegan)
library(RColorBrewer)
library(ggplot2)
phylum <- read.delim('shuf_20k_noM.txt',  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

row.names(phylum) <- c(1:20000)
phylum <- data.frame(t(phylum))
colnames(phylum) <- c(1:20000)
phylum <- phylum[which(rownames(phylum)!="TW095"),]
env <- read.delim('select_bio2.txt', row.names = 1, header=T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
env_all <- data.frame(env[,1:20])
env_all <- env_all[which(rownames(env_all)!="TW095"),]
phylum_hel <- decostand(phylum, method = 'hellinger')
rda_tb_all <- rda(phylum_hel~., env_all, scale = FALSE)
rda_tb.scaling1 <- summary(rda_tb_all, scaling = 2)
#rda_tb.scaling1
label <- read.table("select_taxa4.txt", header=T, stringsAsFactors = F)
DeleteName <-  label[is.na(label$RDA_Region),1]
label$RDA_Region <- as.factor(label$RDA_Region)
rda_tb_forward_r.site <- data.frame(rda_tb.scaling1$sites)[1:2]
rda_tb_forward_r.env <- data.frame(rda_tb.scaling1$biplot)[1:2]
#group=read.table("group.txt",header = TRUE,row.names = 1,)
rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
#rda_tb_forward_r.site <- merge(rda_tb_forward_r.site, group, by = 'row.names')
rda_tb_forward_r.env$sample <- rownames(rda_tb_forward_r.env)

F <- rda_tb_forward_r.site[!rownames(rda_tb_forward_r.site) %in% DeleteName,]
F$ID <- rownames(F)
all <- merge(F,label,by="ID")
pdf("RDA.pdf",width=10,height=8)
p <- ggplot(all, aes(RDA1, RDA2,color=RDA_Region)) +
  geom_point( size=3) +
  scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
  theme(panel.border = element_blank()) +
  labs(x = 'RDA1 (31.7%)', y = 'RDA2 (13.9%)') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb_forward_r.env, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.4, 'cm')), size = 1, color = 'brown',alpha=0.5) +
  geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'brown', size = 6)+
  guides(fill=guide_legend(title=NULL))
print(p)
dev.off()

#各变量贡献图
library(vegan)
library(RColorBrewer)
library(ggplot2)
phylum <- read.delim('shuf_20k_noM.txt',  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

row.names(phylum) <- c(1:20000)
phylum <- data.frame(t(phylum))
colnames(phylum) <- c(1:20000)
phylum <- phylum[which(rownames(phylum)!="TW095"),]
env <- read.delim('select_bio2.txt', row.names = 1, header=T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
env_all <- data.frame(env[,1:20])
env_all <- env_all[which(rownames(env_all)!="TW095"),]
phylum_hel <- decostand(phylum, method = 'hellinger')
rda_tb_all <- rda(phylum_hel~., env_all, scale = FALSE)
rda_tb.scaling1 <- summary(rda_tb_all, scaling = 2)

data <- env_all$temp1
temp1 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp2
temp2 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp3
temp3 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp4
temp4 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp5
temp5 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp6
temp6 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp7
temp7 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp8
temp8 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp9
temp9 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp10
temp10 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$temp11
temp11 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$prec1
prec1 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$prec2
prec2 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$prec3
prec3 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$prec4
prec4 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$prec5
prec5 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$prec6
prec6 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$prec7
prec7 <- rda(phylum_hel~data, scale = FALSE)
data <- env_all$prec8
prec8 <- rda(phylum_hel~data, scale = FALSE)

data <- env_all$elev
elev <- rda(phylum_hel~data, scale = FALSE)

tempName2 <- c(as.numeric(RsquareAdj(elev)[2]),
               as.numeric(RsquareAdj(temp1)[2]),
               as.numeric(RsquareAdj(temp2)[2]),
               as.numeric(RsquareAdj(temp3)[2]),
               as.numeric(RsquareAdj(temp4)[2]),
               as.numeric(RsquareAdj(temp5)[2]),
               as.numeric(RsquareAdj(temp6)[2]),
               as.numeric(RsquareAdj(temp7)[2]),
               as.numeric(RsquareAdj(temp8)[2]),
               as.numeric(RsquareAdj(temp9)[2]),
               as.numeric(RsquareAdj(temp10)[2]),
               as.numeric(RsquareAdj(temp11)[2]),
               as.numeric(RsquareAdj(prec1)[2]),
               as.numeric(RsquareAdj(prec2)[2]),
               as.numeric(RsquareAdj(prec3)[2]),
               as.numeric(RsquareAdj(prec4)[2]),
               as.numeric(RsquareAdj(prec5)[2]),
               as.numeric(RsquareAdj(prec6)[2]),
               as.numeric(RsquareAdj(prec7)[2]),
               as.numeric(RsquareAdj(prec8)[2])
)
name <- colnames(env_all)
data <- data.frame(x=name,y=tempName2)
#data <- read.table("propotion.txt",header=T,stringsAsFactors = F)
rownames(data) <- c("Elevation",paste("BIO",seq(1,19),sep=""))
data$x <- row.names(data)
data$type=NA
data[1,3] <- "Elevation"
data[2:12,3] <- "Temperature"
data[13:20,3] <- "Precipitation"
data_order <- data[order(data$y),]
data_order$name <- factor(data_order$x,levels = unique(data_order$x))
data_order %>%
  mutate(name = fct_reorder(name, desc(y))) %>%
  ggplot(aes(x=name, y=y,fill=type,group=type)) +
  geom_bar(stat="identity", width=0.5) +
  coord_flip() +
  xlab("") + ylab("Propotion variance explained")+
  scale_fill_manual(values=color2)+
  theme_classic()

#分区点线图
#203:/data1/home/yafei/003_Project3/RDA
library(vegan)
library(RColorBrewer)
library(ggplot2)
phylum <- read.delim('shuf_20k_noM.txt',  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

row.names(phylum) <- c(1:20000)
phylum <- data.frame(t(phylum))
colnames(phylum) <- c(1:20000)
phylum <- phylum[which(rownames(phylum)!="TW095"),]
env <- read.delim('select_bio2.txt', row.names = 1, header=T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

env_all <- data.frame(env[,1:20])
env_all <- env_all[which(rownames(env_all)!="TW095"),]

env_ele <- env_all[,1,drop=F]
env_temp <- env_all[,2:12]
env_prec <- env_all[,13:20]
phylum_hel <- decostand(phylum, method = 'hellinger')
rda_tb_all <- rda(phylum_hel~., env_all, scale = FALSE)
rda_tb_temp <- rda(phylum_hel~., env_temp, scale = FALSE)
rda_tb_prec <- rda(phylum_hel~., env_prec, scale = FALSE)
taxa <- read.table("select_taxa4.txt",header=T,stringsAsFactors = F,sep="\t")
taxa_EA_N <- taxa[which(taxa$Region=="Cen_A" | taxa$Region=="NE_A" ),1]
taxa_EA_S <- taxa[which(taxa$Region=="SW_A" | taxa$Region=="SA"|taxa$Region =="Tibet"),1]
taxa_WA <- taxa[which(taxa$Region=="WA"),1]
taxa_Tibet <- taxa[which(taxa$Region=="Tibet"),1]
taxa_SCA <- taxa[which(taxa$Region=="NW_A" | taxa$Region=="CA"),1]
taxa_EU <- taxa[which(taxa$Region=="EU"),1]
alltemp1 <- vector()
allprec1 <- vector()
allele1 <- vector()
alltemp2 <- vector()
allprec2 <- vector()
allele2 <- vector()
taxa_north <- taxa_EA_N[sort(sample(c(1:length(taxa_EA_N)),size=20))]
taxa_south <- taxa_EA_S[sort(sample(c(1:length(taxa_EA_S)),size=20))]
taxa_WA_25 <- taxa_WA[sort(sample(c(1:length(taxa_WA)),size=20))]
taxa_EU_25 <- taxa_EU[sort(sample(c(1:length(taxa_EU)),size=20))]
taxa_SCA_25 <- taxa_SCA[sort(sample(c(1:length(taxa_SCA)),size=20))]
phylum_EA_N <- phylum_hel[taxa_north,]
phylum_EA_S <- phylum_hel[taxa_south,]
phylum_WA <- phylum_hel[taxa_WA_25,]
phylum_SCA <- phylum_hel[taxa_SCA_25,]
phylum_EU <- phylum_hel[taxa_EU_25,]
ele_EA_N <- env_ele[taxa_north,,drop=F]
ele_EA_S <- env_ele[taxa_south,,drop=F]
ele_WA <- env_ele[taxa_WA_25,,drop=F]
ele_SCA <- env_ele[taxa_SCA_25,,drop=F]
ele_EU <- env_ele[taxa_EU_25,,drop=F]
temp_EA_N <- env_temp[taxa_north,,drop=F]
temp_EA_S <- env_temp[taxa_south,,drop=F]
temp_WA <- env_temp[taxa_WA_25,,drop=F]
temp_SCA <- env_temp[taxa_SCA_25,,drop=F]
temp_EU <- env_temp[taxa_EU_25,,drop=F]
prec_EA_N <- env_prec[taxa_north,]
prec_EA_S <- env_prec[taxa_south,]
prec_WA <- env_prec[taxa_WA_25,]
prec_SCA <- env_prec[taxa_SCA_25,]
prec_EU <- env_prec[taxa_EU_25,]
env_EA_N <- env_all[taxa_north,]
env_EA_S <- env_all[taxa_south,]
env_WA <- env_all[taxa_WA_25,]
env_SCA <- env_all[taxa_SCA_25,]
env_EU <- env_all[taxa_EU_25,]
EU_temp_rda <- rda(phylum_EU~., temp_EU, scale = FALSE)
EU_prec_rda <- rda(phylum_EU~., prec_EU, scale = FALSE)
EU_ele_rda <- rda(phylum_EU~., ele_EU, scale = FALSE)
SCA_temp_rda <- rda(phylum_SCA~., temp_SCA, scale = FALSE)
SCA_prec_rda <- rda(phylum_SCA~., prec_SCA, scale = FALSE)
SCA_ele_rda <- rda(phylum_SCA~., ele_SCA, scale = FALSE)
WA_temp_rda <- rda(phylum_WA~., temp_WA, scale = FALSE)
WA_prec_rda <- rda(phylum_WA~., prec_WA, scale = FALSE)
WA_ele_rda <- rda(phylum_WA~., ele_WA, scale = FALSE)
north_temp_rda <- rda(phylum_EA_N~., temp_EA_N, scale = FALSE)
north_prec_rda <- rda(phylum_EA_N~., prec_EA_N, scale = FALSE)
north_ele_rda <- rda(phylum_EA_N~., ele_EA_N, scale = FALSE)
South_temp_rda <- rda(phylum_EA_S~., temp_EA_S, scale = FALSE)
South_prec_rda <- rda(phylum_EA_S~., prec_EA_S, scale = FALSE)
South_ele_rda <- rda(phylum_EA_S~., ele_EA_S, scale = FALSE)

tit <- c("EU_temp", "WA_temp", "CA_temp", "EA_temp", "SA_temp","EU_prec", "WA_prec", "CA_prec", "EA_prec", "SA_prec")
all <- list(EU_temp_rda, WA_temp_rda, SCA_temp_rda, north_temp_rda, South_temp_rda,EU_prec_rda, WA_prec_rda, SCA_prec_rda, north_prec_rda, South_prec_rda)
pdf("RDA_regions.pdf")
for(i in c(1:length(all))){
  exp_by_x <- (as.list(all[[i]]$CCA$eig)$RDA1)/(all[[i]]$tot.chi) * 100
  exp_by_y <- (as.list(all[[i]]$CCA$eig)$RDA2)/(all[[i]]$tot.chi) * 100
  rda_tb.scaling1 <- summary(all[[i]], scaling = 2)
  rda_tb_forward_r.site <- data.frame(rda_tb.scaling1$sites)[1:2]
  rda_tb_forward_r.env <- data.frame(rda_tb.scaling1$biplot)[1:2]
  rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
  rda_tb_forward_r.env$sample <- rownames(rda_tb_forward_r.env)
  F <- rda_tb_forward_r.site
  p <- ggplot(F, aes(RDA1, RDA2)) +
    geom_point(size=3,shape = 21,col="red") +
    #scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
    theme_classic()+
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
    theme(panel.border = element_blank()) +
    xlab(paste('RDA1 (', round(exp_by_x, 2), '%)', sep = '')) +
    ylab(paste('RDA2 (', round(exp_by_y, 2), '%)', sep = '')) +
    geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
    geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    geom_segment(data = rda_tb_forward_r.env, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.3, 'cm')), size = 1, color = 'blue',alpha=0.5) +
    geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'blue', size = 5)+
    guides(fill=guide_legend(title=NULL))+
    ggtitle(tit[i])+
    theme(plot.title = element_text(color="black", size=20, face="bold"),legend.text = element_text(size=20),legend.title=element_blank(),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text = element_blank()) +   ## 删去刻度标签
    theme(axis.ticks = element_blank())   ## 删去刻度线
  print(p)
}
dev.off()

library(gridExtra)

tit <- c("EU_temp", "WA_temp", "CA_temp", "EA_temp", "SA_temp")
all <- list(EU_temp_rda, WA_temp_rda, SCA_temp_rda, north_temp_rda, South_temp_rda)

p <- list()
for(i in c(1:length(all))){
  exp_by_x <- (as.list(all[[i]]$CCA$eig)$RDA1)/(all[[i]]$tot.chi) * 100
  exp_by_y <- (as.list(all[[i]]$CCA$eig)$RDA2)/(all[[i]]$tot.chi) * 100
  rda_tb.scaling1 <- summary(all[[i]], scaling = 2)
  rda_tb_forward_r.site <- data.frame(rda_tb.scaling1$sites)[1:2]
  rda_tb_forward_r.env <- data.frame(rda_tb.scaling1$biplot)[1:2]
  rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
  rda_tb_forward_r.env$sample <- rownames(rda_tb_forward_r.env)
  F <- rda_tb_forward_r.site
  p[[i]] <- ggplot(F, aes(RDA1, RDA2)) +
    geom_point(size=3,shape = 21,col="red") +
    #scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
    theme_classic()+
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
    theme(panel.border = element_blank()) +
    xlab(paste('RDA1 (', round(exp_by_x, 2), '%)', sep = '')) +
    ylab(paste('RDA2 (', round(exp_by_y, 2), '%)', sep = '')) +
    geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
    geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    geom_segment(data = rda_tb_forward_r.env, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.3, 'cm')), size = 1, color = 'blue',alpha=0.5) +
    geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'blue', size = 5)+
    guides(fill=guide_legend(title=NULL))+
    ggtitle(tit[i])+
    theme(plot.title = element_text(color="black", size=20, face="bold"),legend.text = element_text(size=20),legend.title=element_blank(),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text = element_blank()) +   ## 删去刻度标签
    theme(axis.ticks = element_blank())   ## 删去刻度线
}
pdf("RDA_regions_temp.pdf",width=12,height=8)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],nrow=2)
dev.off()

tit <- c("EU_prec", "WA_prec", "CA_prec", "EA_prec", "SA_prec")
all <- list(EU_prec_rda, WA_prec_rda, SCA_prec_rda, north_prec_rda, South_prec_rda)
p <- list()
for(i in c(1:length(all))){
  exp_by_x <- (as.list(all[[i]]$CCA$eig)$RDA1)/(all[[i]]$tot.chi) * 100
  exp_by_y <- (as.list(all[[i]]$CCA$eig)$RDA2)/(all[[i]]$tot.chi) * 100
  rda_tb.scaling1 <- summary(all[[i]], scaling = 2)
  rda_tb_forward_r.site <- data.frame(rda_tb.scaling1$sites)[1:2]
  rda_tb_forward_r.env <- data.frame(rda_tb.scaling1$biplot)[1:2]
  rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
  rda_tb_forward_r.env$sample <- rownames(rda_tb_forward_r.env)
  F <- rda_tb_forward_r.site
  p[[i]] <- ggplot(F, aes(RDA1, RDA2)) +
    geom_point(size=3,shape = 21,col="red") +
    #scale_color_manual(values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
    theme_classic()+
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
    theme(panel.border = element_blank()) +
    xlab(paste('RDA1 (', round(exp_by_x, 2), '%)', sep = '')) +
    ylab(paste('RDA2 (', round(exp_by_y, 2), '%)', sep = '')) +
    geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
    geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
    geom_segment(data = rda_tb_forward_r.env, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.3, 'cm')), size = 1, color = 'blue',alpha=0.5) +
    geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'blue', size = 5)+
    guides(fill=guide_legend(title=NULL))+
    ggtitle(tit[i])+
    theme(plot.title = element_text(color="black", size=20, face="bold"),legend.text = element_text(size=20),legend.title=element_blank(),axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    theme(axis.text = element_blank()) +   ## 删去刻度标签
    theme(axis.ticks = element_blank())   ## 删去刻度线
}

pdf("RDA_regions_prec.pdf",width=12,height=8)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],nrow=2)
dev.off()




