#shell requirement---------
#计算RDA：环境变量和遗传变异
#工作目录：/data2/yafei/Project3/Vmap1.1/Out/VCF/VmapE6/Landrace/Select_taxa
#1. 提取计算迁徙路径的样本的VCF文件
#vcftools --gzvcf D_Land.vcf.gz --keep Select_taxa.txt --maf 0.0000001 --recode --stdout | bgzip -c > Select_taxa/D_Land_Select.vcf.gz

#2. 提取基因区的VCF文件
#bedtools intersect -a Gene.gff3 -b A_Land_Select.vcf.gz -wb > A_Land_Select_gene.vcf
#cat VCF.header A_Land_Select_gene.vcf | bgzip -c > A_Land_Select_gene.vcf.gz

#3. 合并vcf并随机选取3000个位点
#vcf-concat A_Land_Select_gene.vcf.gz B_Land_Select_gene.vcf.gz D_Land_Select_gene.vcf.gz | bgzip -c > All_gene.vcf.gz
#run_pipeline.pl -Xms10g -Xmx200g -vcf All_gene.vcf.gz -sortPositions -export All_gene.hmp.txt -exportType HapmapDiploid

#sed '1d' All.sort.hmp.txt | shuf -n 5000 > shuf_5000.hmp.txt
#head -n 1 All.sort.hmp.txt > shuf.header
#cat shuf.header shuf_5000.hmp.txt > shuf_5000.hmp.txt2
#mv shuf_5000.hmp.txt2 shuf_5000.hmp.txt
#run_pipeline.pl -SortGenotypeFilePlugin -inputFile shuf_5000.hmp.txt -outputFile shuf_5000.sort.hmp.txt -fileType Hapmap
#run_pipeline.pl -Xmx100g -fork1 -h shuf_5000.sort.hmp.txt -export -exportType VCF

#4. 样本分区：EU, WA, CA, SA, NW_A, SW_A, NE_A, SE_A, Tibet

#5. RDA analysis----
#input files:
#env_table: data frame, row: sample name, col: environment variables
#genotype_table: data frame, row: sample name, col: snp site
#code----
library(vegan)
library(RColorBrewer)
library(ggplot2)
setwd("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment")
#input environment variants file and genetic variants file and RDA analysis----
phylum <- read.delim('All_noMiss_0.05_2000.txt',  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
row.names(phylum) <- c(1:2000)
phylum <- data.frame(t(phylum))
env <- read.delim('select_bio.txt', row.names = 1, header=T,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
env_all <- data.frame(env[,1:20])
env_temp <- env_all[,1:12]
env_prec <- env_all[,13:20]
#直接使用原始数据，不做转化。对于群落物种组成数据来讲（因为通常包含很多 0 值），不是很推荐
#rda_result <- rda(phylum~., env, scale = FALSE)
##tb-RDA
#物种数据 Hellinger 预转化（处理包含很多 0 值的群落物种数据时，推荐使用）
phylum_hel <- decostand(phylum, method = 'hellinger')
#使用全部的环境数据
rda_tb_all <- rda(phylum_hel~., env_all, scale = FALSE)
rda_tb_temp <- rda(phylum_hel~., env_temp, scale = FALSE)
rda_tb_prec <- rda(phylum_hel~., env_prec, scale = FALSE)
rda_tb.scaling1 <- summary(rda_tb_all, scaling = 2)
rda_tb.scaling1

#若只关注局部环境数据，除了在原始表格中修改变量个数外，还可直接在 rda() 中指定
#rda_part <- rda(phylum~elevation+one+two+three+four+five+six+seven+eight+nine+ten+eleven+twelve+thirteen+fourteen+fifteen+sixteen+seventeen+eighteen+nineteen, data = env, scale = FALSE)

#plot point-line graph---------

#color = c("#838B8B", "#8470FF", "#D8BFD8", "#FF6349", "#FFD700") 
color <- brewer.pal(n = 5, name = "Accent")[c(1:5)]
DeleteName <-  c("ZN109","TW004","TW005","TW053","TW054","TW060","TW088","TW109","TW110","XI_29")

label <- read.table("select_taxa2.txt", header=T, stringsAsFactors = F)
label$Region <- as.factor(label$Region)
label$cols = label$Region 
label$cols = gsub("EA-N","#D8BFD8",label$cols)
label$cols = gsub("EU","#838B8B",label$cols)
label$cols = gsub("SCA","#97FFFF",label$cols)
label$cols = gsub("WA","#FFD700",label$cols)
label$cols = gsub("EA-S","#FF6349",label$cols)
#label$cols = gsub("AM","#8470FF",label$cols)
#label$cols = gsub("AF","#E7298A",label$cols)
#rda_tb <- rda(phylum_hel~., env_all, scale = FALSE)
#r2 <- RsquareAdj(rda_tb)
#rda_noadj <- r2$r.squared
#rda_adj <- r2$adj.r.squared
#rda_tb_test <- anova(rda_tb, permutations = 999)
#rda_tb_test_axis <- anova(rda_tb, by = 'axis', permutations = 999)
#rda_tb_test_axis$`Pr(>F)` <- p.adjust(rda_tb_test_axis$`Pr(>F)`, method = 'bonferroni')
#pca_eig <- rda_tb$CA$eig
#pca_eig[pca_eig > mean(pca_eig)]
#n <- length(pca_eig)
#bsm <- data.frame(j=seq(1:n), p = 0)
#bsm$p[1] <- 1/n
#for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
#bsm$p <- 100*bsm$p/n
#vif.cca(rda_tb)
#rda_tb_forward_r <- ordiR2step(rda(phylum~1, env, scale = FALSE), scope = formula(rda_tb), R2scope = rda_adj, direction = 'forward')
#rda_tb_forward_r.scaling1 <- summary(rda_tb_forward_r, scaling = 1)
rda_tb_forward_r.site <- data.frame(rda_tb.scaling1$sites)[1:2]
rda_tb_forward_r.env <- data.frame(rda_tb.scaling1$biplot)[1:2]
#group=read.table("group.txt",header = TRUE,row.names = 1,)
rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
#rda_tb_forward_r.site <- merge(rda_tb_forward_r.site, group, by = 'row.names')
rda_tb_forward_r.env$sample <- rownames(rda_tb_forward_r.env)

F <- rda_tb_forward_r.site[!rownames(rda_tb_forward_r.site) %in% DeleteName,]
label$cols <- as.factor(label$cols)
label$Region <- as.factor(label$Region)
p <- ggplot(F, aes(RDA1, RDA2)) +
  geom_point(aes(color = label$Region), size=3) +
  #stat_ellipse(aes(group=label$cols), level = 0.95, show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = color) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
  theme_classic()+
  theme(panel.grid =element_blank()) +   ## 删去网格线
  #theme(axis.text = element_blank()) +   ## 删去刻度标签
  #theme(axis.ticks = element_blank()) +   ## 删去刻度线
  theme(panel.border = element_blank()) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'))+
  #legend.title = (element_blank(), legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'RDA1 (34.12%)', y = 'RDA2 (11.90%)') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb_forward_r.env, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.3, 'cm')), size = 1, color = 'brown',alpha=0.2) +
  geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'brown', size = 5)+
  #scale_colour_discrete(breaks = c("#838B8B","#FFD700", "#97FFFF", "#D8BFD8", "#FF6349"), labels = c('EU','WA','SCA','EA-N','EA-S'))+
  guides(fill=guide_legend(title=NULL))
#geom_label_repel(aes(label =sample, color = group), size = 3, box.padding = unit(0, 'lines'), show.legend = FALSE)
p

#分类样本提取----
taxa <- read.table("select_taxa3.txt",header=T,stringsAsFactors = F)
taxa_EA_N <- taxa[which(taxa$Region=="EA-N"),1]
taxa_EA_S <- taxa[which(taxa$Region=="EA-S"),1]
taxa_EU <- taxa[which(taxa$Region=="EU"),1]
taxa_SCA <- taxa[which(taxa$Region=="SCA"),1]
taxa_WA <- taxa[which(taxa$Region=="WA"),1]

#taxa_CA <- taxa[which(taxa$RDA=="CA"),1]
taxa_EU <- taxa[which(taxa$RDA=="EU"),1]
taxa_N_A <- taxa[which(taxa$RDA=="N_A"),1]
#taxa_NW_A <- taxa[which(taxa$RDA=="NW_A"),1]
#taxa_SA <- taxa[which(taxa$RDA=="SA"),1]
taxa_SCA <- taxa[which(taxa$RDA=="SCA"),1]
#taxa_SE_A <- taxa[which(taxa$RDA=="SE_A"),1]
taxa_SW_A <- taxa[which(taxa$RDA=="SW_A"),1]
#taxa_Tibet <- taxa[which(taxa$RDA=="Tibet"),1]
taxa_WA <- taxa[which(taxa$RDA=="WA"),1]

#taxa_region = list(taxa_EU,taxa_N_A,taxa_SCA,taxa_SW_A,taxa_WA)
taxa_region = list(taxa_EA_N,taxa_EA_S,taxa_EU,taxa_SCA,taxa_WA)
#区域环境变量提取----

#100次重复，计算SE
alltemp <- vector()
allprec <- vector()
for(i in c(1:length(taxa_region))){
  taxaTemp <- vector()
  taxaPrec <- vector()
  x <- 1
  while (x < 100){
    sample_taxa <- taxa_region[[i]][sort(sample(c(1:length(taxa_region[[i]])),size=17))]
    phylum_taxa <- phylum_hel[sample_taxa,]
    temp_taxa <- env_temp[sample_taxa,]
    prec_taxa <- env_prec[sample_taxa,]
    #env_taxa <- env_all[sample_taxa,]
    taxa_temp_rda <- rda(phylum_taxa~., temp_taxa, scale = FALSE)
    #RsquareAdj(taxa_temp_rda)
    taxa_prec_rda <- rda(phylum_taxa~., prec_taxa, scale = FALSE)
    #RsquareAdj(taxa_prec_rda)
    #taxa_all_rda <- rda(phylum_taxa~., env_taxa, scale = FALSE)
    #RsquareAdj(taxa_all_rda)
    #生成画图输入文件----
    #Rsq
    #tempName1 <- c(as.numeric(RsquareAdj(WA_temp_rda)[1]),as.numeric(RsquareAdj(EU_temp_rda)[1]),as.numeric(RsquareAdj(SCA_temp_rda)[1]),as.numeric(RsquareAdj(north_temp_rda)[1]),as.numeric(RsquareAdj(South_temp_rda)[1]))
    #precName1 <- c(as.numeric(RsquareAdj(WA_prec_rda)[1]),as.numeric(RsquareAdj(EU_prec_rda)[1]),as.numeric(RsquareAdj(SCA_prec_rda)[1]),as.numeric(RsquareAdj(north_prec_rda)[1]),as.numeric(RsquareAdj(South_prec_rda)[1]))
    #alltemp1 <- cbind(alltemp1,tempName1)
    #allprec1 <- cbind(allprec1,precName1)
    #Adjust Rsq
    taxaTemp <- append(taxaTemp,as.numeric(RsquareAdj(taxa_temp_rda)[2]))
    taxaPrec <- append(taxaPrec,as.numeric(RsquareAdj(taxa_prec_rda)[2]))
    x <- x+1
  }
  alltemp <- cbind(alltemp,taxaTemp)
  allprec <- cbind(allprec,taxaPrec)
}

#taxa_name <- c("taxa_EU","taxa_N_A","taxa_SCA","taxa_SW_A","taxa_WA")
taxa_name = c("taxa_EA_N","taxa_EA_S","taxa_EU","taxa_SCA","taxa_WA")
colnames(alltemp) <- taxa_name
colnames(allprec) <- taxa_name
AdjRsq <- cbind(apply(alltemp,2,mean,na.omit=T),apply(allprec,2,mean),apply(alltemp,2,sd),apply(allprec,2,sd))
colnames(AdjRsq)<- c("temp_mean","prec_mean","temp_sd","prec_sd")
write.table(AdjRsq, "RDA_AdjRsq_new.txt", row.names = T,sep="\t",col.names = T,quote=F)

#画各区域RDA分解图----
plot(EU_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(EU_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(EU_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

plot(WA_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(WA_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(WA_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

plot(SCA_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(SCA_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(SCA_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

plot(north_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(north_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(north_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

plot(South_west_temp_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1)
points(South_west_temp_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, cex = 1)
text(South_west_temp_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'brown', cex = 1)

#画整体barplot----
#old version ----
#version 1
#data <- read.table("plot_data.txt", header=T, row.names = 1,stringsAsFactors = F)
#data<-t(as.matrix(data))
#barplot(data, beside = TRUE,
#        col = c("lightblue", "mistyrose"),
#        legend = rownames(data), ylim = c(0, 0.3))

# version 2----
color2 <- brewer.pal(n = 4, name = "Accent")
AdjRsq <- read.table("RDA_AdjRsq_new.txt", header=T,sep="\t")
AdjRsq$Region = factor(AdjRsq$Region, levels=c('EU','WA','SCA','EA_N','EA_SW'))
AdjRsq$Type = factor(AdjRsq$Type, levels=c("Temperature","Precipitation"))
ggplot(AdjRsq, aes(x=Region, y=Mean, group=Type,fill=Type)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+
  theme_classic()+
  scale_fill_manual(values = c("#FDC086","#BEAED4")) 

#version 3----
color2 <- brewer.pal(n = 4, name = "Accent")
AdjRsq <- read.table("RDA_AdjRsq.txt", header=T,sep="\t")
AdjRsq$Region = factor(AdjRsq$Region, levels=c('EU','WA','SCA','EA_N','EA_S'))
AdjRsq$Type = factor(AdjRsq$Type, levels=c("Temperature","Precipitation"))

ggplot(AdjRsq, aes(x=Region, y=Mean, group=Type,fill=Type)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+
  theme_classic()+
  scale_fill_manual(values = c("#FDC086","#BEAED4")) 


