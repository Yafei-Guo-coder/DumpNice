library(qqman)
library(tidyverse)
#20210902 FuGWAS----
#准备文件来自：01_GWAS_gwas.sh
#Manhattan plot
setwd("/Users/guoyafei/Documents/04_FuGWAS/07_气孔导度数据/20211007/Manhattan/logP2")
path <- "/Users/guoyafei/Documents/04_FuGWAS/07_气孔导度数据/20211007/Manhattan/logP2/TXT" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F)})

#标注气孔相关基因在曼哈顿上的位置
#shell:yafei@66:/data1/home/yafei/009_GWAS/WEGA_out/stoma/Manhattan/logP2
grep -f Gene_id.txt /data1/home/yafei/009_GWAS/gene/gene_v1.1_Lulab.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$4-1000000"\t"$5+1000000"\t"$6"\t"$7"\t"$8"\t"$9}' | sed '1i ##gff-version 3' > Related_gene_1M.gff3
grep -f Fu_known.gene /data1/home/yafei/009_GWAS/gene/gene_v1.1_Lulab.gff3 | awk '{print $1"\t"$2"\t"$3"\t"$4-1000000"\t"$5+1000000"\t"$6"\t"$7"\t"$8"\t"$9}' | sed '1i ##gff-version 3' | sort -k1,1n -k4,4n > Fu_gene_1M.gff3
for i in `ls *txt`
do
awk '{print $2"\t"$3-1"\t"$3}' $i | sed '1d ' > ${i::-3}bed
done

for i in `ls *bed`; do bedtools intersect -a ../Related_gene_5k.gff3 -b $i -wb; done | awk '{print $10"\t"$12}'|sed 's/\t/-/' > 5k.snp
for i in `ls *bed`; do bedtools intersect -a ../Related_gene.gff3 -b $i -wb; done | awk '{print $10"\t"$12}'|sed 's/\t/-/' > snp
for i in `ls *bed`; do bedtools intersect -a ../Related_gene_1M.gff3 -b $i -wb; done| awk '{print $10"\t"$12}'|sed 's/\t/-/' > 1M.snp
for i in `ls *bed`; do bedtools intersect -a ../Fu_gene_1M.gff3 -b $i -wb; done| awk '{print $10"\t"$12}'|sed 's/\t/-/' > Fu.1M.snp

#基因上下游5M的位点
snp <- read.table("5M.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]
#基因上下游1M的位点
snp <- read.table("1M.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]
#Fu已知基因上下游1M的位点
snp <- read.table("Fu.1M.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]
#基因上下游3M的位点
snp <- read.table("3M.snp",header=F,stringsAsFactors = F)
highsnp <- snp[,1]
#基因上下游5k的位点
highsnp <- c("4-28222453","5-346668758","14-196323275","14-196324770","14-196327488","14-196329328","14-196329549","14-196330321","14-196330993","14-196331076","14-196331701","18-39729024","21-30861571")
#基因区的位点
highsnp <- c("14-196327488","21-30861571")
#-lop的值大于5的位置
snp <- read.table("logP5_50k.snp",header=F,stringsAsFactors = F)
anno <- read.table("/Users/guoyafei/Documents/04_FuGWAS/07_气孔导度数据/20211007/anno.txt",header=F,stringsAsFactors = F)
snpsOfInterest <- snp[,1]
snpsOfAnno <- anno[,1]
pdf("stoma_Fu_gene_1M.pdf",height = 5,width = 15)
for (i in c(1:20)){
  all <- data[[1]]
  colnames(all) <- c("SNP", "CHR", "BP","P")
  #sub <- data[order(data$P),]
  #sub$logP <- -log10(sub$P)
  #sub2 <- sub[which(sub$logP>2),1:4]
  manhattan(all, col = c("#fdbf6f","#BEBADA"), highlight = highsnp,suggestiveline=FALSE,genomewideline=F,logp=T, ylim=c(2,8))
  #manhattan(all, col = c("#fdbf6f","#fdbf6f"), suggestiveline=FALSE, genomewideline=F, logp=T, ylim=c(2,8))
}
dev.off()

#标注-logP的值在5以上的信号位点在曼哈顿上的位置
#shell:yafei@66:/data1/home/yafei/009_GWAS/WEGA_out/stoma/logP5

#提取-lopP5 bed 上下游50k区域的富集基因
#提取-lopP5 bed 上下游50k的区域
for i in `cat txt.names`
do 
awk '{print $2"\t"$3"\t"$4"\t"$7"\t"(-log($7)/log(10))}' $i | awk '{if($5>5 && $4!="NaN") print $0}' |awk '{print $2"\t"$3-50000"\t"$3+50000}' |sed '1d' > logP5/${i::-15}.50k.bed
done

#提取-lopP5 bed的位点
for i in `cat txt.names`
do 
awk '{print $2"\t"$3"\t"$4"\t"$7"\t"(-log($7)/log(10))}' $i | awk '{if($5>5 && $4!="NaN") print $0}' |awk '{print $2"\t"$3}' |sed '1d' | awk '{print $0"\t"NR}' |sed 's/\t/-/'
done > logP5/logP5.pos.txt

sed 's/-/\t/' logP5.pos.txt | awk '{print $1"\t"$2-1"\t"$2}' > logP5.pos.bed
#把-lopP5 bed 上下游50k的区域在1M之内的合并在一起
for i in {1..7}
do
for j in {A,B,D}
do
bedtools merge -d 1000000 -i ${i}${j}.50k.bed > ${i}${j}.50k.merge.bed  
done
done
for i in `ls *50k.merge.bed`
do
awk '{print $0"\t"FILENAME"\t"NR}' $i >> All_50k_region.bed
done

bedtools intersect -a All_50k_region.bed -b logP5.pos.bed -wo|sort -k1,1n -k2,2n >test.txt
#曼哈顿图注释表
bedtools intersect -a /data1/home/yafei/009_GWAS/gene/gene_v1.1_Lulab.gff3 -b All_50k_region.bed -wb | awk 'split($9, array, ";") {print $1"\t"$4"\t"$5"\t"array[1]"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' | sed '1i gene_chr\tgene_start\tgene_end\tgene_id\tsnpBlock_chr\tsnpBlock_start\tsnpBlock_end\tfileName\tsnpBlock_id'> snpBlock_annotation3.txt
#bedtools intersect -a /data1/home/yafei/009_GWAS/gene/gene_v1.1_Lulab.gff3 -b All_50k_NoMerge.bed -wb | awk 'split($9, array, ";") {print $1"\t"$4"\t"$5"\t"array[1]"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' | sed '1i gene_chr\tgene_start\tgene_end\tgene_id\tsnpBlock_chr\tsnpBlock_start\tsnpBlock_end\tfileName\tsnpBlock_id'> snpBlock_annotation3.txt

#mlm bed
for i in `cat txt.names`
do 
awk '{if($4!="NaN")print $3"\t"$4-1"\t"$4}' $i | sed '1d'  > logP5/${i::-15}.mlm.bed 
done

#intersect
for i in {1..7}
do
for j in {A,B,D}
do
bedtools intersect -a ${i}${j}.50k.bed -b ${i}${j}.mlm.bed -wb
done
done | awk '{print $4"\t"$6}'|sort -k1,1n -k2,2n | uniq | sed 's/\t/-/' > logP5_50k.snp


#QQ plot-----
setwd("/Users/guoyafei/Documents/01_个人项目/05_FuGWAS/07_气孔导度数据/20210928/")
data <- read.table("height_all.mlm.txt",header=F,stringsAsFactors = F)
colnames(data) <- c("SNP", "CHR", "BP","P")
qq_dat <- data.frame(obs=-log10(sort(data$P,decreasing=FALSE)),
                     exp=-log10( ppoints(length(data$P))))
pd_qq <- ggplot(data=qq_dat,aes(exp,obs))+
  geom_point(alpha=0.7,color="#7F7F7FFF")+
  geom_abline(color="#D62728FF")+
  xlab("Expected -log10(P-value)")+
  ylab("Observed -log10(P-value)")+
  scale_x_continuous(limits = c(0,7))+
  scale_y_continuous(limits = c(0,7))+
  theme(
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    #axis.line = element_line(size=0.8,color="black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())

png("height_qq.png")
pd_qq
dev.off()

#----

for (i in c(1:20)){
gwasResults <- data[[i]]
colnames(gwasResults) <- c("SNP", "CHR", "BP","P")
don <- gwasResults %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # ！！！！！！Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate( is_annotate=ifelse(SNP %in% snpsOfAnno, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
a<-i+1
pdf(paste(a,".png",sep=""),height = 5,width = 15)
p <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
print(p)
dev.off()
}


