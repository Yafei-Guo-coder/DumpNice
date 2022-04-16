#vcftools --gzvcf Alineage001_Wild_emmer.vcf.gz --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS'
#66:/data1/home/yafei/004_Vmap3/ASMC

for i in `cat population1.txt`
do
vcftools --gzvcf Alineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Alineage001_${i}.frq &
vcftools --gzvcf Blineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Blineage001_${i}.frq &
vcftools --gzvcf Dlineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Dlineage001_${i}.frq &
done

for i in `cat population2.txt`
do
vcftools --gzvcf Alineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Alineage001_${i}.frq &
vcftools --gzvcf Blineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Blineage001_${i}.frq &
done

for i in `cat population3.txt`
do
vcftools --gzvcf Dlineage001_${i}.vcf.gz  --freq --stdout | awk -F":" '{print $1"\t"$2"\t"$3}'| awk '{if($8<0.51) print $1"\t"$1"-"$2"\t"$5"\t"$7"\t"$8"\t"$4; else print $1"\t"$1"-"$2"\t"$7"\t"$5"\t"$6"\t"$4}' | sed '1c CHR\tSNP\tA1\tA2\tMAF\tNCHROBS' > Dlineage001_${i}.frq &
done

for i in `cat clone_gene_triads_Col1.txt`; do grep -w $i /data1/publicData/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3 | grep gene|head -n 1; done |awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sed 's/ID=//'  > clone_gene.gff

#DRC
data <- read.table("/data1/home/yafei/004_Vmap3/ASMC/input/Bread_wheat/Cultivar/A/Alineage001_Cultivar.header.1-1.00.sumOverPair.addGroup.rmNan.neutral.shuf5k.bed",header=F,stringsAsFactors=F)
sub <- data[,7:75]
y <- as.numeric(apply(sub,2,sum))
x <- c(30,60,90,120,150,180,210,240,270,300,330,360,390,420,520,620,720,820,920,1020,1120,1220,1320,1420,1520,1620,1720,1820,1920,2066,2244,2470,2757,3122,3584,4153,4820,5593,6477,7470,8565,9757,11030,12373,13781,15243,16751,18307,19905,21545,23237,24991,26814,28724,30735,32868,35162,37659,40407,43465,46907,50846,55459,61106,68635,79855,96263,124311)
y <- y[1:68]
y <- y/5000
all <- read.table("/data1/home/yafei/004_Vmap3/ASMC/input/Bread_wheat/out/Alineage001_Cultivar.asmc",header=F,stringsAsFactors=F)
test <- all[,4:71]
test <- test[which(test$V4 > 10),]
test[test<0] <- 0
result <- vector()
for(j in c(1:dim(test)[1])){
  b <- vector()
  for(i in c(1:68)){
    b <- append(b, rep(i,times=test[j,i]))
  }
  re <- t.test(y,b)
  result <- append(result,re$p.value)
}
all <- all[which(all$V4 >10),]
all$P <- result
out <- all[,c(2,3,73,74)]
colnames(out) <- c("Chr","Pos","Pos2","P")
write.table(out,"A_Cultivar.out.txt",quote=F,sep="\t",row.names=F)

#manhuttan
library(qqman)
library(tidyverse)
require(gridExtra)
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/ASMC/")
filename <- c("A_Cultivar.out.txt","B_Cultivar.out.txt","D_Cultivar.out.txt","A_Landrace.out.txt","B_Landrace.out.txt","D_Landrace.out.txt")

p <- list()
for (i in seq(1,6)){
    data <- read.table(filename[i], header=T, stringsAsFactors = F)
    gwasResults2 <- data
    colnames(gwasResults2) <- c("CHR", "BP", "SNP","P")
    gwasResults <- gwasResults2
    don <- gwasResults %>% 
      # Compute chromosome size
      group_by(CHR) %>% 
      summarise(chr_len=max(BP)) %>% 
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
      select(-chr_len) %>%
      # Add this info to the initial dataset
      left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
      # Add a cumulative position of each SNP
      arrange(CHR, BP) %>%
      mutate( BPcum=BP+tot) 
    axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
    p[[i]] <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size = 15),
        axis.title.x=element_text(size = 15),
      )
}
pdf("Landrace.pdf",height = 9,width = 9)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=3)
dev.off()
pdf("Cultivar.pdf",height = 9,width = 9)
grid.arrange(p[[4]],p[[5]],p[[6]],nrow=3)
dev.off()

#ihs

library(qqman)
library(tidyverse)
require(gridExtra)
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/ihh12")
filename <- c("A_ihs.5k.txt","B_ihs.5k.txt","D_ihs.5k.txt")

p <- list()
for (i in seq(1,3)){
  data <- read.table(filename[i], header=T, stringsAsFactors = F)
  gwasResults2 <- data
  colnames(gwasResults2) <- c("CHR", "BP", "P")
  gwasResults <- gwasResults2
  don <- gwasResults %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) 
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2)
  p[[i]] <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(size=15),
      axis.text.y=element_text(size=15),
      axis.title.y=element_text(size = 15),
      axis.title.x=element_text(size = 15),
    )
}
pdf("Landrace.pdf",height = 9,width = 9)
grid.arrange(p[[1]],p[[2]],p[[3]],nrow=3)
dev.off()
pdf("Cultivar.pdf",height = 9,width = 9)
grid.arrange(p[[4]],p[[5]],p[[6]],nrow=3)
dev.off()

