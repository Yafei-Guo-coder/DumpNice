#204:/data2/yafei/004_Vmap3/VCF/4M/VCF/Lineage
library(ggplot2)
name <- c("A_cultivar","A_cultivar_noRht","A_cultivar_yesRht","A_domemmer","A_freethresh","A_landrace","A_wildemmer","B_cultivar","B_cultivar_noRht","B_cultivar_yesRht","B_domemmer","B_freethresh","B_landrace","B_wildemmer","D_cultivar","D_cultivar_noRht","D_cultivar_yesRht","D_landrace","D_strangulata")
#tajimaD
filename <- paste(name,".Tajima.D",sep="")
all <- data.frame(CHROM="",	BIN_START="",	N_SNPS="",	TajimaD="",chr="", type="",stringsAsFactors = F)
for(i in c(1:length(filename))){
  data <- read.table(filename[i], header=T, stringsAsFactors = F)
  data$chr <- strsplit(sub('_',':', name[i]),":")[[1]][1]
  data$type <- strsplit(sub('_',':', name[i]),":")[[1]][2]
  all <- rbind(all,data)
}
all <- all[-1,]
filter <- all[which(all$TajimaD != "NaN"),]
filter$TajimaD <- as.numeric(filter$TajimaD)
p <- ggplot(filter, aes(x=type, y=TajimaD,fill=chr),na.rm=TRUE) + 
  #facet_grid(chr~.) +
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  #  scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"))+
  #  scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"), breaks=c(breaks=c("red", "pink", "lightblue","orange","yellow")),labels = c(" SS", " AA", " AABB"," AABBDD"," DD"))+
  #scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"),labels = c(" SS", " AA", " AABB"," AABBDD"," DD"))+
  scale_x_discrete(limits= c("strangulata","wildemmer","domemmer","freethresh","landrace","cultivar","cultivar_noRht","cultivar_yesRht"))+
  theme_bw()+ #背景变为白色
  #B:0.01 D:0.004 AB:0.01
  ylab("TajimaD") + xlab("type") + scale_y_continuous(limits=c(-3,6)) + 
  theme_classic() +
  theme(legend.text = element_text(size = 20),
        legend.key.size=unit(1,'cm'),
        legend.key.width=unit(0.8,'cm'),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 30,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_text(size=20))
pdf("tajimaD_V1.pdf",width =10, height=4)
print(p)
dev.off()

#sites.pi
filename <- paste(name,".sites.pi",sep="")
all <- data.frame(CHROM="",	Start="", Stop="",	PI="", chr="", type="", stringsAsFactors = F)
for(i in c(1:length(filename))){
  data <- read.table(filename[i], header=F, stringsAsFactors = F)
  colnames(data) <- c("CHROM","Start","Stop","PI")
  data$chr <- strsplit(sub('_',':', name[i]),":")[[1]][1]
  data$type <- strsplit(sub('_',':', name[i]),":")[[1]][2]
  all <- rbind(all,data)
}
all <- all[-1,]
#filter <- all[which(all$PI != 0),]
filter <- all[which(all$PI != "NaN"),]
filter$PI <- as.numeric(filter$PI)
p <- ggplot(filter, aes(x=type, y=PI,fill=chr),na.rm=TRUE) + 
  #facet_grid(chr ~ .) +
  geom_boxplot(outlier = FALSE,outlier.shape = NA,outlier.colour = NA,width=0.7,position=position_dodge(0.9),notch=TRUE,notchwidth=0.5,alpha=0.8)+ #绘制箱线图
  #  scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"))+
  #  scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"), breaks=c(breaks=c("red", "pink", "lightblue","orange","yellow")),labels = c(" SS", " AA", " AABB"," AABBDD"," DD"))+
  #scale_fill_manual(name = "Genome", values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#FFD92F"),labels = c(" SS", " AA", " AABB"," AABBDD"," DD"))+
  scale_x_discrete(limits= c("strangulata","wildemmer","domemmer","freethresh","landrace","cultivar","cultivar_noRht","cultivar_yesRht"))+
  theme_bw()+ #背景变为白色
  #B:0.01 D:0.004 AB:0.01
  ylab("PI")+xlab("type")+
  scale_y_continuous(limits=c(0,0.003)) + 
  theme_classic()+
  theme(legend.text = element_text(size = 20),
        legend.key.size=unit(1,'cm'),
        legend.key.width=unit(0.8,'cm'),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",size=20),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size = 30,face="plain"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_text(size=20))
pdf("PI.pdf",width =10, height=4)
print(p)
dev.off()
#计算windows pi

#!/bin/bash
for chr in {001..042}
do
WGS --model depth --type toBed --file /data1/home/yafei/004_Vmap3/Pi/Site/chr${chr}.txt --out chr${chr}_1M.bed --size 1000000 &
done


sleep 3h;
for i in {001..042}
do
sed -i '1d' chr${i}.txt 
done

for chr in {001..042}
do
WGS --model depth --type toBed --file chr${chr}.txt --out chr${chr}_1M.bed --size 1000000 &
done
sleep 10h;
for chr in {1..9}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/publicData/wheat/genotype/VMap/VMap3.2/VMap3.2/chr00${i}_VMap3.2.vcf.gz --file2 chr00${chr}_1M.bed --out chr00${chr}_pi &
  done
done
for chr in {10..42}
do
WGS --chr ${chr} --model diversity --type bedPi --file /data1/publicData/wheat/genotype/VMap/VMap3.2/VMap3.2/chr0${i}_VMap3.2.vcf.gz --file2 chr0${chr}_1M.bed --out chr00${chr}_pi &
done
