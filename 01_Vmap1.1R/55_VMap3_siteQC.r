#site QC
#总体20k
for i in `ls *siteQCfile.txt`; do sed '1d' $i | shuf -n 500 >> site.shuf20k.txt; done
sed -i '1i Chr\tPos\tHeterozygousProportion\tMissingRate\tMaf' site.shuf20k.txt 
#各染色体2k
for i in `ls *siteQCfile.txt`; do sed '1d' $i | shuf -n 5000 | sed '1i Chr\tPos\tHeterozygousProportion\tMissingRate\tMaf' > ${i::-15}.shuf5k.txt; done
mv chr001.shuf5k.txt chr01.shuf5k.txt
mv chr002.shuf5k.txt chr02.shuf5k.txt
mv chr003.shuf5k.txt chr03.shuf5k.txt
mv chr004.shuf5k.txt chr04.shuf5k.txt
mv chr005.shuf5k.txt chr05.shuf5k.txt
mv chr006.shuf5k.txt chr06.shuf5k.txt
mv chr007.shuf5k.txt chr07.shuf5k.txt
mv chr008.shuf5k.txt chr08.shuf5k.txt
mv chr009.shuf5k.txt chr09.shuf5k.txt
count=1
cons=($(seq 1 2 42))
for cons in ${cons[@]}
do
let a=$cons+1
cat chr0${cons}.shuf5k.txt chr0${a}.shuf5k.txt > ${count}_site.shuf5k.txt
let count=${count}+1
done


#总体样本的统计

cat *taxaQCfile.txt| grep -v "Taxa" | sort | datamash -g 1 mean 2 mean 3 | sed '1i Taxa\tHeterozygousProportion\tMissRate' > all_taxa.txt




#各自染色体样本的统计
mv chr001_taxaQCfile.txt chr01_taxaQCfile.txt
mv chr002_taxaQCfile.txt chr02_taxaQCfile.txt
mv chr003_taxaQCfile.txt chr03_taxaQCfile.txt
mv chr004_taxaQCfile.txt chr04_taxaQCfile.txt
mv chr005_taxaQCfile.txt chr05_taxaQCfile.txt
mv chr006_taxaQCfile.txt chr06_taxaQCfile.txt
mv chr007_taxaQCfile.txt chr07_taxaQCfile.txt
mv chr008_taxaQCfile.txt chr08_taxaQCfile.txt
mv chr009_taxaQCfile.txt chr09_taxaQCfile.txt


count=1
cons=($(seq 1 2 42))
for cons in ${cons[@]}
do
let a=$cons+1
cat chr0${cons}_taxaQCfile.txt chr0${a}_taxaQCfile.txt | grep -v "Taxa" | sort | datamash -g 1 mean 2 mean 3 | sed '1i Taxa\tHeterozygousProportion\tMissRate' > ${count}_taxa.txt
let count=${count}+1
done

#画图
library(ggplot2)
library(gridExtra)
p <- list()
for (i in c(1:21)){
  name <- paste("chr",i,sep="")
  filename <- paste(i,"_site.shuf5k.txt",sep="")
  data <- read.table(filename,header=T,stringsAsFactors = F)
  data$HeterozygousProportion <- as.numeric(data$HeterozygousProportion)
  data$MissingRate <- as.numeric(data$MissingRate)
  data$Maf <- as.numeric(data$Maf)
  p[[i]] <- ggplot(data=data,aes(x = HeterozygousProportion))+
    geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
    xlab("HeterozygousProportion")+
    ylab("")+
    scale_x_continuous(limits = c(0,1))+
    #scale_y_continuous(limits = c(0,7))+
    ggtitle(name)+
    theme(
      plot.title = element_text(color="red", size=20, face="bold.italic"),
      axis.title = element_text(size=12,face="bold"),
      axis.text = element_text(face="bold",size=8,color = "black"),
      axis.ticks= element_line(size=0.8,colour = "black"),
      panel.grid =element_blank(),
      panel.border = element_rect(fill=NA,size = 0.8),
      panel.background = element_blank())
}

pdf("A_site_het.pdf",width = 9,height = 9)
grid.arrange(p[[1]],p[[4]],p[[7]],p[[10]],p[[13]],p[[16]],p[[19]],nrow=3)
dev.off()
pdf("B_site_het.pdf",width = 9,height = 9)
grid.arrange(p[[2]],p[[5]],p[[8]],p[[11]],p[[14]],p[[17]],p[[20]],nrow=3)
dev.off()
pdf("D_site_het.pdf",width = 9,height = 9)
grid.arrange(p[[3]],p[[6]],p[[9]],p[[12]],p[[15]],p[[18]],p[[21]],nrow=3)
dev.off()


#Miss
p <- list()
for (i in c(1:21)){
  name <- paste("chr",i,sep="")
  filename <- paste(i,"_site.shuf5k.txt",sep="")
  data <- read.table(filename,header=T,stringsAsFactors = F)
  data$HeterozygousProportion <- as.numeric(data$HeterozygousProportion)
  data$MissingRate <- as.numeric(data$MissingRate)
  data$Maf <- as.numeric(data$Maf)
  p[[i]] <- ggplot(data=data,aes(x = MissingRate))+
    geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
    xlab("MissingRate")+
    ylab("")+
    scale_x_continuous(limits = c(0,1))+
    #scale_y_continuous(limits = c(0,7))+
    ggtitle(name)+
    theme(
      plot.title = element_text(color="red", size=20, face="bold.italic"),
      axis.title = element_text(size=12,face="bold"),
      axis.text = element_text(face="bold",size=8,color = "black"),
      axis.ticks= element_line(size=0.8,colour = "black"),
      panel.grid =element_blank(),
      panel.border = element_rect(fill=NA,size = 0.8),
      panel.background = element_blank())
}

pdf("A_site_mis.pdf",width = 9,height = 9)
grid.arrange(p[[1]],p[[4]],p[[7]],p[[10]],p[[13]],p[[16]],p[[19]],nrow=3)
dev.off()
pdf("B_site_mis.pdf",width = 9,height = 9)
grid.arrange(p[[2]],p[[5]],p[[8]],p[[11]],p[[14]],p[[17]],p[[20]],nrow=3)
dev.off()
pdf("D_site_mis.pdf",width = 9,height = 9)
grid.arrange(p[[3]],p[[6]],p[[9]],p[[12]],p[[15]],p[[18]],p[[21]],nrow=3)
dev.off()

#Maf
p <- list()
for (i in c(1:21)){
  name <- paste("chr",i,sep="")
  filename <- paste(i,"_site.shuf5k.txt",sep="")
  data <- read.table(filename,header=T,stringsAsFactors = F)
  data$HeterozygousProportion <- as.numeric(data$HeterozygousProportion)
  data$MissingRate <- as.numeric(data$MissingRate)
  data$Maf <- as.numeric(data$Maf)
  p[[i]] <- ggplot(data=data,aes(x = Maf))+
    geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
    xlab("Maf")+
    ylab("")+
    scale_x_continuous(limits = c(0,0.5))+
    #scale_y_continuous(limits = c(0,7))+
    ggtitle(name)+
    theme(
      plot.title = element_text(color="red", size=20, face="bold.italic"),
      axis.title = element_text(size=12,face="bold"),
      axis.text = element_text(face="bold",size=8,color = "black"),
      axis.ticks= element_line(size=0.8,colour = "black"),
      panel.grid =element_blank(),
      panel.border = element_rect(fill=NA,size = 0.8),
      panel.background = element_blank())
}

pdf("A_site_maf.pdf",width = 9,height = 9)
grid.arrange(p[[1]],p[[4]],p[[7]],p[[10]],p[[13]],p[[16]],p[[19]],nrow=3)
dev.off()
pdf("B_site_maf.pdf",width = 9,height = 9)
grid.arrange(p[[2]],p[[5]],p[[8]],p[[11]],p[[14]],p[[17]],p[[20]],nrow=3)
dev.off()
pdf("D_site_maf.pdf",width = 9,height = 9)
grid.arrange(p[[3]],p[[6]],p[[9]],p[[12]],p[[15]],p[[18]],p[[21]],nrow=3)
dev.off()

library(ggplot2)
library(gridExtra)
pdf("all_site.pdf",width = 9,height = 9)

data <- read.table("site.shuf20k.txt",header=T,stringsAsFactors = F)
data$HeterozygousProportion <- as.numeric(data$HeterozygousProportion)
data$MissingRate <- as.numeric(data$MissingRate)
data$Maf <- as.numeric(data$Maf)
ggplot(data=data,aes(x = HeterozygousProportion))+
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  xlab("HeterozygousProportion")+
  ylab("")+
  scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,7))+
  theme(
    plot.title = element_text(color="red", size=20, face="bold.italic"),
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())
ggplot(data=data,aes(x = MissingRate))+
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  xlab("MissingRate")+
  ylab("")+
  scale_x_continuous(limits = c(0,1))+
  #scale_y_continuous(limits = c(0,7))+
  
  theme(
    plot.title = element_text(color="red", size=20, face="bold.italic"),
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())

ggplot(data=data,aes(x = Maf))+
  geom_histogram(binwidth = 0.01, fill = "lightblue", colour = "black")+
  xlab("Maf")+
  ylab("")+
  scale_x_continuous(limits = c(0,0.5))+
  #scale_y_continuous(limits = c(0,7))+
  theme(
    plot.title = element_text(color="red", size=20, face="bold.italic"),
    axis.title = element_text(size=12,face="bold"),
    axis.text = element_text(face="bold",size=8,color = "black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.border = element_rect(fill=NA,size = 0.8),
    panel.background = element_blank())	


dev.off()