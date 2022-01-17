#work directory:
#/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/06_Structure/sample_group.txt
#/data1/home/yafei/003_Project3/bayenv/13pop/candidate_gene
#/data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225
#Rht-B1:
#TraesCS4B02G043100
#chr21 start:30861268 end:30863723
#ppd-D1
#TraesCS2D01G079600
#chr11 start:33952048	end:33956269
cat *cloned.gene | sort | uniq -c | awk '{print $2"\t"$3}' > top5.candidate.txt
awk '{print $2}' top5.candidate.txt  > test
grep -f test gene_v1.1_Lulab.gff3 |awk '{print $1"\t"$4"\t"$5"\t"$9}' |awk -F";" '{print $1}' |sed 's/ID=//' > test2
awk 'NR==FNR{a[$2]=$1;b[$2]=$2} NR!=FNR {if($4 in b) {print $0"\t"a[$4]}}' top5.candidate.txt test2 | awk '{print $1"\t"$2-5000"\t"$3+5000"\t"$4"\t"$5}' |sed '1i Chr\tstart-5k\tend+5k\tID\tName' > gene_region.txt

cat gene_region.txt |while read chr from to ID Name
do
  if [ $chr != "Chr" ];then
    bcftools filter /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225/chr${chr}.E6_Landrace_locate.vcf.gz --regions ${chr}:${from}-${to} > ${ID}-${Name}.vcf
  fi
done
for i in `ls *vcf`
do

vcftools --vcf $i --012
paste out.012.pos <(cat out.012 | datamash transpose | sed '1d' ) > ${i::-4}.vcf.txt
cat out.012.indv | datamash transpose | awk '{print "File\tChr\tPos\t"$0}' > indi.txt
rm out.012*
done
for i in `ls *vcf.txt`
do
awk '{print FILENAME"\t"$0}' $i  | sed 's/.vcf.txt//g' > ${i::-4}.vcf.txt2
done
cat indi.txt *vcf.txt2 > all.pos.txt
rm TraesCS*

library(reshape)
library(ggplot2)
library(RColorBrewer)
require (rworldmap)
setwd("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/06_Structure/")
data <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/01_RDA_plot/225.taxa.txt", header=T,row.names = 1, sep="\t", stringsAsFactors = F)
taxa <- read.table("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/06_Structure/20210829/Cluster_Location/cluster_70.txt",header=T,row.names = 1, stringsAsFactors = F)
for(i in names(table(data$File))){
  sub <- data[which(data$File==i),]
  file <- paste(i,".pdf",sep="")
  pdf(file)
  for( j in c(1:dim(sub)[1])){
    tit <- sub[j,3]
    sample <- sub[j,4:228]
    loc <- taxa[names(sample),4:5]
    loc$type <- as.numeric(sample[1,])
    loc$value <- 1
    wheat_reshape <- cast(loc,Latitude_70+Longitude_70~type) 
    wheat_reshape2 <- as.data.frame(wheat_reshape)
    name <- names(table(loc$type))
    if((j-1)%%4 !=0){
      mapPies(wheat_reshape2,xlim=c(-10,130),ylim=c(10,70),nameX="Longitude_70",nameY="Latitude_70",nameZs=name,symbolSize=1.5,
              zColours=brewer.pal(8, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9],main="tit")
      legend(25,95,box.lty=0,bg="transparent","108 landrace GrowHabbit", col="black")
    }else{
      par(mfrow=c(2,2),oma=c(1,1,1,1), mar=c(0,1,0,1), cex=1)
      mapPies(wheat_reshape2,xlim=c(-10,130),ylim=c(10,70),nameX="Longitude_70",nameY="Latitude_70",nameZs=name,symbolSize=1.5,
              zColours=brewer.pal(12, "Set3")[c(2,3,4,5)],barOrient='vert',oceanCol="white",landCol=brewer.pal(9,"Pastel1")[9], main="tit")
      legend(25,95,box.lty=0,bg="transparent",tit, col="black")
    }
  }
  dev.off()
}
