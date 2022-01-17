#XP-CLR
#工作路径：
#xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth
#XPCLR的流程：XPCLR -xpclr ../groupSouth/groupSouthChr${i}.geno ../groupSCA/groupSCAChr${i}.geno ../groupSCA/groupSCAChr${i}.snp South_SCA_10kchr${i} -w1 0.005 500 10000 $i -p1 0.95 &
#这是参数 -w1 0.005 500 10000 -p1 0.95
#下面的是做smooth标准化
#该目录下一共有 9 个结果文件夹，分别是
#XPCLRresult_EU_South
#XPCLRresult_North2_South
#XPCLRresult_Strang_WA
#XPCLRresult_Tibet_South
#XPCLRresult_WA_EU
#XPCLRresult_WA_South
#XPCLRresult_neg_North1_North2
#XPCLRresult_neg_WA_North1
#XPCLRresult_neg_WA_North2
neg_WA_1_2_10kchr6.xpclr.txt
#其中三个是negative contral.

#!/usr/bin/Rscript.R
library(GenWin)
library(dplyr)
Folder <- c("EU_South","North2_South","Strang_WA","Tibet_South","WA_EU","WA_South","neg_North1_North2","neg_WA_North1","neg_WA_North2")
for (num in c(1:9)){
  filePath <- paste("/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/XPCLRresult_",Folder[num],sep="")
  setwd(filePath)
  for(i in c(1:42)){
    file=paste(Folder[num],"_10kchr",i,".xpclr.txt",sep="")
    Data1=read.table(file,sep=" ")
    Data2<-na.omit(Data1)
    Data <- filter(Data2, Data2[,6] != Inf)
    Z=matrix(, nrow = nrow(Data),ncol=1)
    figure=paste("pdfnom_c",i,".pdf",sep="")
    pdf(figure,width=12,height=8)
    for(j in 1:nrow(Data)){
      Z[j,1]=(Data[j,6]-mean(Data$V6))/sd(Data$V6)
    }
    NORM=splineAnalyze(Y=Z,map=Data$V4,smoothness = 2000,plotRaw=F,plotWindows = T,method = 2)
    dev.off()
    normScore=NORM$windowData
    outFile=paste("/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/",Folder[num],"_smooth",i,".txt",sep="")
    write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
  }
}
#Rscript smooth.R
#GenWin(smooth)的结果是这样的"WindowStart" "WindowStop" "SNPcount" "MeanY" "Wstat"，没有染色体，要加上染色体的信息
#工作路径：xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/XPCLRresultSouth_SCA/smooth
#接下来使用02_Xp-clr_SelectGene.sh脚本给smooth结果添加染色体号并进行排序，并且合并成A，B，D lineage.

#permutation
#204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3_permutation
library(GenWin)
library(dplyr)
Folder <- c("WA_CA", "WA_EA", "WA_EU", "WA_SA")
for (num in c(1:4)){
  for (count in c(1:10)){
    filePath <- paste("/data1/home/yafei/003_Project3/XP-CLR/permutation/",Folder[num],"/Xp",count,sep="")
    setwd(filePath)
    for(i in c(1:42)){
      file=paste("xpclr_",Folder[num],i,".xpclr.txt",sep="")
      Data1=read.table(file,sep=" ")
      Data2<-na.omit(data)
      Data <- Data2[which(Data2[,6] != Inf),]
      Z=matrix(, nrow = nrow(Data),ncol=1)
      #figure=paste("pdfnom_c",i,".pdf",sep="")
      #pdf(figure,width=12,height=8)
      m <- mean(Data$V6)
      s <- sd(Data$V6)
      for(j in 1:nrow(Data)){
        Z[j,1]=(Data[j,6]-m)/s
      }
      NORM=splineAnalyze(Y=Z,map=Data$V4,smoothness = 2000,plotRaw=F,plotWindows = F,method = 4)
      #dev.off()
      normScore=NORM$windowData
      outFile=paste("/data1/home/yafei/003_Project3/XP-CLR/permutation/smooth/",Folder[num],"/Xp",count,"_",Folder[num],"_smooth",i,".txt",sep="")
      write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
    }
  }
}

#merge
Xp3_WA_CA_smooth22.txt

Name=(WA_CA, WA_EA, WA_EU, WA_SA)

for num in {1..10}
do
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
sed '1d' Xp${num}_WA_CA_smooth${i}.txt | awk '{print "'$i'""\t"$0}'  
done |sed '/NA/d' | sort -k6,6g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > Xp${num}_WA_CA_smooth_A.txt
for i in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
sed '1d' Xp${num}_WA_CA_smooth${i}.txt | awk '{print "'$i'""\t"$0}'  
done |sed '/NA/d' | sort -k6,6g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > Xp${num}_WA_CA_smooth_B.txt
for i in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
sed '1d' Xp${num}_WA_CA_smooth${i}.txt  | awk '{print "'$i'""\t"$0}'  
done |sed '/NA/d' | sort -k6,6g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > Xp${num}_WA_CA_smooth_D.txt
done

mkdir lineage
mv *A.txt *B.txt *D.txt lineage
cd lineage
mkdir Top5%
for i in `ls *txt`; do  wc -l $i; done | awk '{print $2"\t"$1"\t"$1*0.01"\t"$1*0.05}' | awk -F"[.|\t]" '{print $1"\t"$3"\t"$4"\t"$6}' | awk '{print "tail -n "$4,$1".txt > Top5%/"$1".top5.bed"}'



