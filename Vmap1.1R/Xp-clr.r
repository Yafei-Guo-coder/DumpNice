#XP-CLR
#工作路径：xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/
#XPCLR的流程：XPCLR -xpclr ../groupSouth/groupSouthChr${i}.geno ../groupSCA/groupSCAChr${i}.geno ../groupSCA/groupSCAChr${i}.snp South_SCA_10kchr${i} -w1 0.005 500 10000 $i -p1 0.95 &
#这是参数 -w1 0.005 500 10000 -p1 0.95
#下面的是做smooth标准化

#!/usr/bin/Rscript.R
setwd("/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/XPCLRresultNorth_SCA")
library(GenWin)
library(dplyr)
for(i in c(1:42)){
  file=paste("North_SCA_10kchr",i,".xpclr.txt",sep="")
  Data1=read.table(file,sep=" ")
  Data2<-na.omit(Data1)
  Data <- filter(Data2, Data2[,6] != Inf)
  Z=matrix(, nrow = nrow(Data),ncol=1)
  figure=paste("pdfnom_c",i,".pdf",sep="")
  pdf(figure,width=12,height=8)
  for(j in 1:nrow(Data)){
    Z[j,1]=(Data[j,6]-mean(Data$V6))/sd(Data$V6)
  }
  NORM=splineAnalyze(Y=Z,map=Data$V4,smoothness = 2000,plotRaw=T,plotWindows = T,method = 4)
  dev.off()
  normScore=NORM$windowData
  outFile=paste("/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/XPCLRresultNorth_SCA/smooth",i,".txt",sep="")
  write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
}

setwd("/data1/home/xuebo/Projects/Evo/xpclr/lineage/newfileV9/Agroup/XPCLRresult")
library(GenWin)
library(dplyr)
for(i in c(1,2,7,8,13,14,19,20,25,26,31,32,37,38)){
  file=paste("SpeciesA_100kchr",i,".xpclr.txt",sep="")
  Data1=read.table(file,sep=" ")
  Data2<-na.omit(Data1)
  Data <- filter(Data2, Data2[,6] != Inf)
  Z=matrix(, nrow = nrow(Data),ncol=1)
  figure=paste("pdfnom_c",i,".pdf",sep="")
  pdf(figure,width=12,height=8)
  for(j in 1:nrow(Data)){
    Z[j,1]=(Data[j,6]-mean(Data$V6))/sd(Data$V6)
  }
  NORM=splineAnalyze(Y=Z,map=Data$V4,smoothness = 2000,plotRaw=T,plotWindows = T,method = 4)
  dev.off()
  normScore=NORM$windowData
  outFile=paste("/data1/home/xuebo/Projects/Evo/xpclr/lineage/newfileV9/Agroup/XPCLRresult/normalozation_100k/normXPCLRscore_chr",i,".txt",sep="")
  write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
}

#GenWin的结果: "WindowStart" "WindowStop" "SNPcount" "MeanY" "Wstat"，没有染色体，要加上染色体的信息
#工作路径：xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/XPCLRresultSouth_SCA/smooth

#smooth结果添加染色体号并进行排序
#!/usr/bin/shell
for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
do
sed '1d' smooth$i.txt | awk '{print "'$i'""\t"$0}'  
done |sed '/NA/d' | sort -k5,5n -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > smooth_A.txt

for i in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
do
sed '1d' smooth$i.txt | awk '{print "'$i'""\t"$0}'  
done |sed '/NA/d' | sort -k5,5n -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > smooth_B.txt

for i in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
do
sed '1d' smooth$i.txt | awk '{print "'$i'""\t"$0}'  
done |sed '/NA/d' | sort -k5,5n -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > smooth_D.txt

#42条染色体合并成21条
bash bashA.sh 
#目标文件：/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/XPCLRresultSouth_SCA/smooth/smooth_A_change.txt
#画21条染色体的曼哈顿图

#每个亚基因组取5%看信号（42条染色体）

#XPCLRresultNorth_SCA
wc smooth_A.txt #99964*0.01=1000(3.76)
wc smooth_B.txt #125422*0.01=1254(3.98)
wc smooth_D.txt #107516*0.01=1075(3.61)

tail -n 1000 smooth_A.txt | sort -k1,1n -k2,2n | awk '{print $1"\t"$2"\t"$3}' | awk '{if($2<0) print $1"\t0\t"$3; else print $0}' > North_top1_A.txt
tail -n 1254 smooth_B.txt | sort -k1,1n -k2,2n | awk '{print $1"\t"$2"\t"$3}' | awk '{if($2<0) print $1"\t0\t"$3; else print $0}' > North_top1_B.txt
tail -n 1075 smooth_D.txt | sort -k1,1n -k2,2n | awk '{print $1"\t"$2"\t"$3}' | awk '{if($2<0) print $1"\t0\t"$3; else print $0}' > North_top1_D.txt

#定位基因
bedtools intersect -a gene_v1.1_Lulab.gff3 -b North_top1_A.txt -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > North_top1_A.gene
bedtools intersect -a gene_v1.1_Lulab.gff3 -b North_top1_B.txt -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > North_top1_B.gene
bedtools intersect -a gene_v1.1_Lulab.gff3 -b North_top1_D.txt -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > North_top1_D.gene

#功能富集
#B：GO:0009627	P	systemic acquired resistance	28	3262	1.5e-06	0.0079
#D：GO:0005319	F	lipid transporter activity	7	181	1.1e-05	0.0072

#XPCLRresultSouth_SCA
wc smooth_A.txt #94773*0.01=948(3.849)
wc smooth_B.txt #119459*0.01=1195(3.73)
wc smooth_D.txt #115796*0.01=1158(3.58)

tail -n 948  smooth_A.txt | sort -k1,1n -k2,2n | awk '{print $1"\t"$2"\t"$3}' | awk '{if($2<0) print $1"\t0\t"$3; else print $0}' > South_top1_A.txt
tail -n 1195 smooth_B.txt | sort -k1,1n -k2,2n | awk '{print $1"\t"$2"\t"$3}' | awk '{if($2<0) print $1"\t0\t"$3; else print $0}' > South_top1_B.txt
tail -n 1158 smooth_D.txt | sort -k1,1n -k2,2n | awk '{print $1"\t"$2"\t"$3}' | awk '{if($2<0) print $1"\t0\t"$3; else print $0}' > South_top1_D.txt

#定位基因
bedtools intersect -a gene_v1.1_Lulab.gff3 -b South_top1_A.txt -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > South_top1_A.gene
bedtools intersect -a gene_v1.1_Lulab.gff3 -b South_top1_B.txt -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > South_top1_B.gene
bedtools intersect -a gene_v1.1_Lulab.gff3 -b South_top1_D.txt -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > South_top1_D.gene

#定位基因
sed '1d' South_A_change.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n | head
sed '1d' South_A_change.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n > South_A_change.bed
sed '1d' South_B_change.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n > South_B_change.bed
sed '1d' South_D_change.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n > South_D_change.bed
cd ../../XPCLRresultNorth_SCA/smooth/
  sed '1d' North_A_change.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n > North_A_change.bed
sed '1d' North_B_change.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n > North_B_change.bed
sed '1d' North_D_change.txt | awk '{print $1"\t"$2"\t"$3}' | sort -k1,1n -k2,2n > North_D_change.bed
bedtools intersect -a pattern -b North_A_change.bed -wa 

#定位位点
#South_top5_A.gene
25	432228410	432232121	ID=TraesCS5A02G216000
#South_top5_D.gene
24	47379437	47393597	ID=TraesCS4D02G341700
41	12457171	12461806	ID=TraesCS7D02G026000
5	11451423	11459353	ID=TraesCS1D02G029100

#North_top5_A.gene
31	429395155	429398129	ID=TraesCS6A02G227900
7	81470715	81484373	ID=TraesCS2A02G135000  # 81479000	81511000	3	3.242731557295	5.61657581254189
7	80586966	80592975	ID=TraesCS2A02G134000
#North_top5_B.gene
34	253300670	253308145	ID=TraesCS6B02G440200
#North_top5_D.gene
24	47379437	47393597	ID=TraesCS4D02G341700
41	12457171	12461806	ID=TraesCS7D02G026000

#提取manhatan图中显著的基因位置
sed '1d' North_all_change.txt |awk '{print NR"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' | sed '1i SNP\tChr\tWindowStart\tSNPcount\tMeanY\tWstat' > North_manhatan.txt

#North
sed '1d' North_top5_A.txt | awk '{print $1"\t"$2"\t"$3}' > North_top5_A.bed
sed '1d' North_top5_B.txt | awk '{print $1"\t"$2"\t"$3}' > North_top5_B.bed
sed '1d' North_top5_D.txt | awk '{print $1"\t"$2"\t"$3}' > North_top5_D.bed
#构造gene.bed

bedtools intersect -a North_top5_A.bed -b gene.bed2 -wa | awk '{print $1"\t"$2+1000000"\t"$3-1000000}'> A.pos
grep -f A.pos smooth_A.txt| awk '{print $4"\t"$5}' > posA
grep -f posA North_manhatan.txt >> pos

bedtools intersect -a North_top5_B.bed -b gene.bed -wa | awk '{print $1"\t"$2+1000000"\t"$3-1000000}'> B.pos
grep -f B.pos smooth_B.txt| awk '{print $4"\t"$5}' > posB
grep -f posB North_manhatan.txt 

bedtools intersect -a North_top5_D.bed -b gene.bed -wa | awk '{print $1"\t"$2+1000000"\t"$3-1000000}'> D.pos
grep -f D.pos smooth_D.txt| awk '{print $4"\t"$5}' > posD
grep -f posD North_manhatan.txt 

#South
sed '1d' South_top5_A.txt | awk '{print $1"\t"$2"\t"$3}' > South_top5_A.bed
sed '1d' South_top5_B.txt | awk '{print $1"\t"$2"\t"$3}' > South_top5_B.bed
sed '1d' South_top5_D.txt | awk '{print $1"\t"$2"\t"$3}' > South_top5_D.bed
#构造gene.bed

bedtools intersect -a South_top5_A.bed -b gene.bed -wa | awk '{print $1"\t"$2+1000000"\t"$3-1000000}'> A.pos
grep -f A.pos smooth_A.txt| awk '{print $4"\t"$5}' > posA
grep -f posA South_manhatan.txt >> pos

bedtools intersect -a South_top5_B.bed -b gene.bed -wa | awk '{print $1"\t"$2+1000000"\t"$3-1000000}'> B.pos
grep -f B.pos smooth_B.txt| awk '{print $4"\t"$5}' > posB
grep -f posB South_manhatan.txt 

bedtools intersect -a South_top5_D.bed -b gene.bed -wa | awk '{print $1"\t"$2+1000000"\t"$3-1000000}'> D.pos
grep -f D.pos smooth_D.txt| awk '{print $4"\t"$5}' > posD
grep -f posD South_manhatan.txt 
