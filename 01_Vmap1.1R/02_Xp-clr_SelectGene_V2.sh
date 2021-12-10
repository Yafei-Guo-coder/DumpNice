#工作路径：xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%
#EU_South*
#North2_South*
#Strang_WA*
#Tibet_South*
#WA_EU*
#WA_South*
#neg_North1_North2*
#neg_WA_North1*
#neg_WA_North2*

#GenWin(smooth)的结果是这样的"WindowStart" "WindowStop" "SNPcount" "MeanY" "Wstat"，没有染色体，要加上染色体的信息
#smooth结果添加染色体号并进行排序，并且合并成A，B，D lineage.

Name=(EU_South North2_South Strang_WA Tibet_South WA_EU WA_South neg_North1_North2 neg_WA_North1 neg_WA_North2)
for num in {0..8}
do
  for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
    do
      sed '1d' ${Name[$num]}_smooth${i}.txt | awk '{print "'$i'""\t"$0}'  
    done |sed '/NA/d' | sort -k5,5g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > ${Name[$num]}_smooth_A.txt
  for i in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
    do
      sed '1d' ${Name[$num]}_smooth${i}.txt | awk '{print "'$i'""\t"$0}'  
    done |sed '/NA/d' | sort -k5,5g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > ${Name[$num]}_smooth_B.txt
  for i in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
    do
      sed '1d' ${Name[$num]}_smooth${i}.txt | awk '{print "'$i'""\t"$0}'  
    done |sed '/NA/d' | sort -k5,5g -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > ${Name[$num]}_smooth_D.txt
done

#42条染色体合并成21条(因为gff文件是42条染色体,此步可以跳过)
#bash 42-21chr.sh  
#目标文件：/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/XPCLRresultSouth_SCA/smooth/smooth_A_7.txt

#每个亚基因组取5%看信号(42条染色体): bash top.sh
#File_name	Line_num	Top1%_num	Top5%_num
for i in `ls *txt`; do  wc -l $i; done | awk '{print $2"\t"$1"\t"$1*0.01"\t"$1*0.05}' | awk -F"[.|\t]" '{print $1"\t"$3"\t"$4"\t"$6}' | awk '{print "tail -n "$3,$1".txt > Top1%/"$1".top1.bed"}'
#阈值3-4
for i in `ls *txt`; do  wc -l $i; done | awk '{print $2"\t"$1"\t"$1*0.01"\t"$1*0.05}' | awk -F"[.|\t]" '{print $1"\t"$3"\t"$4"\t"$6}' | awk '{print "tail -n "$4,$1".txt > Top5%/"$1".top5.bed"}'
#阈值1-2


#定位已注释基因:gff
for i in `ls *bed`
do
bedtools intersect -a ../gene_v1.1_Lulab.gff3 -b $i -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > ${i::-3}gff.gene
done

#定位已克隆基因:cloned gene 以及定位抗病基因:nlr gene
for i in `ls *gff.gene`
do
awk -F"ID=" '{print $2}' $i > ${i::-8}gene2
done
for i in `ls *.gene2`
do
grep -w -f $i ../all_cloned_gene.txt > ${i::-5}cloned.gene
#grep -w -f $i ../nlr_gene.txt > nlr/${i::-5}cloned.gene
done
rm *gene2

ls *cloned.gene |xargs -n1 > cloned_gene.txt
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done | awk '{print $1}' | awk -F"_smooth" '{print $1}'|uniq > file_prefix.txt
#change file format to plot heatmap(A,B,D lineage seperate)

#A lineage
ls *A.top5.cloned.gene |xargs -n1 > A_cloned_gene.txt
sed 's/$/_smooth_A.top5.cloned.gene/' file_prefix.txt > A_file.txt
for i in `cat A_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > A_gene.txt

#B lineage
ls *B.top5.cloned.gene |xargs -n1 > B_cloned_gene.txt
sed 's/$/_smooth_B.top5.cloned.gene/' file_prefix.txt> B_file.txt
for i in `cat B_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > B_gene.txt

#D lineage 
ls *D.top5.cloned.gene |xargs -n1 > D_cloned_gene.txt
sed 's/$/_smooth_D.top5.cloned.gene/' file_prefix.txt > D_file.txt
for i in `cat D_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > D_gene.txt

#R
v <- c("A","B","D")
for ( i in c(1:3)){
  filename <- paste(v[i],"_file.txt",sep="")
  genename <- paste(v[i],"_gene.txt",sep="")
  file <- read.table(filename,header=F,stringsAsFactors=F)
  gene <- read.table(genename,header=F,stringsAsFactors=F)
  files <- file[,1]
  genes <- gene[,1]
  x=vector()
  for (j in files){
    a<-paste(j,genes,sep=" ")
    x <- c(x,a)
  }
  out <- paste(v[i],"_file_gene_mode.txt",sep="")
  write.table(x,out,quote=F,row.names=F,col.names=F,sep="\t")
}

#shell
#A lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep A.top5.cloned.gene > A_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' A_postive_file_gene_mode.txt A_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > A_heatmap_format1.txt
#B lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep B.top5.cloned.gene > B_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' B_postive_file_gene_mode.txt B_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > B_heatmap_format1.txt
#D lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep D.top5.cloned.gene > D_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' D_postive_file_gene_mode.txt D_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > D_heatmap_format1.txt
#R
library(reshape)
name <- read.table("file_prefix.txt",header=F,stringsAsFactors=F)
colnames(name) <- "file"
v <- c("A","B","D")
num <- c(1,0.6,0.2)
for ( i in c(1:3)){
  filename <- paste(v[i],"_heatmap_format1.txt",sep="")
  file <- read.table(filename,header=F,stringsAsFactors=F)
  sub <- file[,c(1,3,4)]
  colnames(sub) <- c("file","gene","type")
  sub[which(sub$type==1),3] <- num[i]
  sub[which(sub$type==0),3] <- num[i]-0.2
  mt <- melt(sub,id=c("file","gene"))
  st <- cast(mt,file~variable+gene)
  name<-merge(name,st,by="file") 
}
write.table(name,"heatmap_format2.txt",quote=F,row.names=F,col.names=T,sep="\t")

#Go_analysis(invalid:应该用全部的信号基因)
for i in `ls *gff.gene`
do
awk -F"=" '{print $2}' $i |sort | uniq > Go/${i::-8}txt
done

#新加了wheatOmic的小麦中的已知基因
#工作路径：204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Gene_Anno
#更新了../all_cloned_gene.txt的已知基因
#../GeneLulab1_1_LC.gff3:low confidence的基因注释文件
#../gene_v1.1_Lulab.gff3:lhigh confidence的基因注释文件
#known.gene.txt:517个已知基因（BGC1-4A	TraesCS4A02G284000)
#定位已克隆基因:cloned gene 以及定位抗病基因:nlr gene
for i in `ls *gff.gene`
do
awk -F"ID=" '{print $2}' $i > ${i::-8}gene2
done
for i in `ls *.gene2`
do
grep -w -f $i known.gene.txt > ${i::-5}cloned.gene
done
rm *gene2

ls *cloned.gene |xargs -n1 > cloned_gene.txt
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done | awk '{print $1}' | awk -F"_smooth" '{print $1}'|uniq > file_prefix.txt
#change file format to plot heatmap(A,B,D lineage seperate)

#A lineage
ls *A.top5.cloned.gene |xargs -n1 > A_cloned_gene.txt
sed 's/$/_smooth_A.top5.cloned.gene/' file_prefix.txt > A_file.txt
for i in `cat A_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > A_gene.txt

#B lineage
ls *B.top5.cloned.gene |xargs -n1 > B_cloned_gene.txt
sed 's/$/_smooth_B.top5.cloned.gene/' file_prefix.txt> B_file.txt
for i in `cat B_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > B_gene.txt

#D lineage 
ls *D.top5.cloned.gene |xargs -n1 > D_cloned_gene.txt
sed 's/$/_smooth_D.top5.cloned.gene/' file_prefix.txt > D_file.txt
for i in `cat D_cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $2}' | sort | uniq > D_gene.txt

#R
v <- c("A","B","D")
for ( i in c(1:3)){
  filename <- paste(v[i],"_file.txt",sep="")
  genename <- paste(v[i],"_gene.txt",sep="")
  file <- read.table(filename,header=F,stringsAsFactors=F)
  gene <- read.table(genename,header=F,stringsAsFactors=F)
  files <- file[,1]
  genes <- gene[,1]
  x=vector()
  for (j in files){
    a<-paste(j,genes,sep=" ")
    x <- c(x,a)
  }
  out <- paste(v[i],"_file_gene_mode.txt",sep="")
  write.table(x,out,quote=F,row.names=F,col.names=F,sep="\t")
}
#shell
#A lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep A.top5.cloned.gene > A_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' A_postive_file_gene_mode.txt A_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > A_heatmap_format1.txt
#B lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep B.top5.cloned.gene > B_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' B_postive_file_gene_mode.txt B_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > B_heatmap_format1.txt
#D lineage
for i in `cat cloned_gene.txt`; do awk '{print "'$i'""\t"$0}' $i; done| awk '{print $1" "$2}' | grep D.top5.cloned.gene > D_postive_file_gene_mode.txt
awk 'NR==FNR{a[$0]=$0}NR!=FNR{if($0 in a) {print a[$0]"\t1"} else print $0"\t"0}' D_postive_file_gene_mode.txt D_file_gene_mode.txt | awk -F"_smooth_" '{print $1"\t"$2}'|sed 's/ /\t/' > D_heatmap_format1.txt
#R
library(reshape)
name <- read.table("file_prefix.txt",header=F,stringsAsFactors=F)
colnames(name) <- "file"
v <- c("A","B","D")
num <- c(1,0.6,0.2)
for ( i in c(1:3)){
  filename <- paste(v[i],"_heatmap_format1.txt",sep="")
  file <- read.table(filename,header=F,stringsAsFactors=F)
  sub <- file[,c(1,3,4)]
  colnames(sub) <- c("file","gene","type")
  sub[which(sub$type==1),3] <- num[i]
  sub[which(sub$type==0),3] <- num[i]-0.2
  mt <- melt(sub,id=c("file","gene"))
  st <- cast(mt,file~variable+gene)
  name<-merge(name,st,by="file") 
}
write.table(name,"heatmap_format2.txt",quote=F,row.names=F,col.names=T,sep="\t")







