204:/data1/home/yafei/003_Project3/bayenv/13pop/permulation/Top5
204:/data1/home/yafei/003_Project3/bayenv/13pop/XP-CLR/Top5/
Folder=(CA_EA CA_SA EA_SA EU_CA EU_EA EU_SA Tibet_SA)
for num in {0,1,2,3,5,6}
do
{
	mkdir ${Folder[$num]}/All ${Folder[$num]}/Prec ${Folder[$num]}/Temp
	cp /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2_Top5/All/*log.bayenv.top5.bed ${Folder[$num]}/All
	cp /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2_Top5/Prec/*prec.bayenv.top5.bed ${Folder[$num]}/Prec
	cp /data1/home/yafei/003_Project3/bayenv/13pop/Add544Gene_V2_Top5/Temp/*temp.bayenv.top5.bed ${Folder[$num]}/Temp
	sed  "s/Folder\[1\]/${Folder[$num]}/" All_gene.sh > ${Folder[$num]}_All_gene.sh
	sed  "s/Folder\[1\]/${Folder[$num]}/" Prec_gene.sh > ${Folder[$num]}_Prec_gene.sh
	sed  "s/Folder\[1\]/${Folder[$num]}/" Temp_gene.sh > ${Folder[$num]}_Temp_gene.sh
	bash ${Folder[$num]}_All_gene.sh
	bash ${Folder[$num]}_Prec_gene.sh
	bash ${Folder[$num]}_Temp_gene.sh
	}
done

#All_gene.sh
#路径
#204:/data1/home/yafei/003_Project3/bayenv/13pop/XP-CLR/Top5
#permulation: 
#203:/data1/home/yafei/003_Project3/XP-CLR/permutation
#204:/data1/home/yafei/003_Project3/bayenv/13pop/permulation
cd All
mkdir bayenv_xpclr_site
mkdir bayenv_xpclr_region
cp ../*.bed ./
#all
for i in `ls *smooth_A.top5.bed`
do
bedtools intersect -b $i -a Alog.bayenv.top5.bed  -wa |sort | uniq > bayenv_xpclr_site/$i
done
for i in `ls *smooth_B.top5.bed`
do
bedtools intersect -b $i -a Blog.bayenv.top5.bed -wa|sort | uniq > bayenv_xpclr_site/$i
done
for i in `ls *smooth_D.top5.bed`
do
bedtools intersect -b $i -a Dlog.bayenv.top5.bed -wa |sort | uniq > bayenv_xpclr_site/$i
done

#all
for i in `ls *smooth_A.top5.bed`
do
bedtools intersect -a $i -b Alog.bayenv.top5.bed  -wa |sort | uniq > bayenv_xpclr_region/$i
done
for i in `ls *smooth_B.top5.bed`
do
bedtools intersect -a $i -b Blog.bayenv.top5.bed -wa|sort | uniq > bayenv_xpclr_region/$i
done
for i in `ls *smooth_D.top5.bed`
do
bedtools intersect -a $i -b Dlog.bayenv.top5.bed -wa |sort | uniq > bayenv_xpclr_region/$i
done

#--------------------------------------site------------------------------
cd /data1/home/yafei/003_Project3/bayenv/13pop/XP-CLR/Top5/All/bayenv_xpclr_site
mkdir nlr
for i in `ls *bed`
do
bedtools intersect -a /data1/home/yafei/003_Project3/bayenv/13pop/gene_v1.1_Lulab.gff3 -b $i -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > ${i::-3}gff.gene
done

#定位已克隆基因:cloned gene 以及定位抗病基因:nlr gene
for i in `ls *gff.gene`
do
awk -F"ID=" '{print $2}' $i > ${i::-8}gene2
done
for i in `ls *.gene2`
do
grep -w -f $i /data1/home/yafei/003_Project3/bayenv/13pop/cloned.gene.txt > ${i::-5}cloned.gene
grep -w -f $i /data1/home/yafei/003_Project3/bayenv/13pop/nlr_gene.txt > nlr/${i::-5}cloned.gene
done
rm *gene2
#---------------------------------------clone_gene------------------
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
Rscript /data1/home/yafei/003_Project3/bayenv/13pop/gene_mode.r
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
#R
Rscript /data1/home/yafei/003_Project3/bayenv/13pop/heatmap_format2.r
#---------------------------------------nlr------------------
cd nlr
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
Rscript /data1/home/yafei/003_Project3/bayenv/13pop/gene_mode.r
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
#R
Rscript /data1/home/yafei/003_Project3/bayenv/13pop/heatmap_format2.r

#--------------------------------------region------------------------------
cd /data1/home/yafei/003_Project3/bayenv/13pop/XP-CLR/Top5/All/bayenv_xpclr_region
mkdir nlr
for i in `ls *bed`
do
bedtools intersect -a /data1/home/yafei/003_Project3/bayenv/13pop/gene_v1.1_Lulab.gff3 -b $i -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > ${i::-3}gff.gene
done

#定位已克隆基因:cloned gene 以及定位抗病基因:nlr gene
for i in `ls *gff.gene`
do
awk -F"ID=" '{print $2}' $i > ${i::-8}gene2
done
for i in `ls *.gene2`
do
grep -w -f $i /data1/home/yafei/003_Project3/bayenv/13pop/cloned.gene.txt > ${i::-5}cloned.gene
grep -w -f $i /data1/home/yafei/003_Project3/bayenv/13pop/nlr_gene.txt > nlr/${i::-5}cloned.gene
done
rm *gene2
#---------------------------------------clone_gene------------------
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
Rscript /data1/home/yafei/003_Project3/bayenv/13pop/gene_mode.r
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
#R
Rscript /data1/home/yafei/003_Project3/bayenv/13pop/heatmap_format2.r
#---------------------------------------nlr------------------
cd nlr
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
Rscript /data1/home/yafei/003_Project3/bayenv/13pop/gene_mode.r
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

Rscript /data1/home/yafei/003_Project3/bayenv/13pop/heatmap_format2.r



cd /data1/home/yafei/003_Project3/bayenv/13pop/XP-CLR/Top5/All

for i in `ls *smooth_A.top5.bed`; do wc -l Alog.bayenv.top5.bed; awk 'BEGIN{count=0} {count=count+$3-$2} END{print count}' $i; echo $i;echo "4934891648"; done | xargs -n5 | grep -v total >> A
#wc -l Alog.bayenv.top5.bed >>A
cd bayenv_xpclr_site/
for i in `ls *smooth_A.top5.bed`; do wc -l $i;done |awk '{print $1}' |grep -v total >> A

cd ../
for i in `ls *smooth_B.top5.bed`; do wc -l Blog.bayenv.top5.bed ;awk 'BEGIN{count=0} {count=count+$3-$2} END{print count}' $i; echo $i;echo "5180314468"; done | xargs -n5 | grep -v total >> B
#wc -l Blog.bayenv.top5.bed >>B
cd bayenv_xpclr_site/
for i in `ls *smooth_B.top5.bed`; do wc -l $i;done|awk '{print $1}' |grep -v total >> B

cd ../
for i in `ls *smooth_D.top5.bed`; do wc -l Dlog.bayenv.top5.bed; awk 'BEGIN{count=0} {count=count+$3-$2} END{print count}' $i; echo $i;echo "3951074735"; done | xargs -n5 | grep -v total >> D
#wc -l Dlog.bayenv.top5.bed >>D
cd bayenv_xpclr_site/
for i in `ls *smooth_D.top5.bed`; do wc -l $i;done |awk '{print $1}'|grep -v total >> D

paste ../A A -d " " >> out
paste ../B B -d " " >> out
paste ../D D -d " " >> out
cat out | awk '{print $4"\t" $1/$5"\t"$6/$3"\t"$6/$3/($1/$5)}' > /data1/home/yafei/003_Project3/bayenv/13pop/XP-CLR/Top5/All/out.lineage.txt
cat out | sed 's/_smooth_A.top5.bed//' | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6}' |sed 's/_smooth_B.top5.bed//'|sed 's/_smooth_D.top5.bed//' |sort -k3| datamash -g 3 sum 1 sum 2 sum 4 sum 5 | awk '{print $1"\t"$2/$4"\t"$5/$3"\t"$5/$3/($2/$4)}' > /data1/home/yafei/003_Project3/bayenv/13pop/XP-CLR/Top5/All/out.all.txt
rm A B D
cd ../
rm A B D

#heatmap
library(pheatmap)
library(RColorBrewer)
display.brewer.all()
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/10_Gene")
##Top5 VIP gene & Top1 NLR gene
data2 <- read.table("heatmap_format2.txt",
                    header=T,
                    row.names= 1, stringsAsFactors=F,sep="\t")
data <- as.matrix(data2)
x=vector()
for ( i in c(1:dim(data)[2])){
  x <- c(x,strsplit(colnames(data), "type_")[[i]][2])
}
colnames(data)<-x
#region <- c("Strang_WA", "WA_EU", "WA_CA", "WA_NW_A", "WA_NE_A", "WA_SA", "WA_SW_A", "WA_Tibet", "WA_SE_A", "CA_SA", "SA_SW_A", "SA_Tibet", "NW_A_NE_A", "SW_A_NE_A", "SE_A_NE_A")
#region <- c("neg_North1_North2", "neg_WA_North1", "neg_WA_North2", "WA_EU", "Strang_WA","WA_South","Tibet_South", "North2_South")
region <- c("neg_WA_1_2","Tibet_South","EU_North1","EU_North2","EU_South","neg_North1_North2","neg_WA_North1","neg_WA_North2","North2_South","South_North1","WA_EU","WA_South")
data <- data[region,]
colsA = c("#F7FCFD","#8C96C6","#F7FCF5","#74C476","#FFFFCC","#FD8D3C")
pheatmap(data,cluster_rows = F,cluster_cols = T,border_color=NA,color = colsA,fontsize=8)
