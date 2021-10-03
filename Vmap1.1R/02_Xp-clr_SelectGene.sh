#工作路径：xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V2/smooth
#CA_NW_A_smooth*
#CA_SA_smooth*
#NW_A_NE_A_smooth*
#SA_SW_A_smooth*
#SA_Tibet_smooth*
#SE_A_NE_A_smooth*
#SW_A_NE_A_smooth*
#Strang_WA_smooth*
#WA_CA_smooth*
#WA_EU_smooth*
#WA_NE_A_smooth*
#WA_NW_A_smooth*
#WA_SA_smooth*
#WA_SE_A_smooth*
#WA_SW_A_smooth*
#WA_Tibet_smooth*
#GenWin(smooth)的结果是这样的"WindowStart" "WindowStop" "SNPcount" "MeanY" "Wstat"，没有染色体，要加上染色体的信息
#smooth结果添加染色体号并进行排序，并且合并成A，B，D lineage.

Name=(CA_NW_A CA_SA NW_A_NE_A SA_SW_A SA_Tibet SE_A_NE_A SW_A_NE_A Strang_WA WA_CA WA_EU WA_NE_A WA_NW_A WA_SA WA_SE_A WA_SW_A WA_Tibet)
for file in Name
do
  for i in {1,2,7,8,13,14,19,20,25,26,31,32,37,38}
    do
      sed '1d' $file_smooth$i.txt | awk '{print "'$i'""\t"$0}'  
    done |sed '/NA/d' | sort -k5,5n -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > $file_smooth_A.txt
  for i in {3,4,9,10,15,16,21,22,27,28,33,34,39,40}
    do
      sed '1d' $file_smooth$i.txt | awk '{print "'$i'""\t"$0}'  
    done |sed '/NA/d' | sort -k5,5n -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > $file_smooth_B.txt
  for i in {5,6,11,12,17,18,23,24,29,30,35,36,41,42}
    do
      sed '1d' $file_smooth$i.txt | awk '{print "'$i'""\t"$0}'  
    done |sed '/NA/d' | sort -k5,5n -k1,1n -k2,2n | sed '1i Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > $file_smooth_D.txt
done

#42条染色体合并成21条(因为gff文件是42条染色体,此步可以跳过)
#bash 42-21chr.sh  
#目标文件：/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/XPCLRresultSouth_SCA/smooth/smooth_A_7.txt

#42条染色体smooth结果
#目标文件目录：/data2/xuebo/Projects/Speciation/xpclr/Selection_V2/smooth

#每个亚基因组取5%看信号(42条染色体): bash top.sh
#File_name	Line_num	Top1%_num	Top5%_num
EA_A_7.txt	98813	988	4941
EA_B_7.txt	105670	1057	5284
EA_D_7.txt	104488	1045	5224
EU_A_7.txt	116545	1165	5827
EU_B_7.txt	137135	1371	6857
EU_D_7.txt	123002	1230	6150
North_A_7.txt	99964	1000	4998
North_B_7.txt	125422	1254	6271
North_D_7.txt	107516	1075	5376
SCA_A_7.txt	101636	1016	5082
SCA_B_7.txt	132997	1330	6650
SCA_D_7.txt	137558	1376	6878
South_A_7.txt	94773	948	4739
South_B_7.txt	119459	1195	5973
South_D_7.txt	115796	1158	5790
WA_A_7.txt	139411	1394	6971
WA_B_7.txt	138839	1388	6942
WA_D_7.txt	73047	730	3652

#定位基因: bash gene.sh
bedtools intersect -a gene_v1.1_Lulab.gff3 -b North_top1_A.txt -wa | awk '{print $1"\t"$4"\t"$5"\t"$9}' | awk -F";" '{print $1}' | sort | uniq > North_top1_A.gene

#已知基因的列表

#定位位点

#画top5% nlr基因的曼哈顿图
204@xuebo:/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/smooth/Top5%/gene/Manhattan/gene
for i in `ls *A_gene.txt`; do grep -f ../../../nlr_gene.txt $i; echo $i; done

