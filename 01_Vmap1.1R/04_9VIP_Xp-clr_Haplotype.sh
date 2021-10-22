#具体分析一些样本的单倍型，总共分成3类
#1. 保护祖先种的意义: 祖先种中存在的变异，在传播的过程中丢失。
  #eg(抗病): Sr33(TraesCS1D02G029100), Sr45(TraesCS1D02G040400), Tsn1(TraesCS5B02G059000)
#2. 保护地区特有品种的意义: 祖先种中没有，在传播的过程中适应区域内环境，形成新的变异。
  #eg(氮吸收，株高): NGR5_1(TraesCS1A02G242800), GID-A1(TraesCS1A02G255100), Rht-A1(TraesCS4A02G271000), Rht-B1(TraesCS4B02G043100)
#3. 单倍型在区域变化的例子。
  #eg(开花，光敏): TaFT8-2_1(TraesCS2A02G485000), TaPHYC_1(TraesCS5A02G391300)

#working directory
204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6/Xp-CLR_V2/VIP_genes
#提取祖先的对应位点vcf
#VIP_gene.txt(Old)
TraesCS1D02G029100
TraesCS1D02G040400
TraesCS5B02G059000
TraesCS1A02G242800
TraesCS1A02G255100
TraesCS4A02G271000
TraesCS4B02G043100
TraesCS2A02G485000
TraesCS5A02G391300
#VIP_gene.txt(New)
Ppd-A1
Ppd-1
VRN2-2

#提取带祖先的样本
grep -w -f VIP_gene.txt 87gene_5k.txt > VIPgene_5k.txt
bash getVcf.sh VIPgene_5k.txt > VIPgene_5k.log
#提取在六倍体里面分离的位点

for i in `cat VIP_gene.txt`; do sed '/#/d' *$i*.vcf | awk '{print $1"\t"$2}'; done > VIP_gene.pos

for i in `ls *vcf`; do vcftools --vcf $i --positions VIP_gene.pos --recode --out ${i::-11}.pos; done
#1.传到本地用Java Migration/Haplotype将vcf转换成单倍型txt格式，使用07_VCF_Haplotype_Visual.r进行可视化。
#本地路径：/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/VIP_gene/
#2.使用04_9VIP_Gene_Seq_anno.r进行变异注释。




