#E6 VCF
204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6/VCF
204:xuebo:/data2/xuebo/Projects/Speciation/E6/Landrace_locate_225
#做迁徙和环境适应性路径
#42chr
204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6
203:yafei:/data1/home/yafei/008_Software/snpEff/data2
#lineages
204:yafei:/data2/yafei/003_Project3/Vmap1.1/Out_E6/VCF/VmapE6

#XP-clr + smooth
204:xuebo:/data2/xuebo/Projects/Speciation/xpclr/Selection_V2/smooth/smooth_result/Top5%
203:xuebo:/data1/home/xuebo/Projects/Speciation/xpclr/Selection_V2
#VIP gene分析
204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6/Xp-CLR_V2/VIP_genes
204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6/Xp-CLR_V2
#VCF 变异注释分析
203:yafei:/data1/home/yafei/008_Software/snpEff/

SNPeff----52个GA通路相关基因的变异分析
参考网站：https://www.jianshu.com/p/a6e46d0c07ee(SnpEff使用方法)
https://www.jianshu.com/p/f898ffc4ef48(SnpEff结果解读)
fasta文件: 203: /data1/home/yafei/Project3/HaploType/Fulab/Fa/
vcf文件: 203: /data1/home/yafei/Project3/HaploType/Fulab/Vcf/
软件位置: 203: /data1/home/yafei/Software/snpEff/
结果文件: 203: /data1/home/yafei/Software/snpEff/data/

#XP-CLR negative contral
/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/NegativeContral
#IBS with CS
203:/data2/yafei/003_project3/Project3/CS_Vmap
计算FST
203:/data2/yafei/003_project3/Project3/FST_group

Fst_population
203:/data2/yafei/003_project3/Project3/FST_group/VmapData/FST_selection/Fst_pop

#VMap3数据集
filter1: reliable library
filter2: hetThresh = 0.05; nonmissingThresh = 0.8; macThresh = 2; biallele
java -Xmx300g -jar /data2/yafei/004_Vmap3/Fastcall2/03_Jar/PrivatePlantGenetics.jar
204:/data2/yafei/004_Vmap3/Fastcall2/02_Output/vcf_AB1/Filter2/ 

#VMap3基本统计
204:/data2/yafei/004_Vmap3/Fastcall2/02_Output/vcf_AB1/Filter2/5k.tree/structure
