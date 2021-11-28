#---------------------------Dlineage_PopDiv_no.est---------------------------
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max
//all Ns are in number of haploid individuals
1  pop1		logunif	100	 100000	 output
1  pop2		logunif	100	 100000	 output
0  RSANC	logunif	0.1	 100	 output
0  TPROP	logunif 0.0001	 0.5	 output

[RULES]

[COMPLEX PARAMETERS]
#---------------------------------- Dlineage_PopDiv_no.tpl -------------------------------
//Parameters for the coalescence simulation program : fastsimcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
pop1
pop2
//Haploid samples sizes 
40
40
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2 historical event
10000 1 0 1 RSANC 0 0
6500 0 1 0 RSANC 0.0002 0
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 6.5e-9 OUTEXP

#run
./fsc26 -t Dlineage_PopDiv_diff.tpl -e Dlineage_PopDiv_diff.est -d -0 -n 100000 -L 40 -s 0 -M -q -c80
-d:dsfs
-0:不考虑 SFS 似然计算的单态位点
-C:在似然计算中考虑的最小观察到的 SFS 条目计数（默认值 = 1，但值可以 < 1。例如 0.5）
-n:要执行的模拟次数，也应用于参数估计
-L:在 lhood 最大化期间要执行的循环（ECM 循环）数。 默认为 20
-s:输出 DNA 作为 SNP 数据，并指定最大数量。 要输出的 SNP（使用 0 输出所有 SNP）。
-M:从迭代之间的 SFS 值通过 maxlhood 执行参数估计
-q:quiet
-c:用于参数估计的openMP线程数（默认=1，max=numBatches，使用0让openMP选择最优值）
  
-i:参数文件的名字
-E:从参数先验中抽取的次数（可选）, 在模板文件中替换列出的参数值
-g:使用基因型数据生成 arlequin 项目
-S:输出整个 DNA 序列，包括单态位点
-I:根据无限位点 (IS) 突变模型生成 DNA 突变
-m:msfs

Historical events can be used to: 
• Change the size of a given population. 
• Change the growth rate of a given population. 
• Change the migration matrix to be used between populations. 
• Move a fraction of the genes of a given population to another population. This amounts to 
implementing a (stochastic) admixture or introgression event. 
• Move all genes from a population to another population. This amounts to fusing two 
populations into one, looking backward in time, or implementing a populations’ fission 
looking forward in time. 
• One or more of these events at the same time

