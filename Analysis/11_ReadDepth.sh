#-------------------------------------------------------------------计算reads depth.sh--------------------------------------------------------------------------------------
#工作目录：204:/data2/xuebo/Projects/Speciation/More_accessions/hapScan/
#位点depth
samtools depth -b bed_file sample.bam > sample.depth  #bed 用来指定统计区间，运行后输出指定区间每一个碱基的测序深度。
#区间depth
bedtools makewindows -g genome.txt -w 10000000 -s 2000000 > windows.bed
#bedtools makewindows用来自动生成划窗区间。-g genome.txt是要划分的基因组，格式为两列：染色体、染色体长度；-w 10000000为窗口大小为10M；-s 1000000为步长为1M，即窗口在染色体上每次向右平移1M的距离；windows.bed为输出文件，格式为三列：染色体、区间开始位点、区间结束位点。
bedtools coverage -a windows.bed -b xxx.sort.bam > xxx.depth.txt      
#bedtools coverage对划分好的每个滑动窗口进行reads数（depth)的统计。-a windows为上一步划分好的区间；-b xxx.sort.bam为测序数据mapping到参考基因组的比对文件；xxx.depth.txt为统计结果的输出文件，格式为7列：染色体、区间起始位点、区间结束位点、该区间内的reads数、该区间内的碱基数、区间大小、该区间的平均覆盖度。       
samtools bedcov bed_file samplename.bam > sample.bedcov
#输出的文件中计算了bed文件每一个区间的碱基总数，这里并不是reads的条数
