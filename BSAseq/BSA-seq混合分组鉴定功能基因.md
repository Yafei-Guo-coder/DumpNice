## BSA-seq混合分组鉴定功能基因
### 1.使用QTLseq软件by YuSugihara
#### （1）使用conda安装QTLseq。
```bash
$ conda install -c bioconda qtlseq
```
#### （2）运行qtlseq命令
qtlseq能够指定多种输入文件，包括fastaq，bam和vcf文件。具体命令和软件包信息可以参考：https://github.com/YuSugihara/QTL-seq

此处我们以bam文件为输入文件，运行qtlseq软件：
```bash
$ qtlseq -r TAIR10.fa -p 8238.rmdup.bam -b1 HR.rmdup.bam -b2 NHR.rmdup.bam -n1 25 -n2 25 -w 1000 -s 10 -o BSA_result
```
`-p`指定亲本的bam文件；
`-b1`指定混池1的bam文件；
`-b2`指定混池2的bam文件；


计算结果会分成五个文件夹，分别存放于`BSA_result`文件夹下，文件夹的结构如下：
```
├── 10_ref
│   ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -> /mnt/zhou/hangyuan/BSA-seq/P20210613_AtBSA-raw.bam.results_20210621/03.bamData/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
│   ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.amb
│   ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.ann
│   ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.bwt
│   ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai
│   ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.pac
│   └── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.sa
├── 20_bam
│   ├── bulk.bam
│   ├── bulk.bam.bai
│   ├── cultivar.bam
│   └── cultivar.bam.bai
├── 30_vcf
│   ├── mutmap.vcf.gz
│   └── mutmap.vcf.gz.tbi
├── 40_mutmap
│   ├── mutmap_plot.png
│   ├── sliding_window.p95.tsv
│   ├── sliding_window.p99.tsv
│   ├── sliding_window.tsv
│   ├── snp_index.p95.tsv
│   ├── snp_index.p99.tsv
│   └── snp_index.tsv
└── log
    ├── bcftools.log
    ├── bwa.log
    ├── samtools.log
    └── tabix.log
```
文件夹中最重要的信息保存于40_mutmap(qtlseq)，包含了软件自带绘图包画出来的图，但是比较丑。超过p95和p99的SNP被保存于`snp_index.p95.tsv`和`snp_index.p99.tsv`文件，可以用作区间来定位候选基因。

#### （2）运行qtlplot命令，使用vcf文件单独画图。
当然如果已经有了vcf文件，也可以直接运行qtlplot计算vcf文件中各个snp位点的snp-index。此处`-w`指的是window size，单位是kb，1000kb代表1M。`-s`代表步长，单位也是kb，10kb代表以10k为步长滑窗。输出结果存放于新建的`newindow`文件夹中。
```bash
$ cd 30_vcf
$ qtlplot -v qtlseq.vcf.gz -o newindow -n1 25 -n2 25 -w 1000 -s 10
```
#### （3）snp-index曼哈顿图美化
QTLseq软件包绘制的散点图可以单独美化。方式如下：
下载`snp_index.tsv`和`sliding_window.tsv`，将两个文件改名为txt文件。输入到R脚本`BSAseq.R`中，前者用于绘制结果图中的散点，后者用于绘制平均线和阈值线。R脚本需要根据输入文件的染色体数修改部分参数，例如配色。之后只需要一行一行跑完R脚本即可。

-----

### 2.使用R语言包QTLseqr
相较于QTLseq软件方法，本方法需要手动进行call snp步骤，将得到的vcf文件使用gatk特定命令转化为table文件，并作为输入文件导入r包。过程略微复杂，能够得到线性绘图结果，但是画对于画散点图并不友好。

#### （1）参考基因组建立索引
使用GATK之前，需要先建立参考基因组索引文件.dict和.fai。`.dict`中包含了基因组中contigs的名字，也就是一个字典；`.fai`也就是fasta index file，索引文件，可以快速找出参考基因组的碱基，由samtools faidx构建。因此在存放有参考基因组的文件夹中，使用shell脚本`gatk_step1_index.sh`完成：
```bash
$ sh gatk_step1_index TAIR10.fasta

$ mv TAIR10.fa.dict TAIR10.dict #字典文件改名
```
注意！运行完脚本之后应当**手动将**`fa.dict`**文件改名为**`.dict`，否则后继call snp步骤会报错。

#### （2）原始测序数据质控，bwa比对，samtools排序及GATK call snp
直接运行脚本`gatk_step2_callvar.sh`，可自动执行全部步骤
```bash
$ sh gatk_step2_callvar.sh sample_1.fastq sample_2.fastq out_dictionary
```
计算结果会分成四个文件夹，分别存放于`out_dictionary`文件夹下，命名为bwa，cleanfq，gatk。文件夹的结构如下：
```
.
└── 8238
    ├── bwa
    │   ├── 8238.bam
    │   ├── 8238.log.mark
    │   ├── 8238.markup_metrics.txt
    │   ├── 8238.sorted.bam
    │   ├── 8238.sorted.bam.bai
    │   ├── 8238.sorted.markup.bam
    │   └── 8238.sorted.markup.bam.bai
    ├── cleanfq
    │   ├── 8238.html
    │   ├── 8238.json
    │   ├── 8238_1.fastq
    │   └── 8238_2.fastq
    └── gatk
        ├── 8238.HC.gvcf.gz
        └── 8238.HC.gvcf.gz.tbi 
```

#### （3）GATK合并GVCF文件
手动移动上步得到的gvcf文件及其对应的tbi文件到新的文件夹，执行`gatk_step3_gvcfmerg.sh`脚本。混池测序需要将两个极端性状池call出的vcf文件放到一起。
```bash
$ sh gatk_step3_gvcfmerg.sh
```
**注意！此脚本使用前需要提前手动修改脚本中待合并的文件！**

#### （4）SNP过滤
将脚本`gatk_step4_filter_SNP.sh`存放入合并后的vcf文件所在的文件夹中。运行此脚本以过滤SNP，**过滤参数如需自定义需要提前修改脚本**。

脚本4内置了gatk VariantsToTable命令，可以将vcf文件转化为r语言QTLseqr包所需的输入文件。
```bash
$ sh gatk_step4_filter_SNP.sh
```
最终我们得到了文件`BSA.filter.table`，输入进R语言完成后即步骤。后继使用r脚本`QTLseq_BSA.R`完成，输入进去一行行往下跑就能出图。