2018.11.5 周一
#常用命令总结
ll -rht
df -h
free -g 
history 对输出的命令进行展示
ps -aux|grep cp |grep mnt|awk '{print $2}'|xargs kill -9
ps -aux|grep bwa|awk '{print $2}'|xargs kill -9
ps -aux|grep vcftools|awk '{print $2}'|xargs kill -9
ps -aux|grep yuhong|awk '{print $2}'|xargs kill -9

kill -9 269378  # KILL    9    强制终止
iostat -k 1  #查看传输速度
blkid #查看文件存储位置
#查看命令运行开始的时间，2步操作
export HISTTIMEFORMAT='%F %T ' 
history | grep 'bwa'
ls -l | grep '^' | wc -l



#1 先解决了硬盘拷贝问题，用exfat格式； 再解决了数据上传问题，高老师帮助。
命令杀彻底,如下：  
ps -aux|grep cp |grep mnt|awk '{print $2}'|xargs kill -9
nohup cp -Rf /mnt/usb/ABD001/ ./ &  #可以用jobs查看
nohup cp -Rf /mnt/usb/ABD002/ ./ &
chown -R sharedData:sharedData Jiao/ #把拷贝的数据权限改成shareData：

#2 比对不要超过11个CPU，先计划一下怎么跑？

#路径大全
/data2/sharedData/Jiao/ABD001
/data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz  #索引库文件目录 
/data2/aoyue/output/bamfile/  #输出文件目录

#第一批数据：
#原始数据存放位置
aoyue@lulab1:/data2/sharedData/Jiao/ABD
#*fixmate.pos.bam文件存放位置
aoyue@lulab1:/data2/aoyue/output/bamfile
#*rmdup.bam文件以及rmdup.bam.bai文件存放位置
aoyue@lulab1:/data2/aoyue/output/bamfile

#第二批数据：
#原始数据存放位置
aoyue@lulab1:/data2/sharedData/Jiao/ABD2
aoyue@lulab2:/data2/aoyue/ABD2
#*fixmate.pos.bam文件存放位置
aoyue@lulab2:/data2/aoyue/output/bamfile
#*rmdup.bam文件以及rmdup.bam.bai文件存放位置
aoyue@lulab2:/data2/aoyue/output/bamfile

####################################  第一批数据的测试及运行  ################################
####################################  第一批数据的测试及运行  ################################
####################################  第一批数据的测试及运行  ################################

####################################  BWA  ##############################################
####################################  BWA  ##############################################
####################################  BWA  ##############################################

###  测试BWA  ###
###  测试BWA  ###
###  测试BWA  ###
#先做一个文件测试：
bwa mem -t 30 -R '@RG\tID:GRC-L1\tPL:illumina\tSM:GRC-L1\tLB:GRC-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD001/GRC-L1_1.fq.gz /data2/sharedData/Jiao/ABD001/GRC-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/GRC-L1.pe.bam && echo "** bwa mapping done **" &
bwa mem -t 30 -R '@RG\tID:HM\tPL:illumina\tSM:HM\tLB:HM' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD001/HM_1.fq.gz /data2/sharedData/Jiao/ABD001/HM_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HM.pe.bam && echo "** bwa mapping done **" &
日志信息：   
[aoyue@lulab1 bamfile]$ bwa mem -t 30 -R '@RG\tID:GRC-L1\tPL:illumina\tSM:GRC-L1\tLB:GRC-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD001/GRC-L1_1.fq.gz /data2/sharedData/Jiao/ABD001/GRC-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/GRC-L1.pe.bam && echo "** bwa mapping done **" &
[1] 274921
[aoyue@lulab1 bamfile]$ bwa mem -t 30 -R '@RG\tID:HM\tPL:illumina\tSM:HM\tLB:HM' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD001/HM_1.fq.gz /data2/sharedData/Jiao/ABD001/HM_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HM.pe.bam && echo "** bwa mapping done **" &
[2] 274924
#再做几个文件测试：
bwa mem -t 30 -R '@RG\tID:HNSH\tPL:illumina\tSM:HNSH\tLB:HNSH' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD001/HNSH_1.fq.gz /data2/sharedData/Jiao/ABD001/HNSH_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HNSH.pe.bam && echo "** bwa mapping done **" &
[1] 331425
bwa mem -t 100 -R '@RG\tID:HRV-L1\tPL:illumina\tSM:HRV-L1\tLB:HRV-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD/HRV-L1_1.fq.gz /data2/sharedData/Jiao/ABD/HRV-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HM.pe.bam && echo "** bwa mapping done **" &
bwa mem -t 100 -R '@RG\tID:HRV-L1\tPL:illumina\tSM:HRV-L1\tLB:HRV-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD/HRV-L1_1.fq.gz /data2/sharedData/Jiao/ABD/HRV-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HRV-L1.pe.bam && echo '** bwa mapping done **' 
[1] 357936
#出了问题
[M::mem_process_seqs] Processed 2000000 reads in 3276.598 CPU sec, 110.731 real sec
[M::mem_process_seqs] Processed 2000000 reads in 3363.172 CPU sec, 112.584 real sec
bash: 行 9: 274922 断开的管道         bwa mem -t 30 -R '@RG\tID:GRC-L1\tPL:illumina\tSM:GRC-L1\tLB:GRC-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD001/GRC-L1_1.fq.gz /data2/sharedData/Jiao/ABD001/GRC-L1_2.fq.gz
     274923 已终止               | samtools view -S -b - > /data2/aoyue/output/bamfile/GRC-L1.pe.bam
[M::process] read 2000000 sequences (300000000 bp)...
#HM这个样已完成
Version: 0.7.17-r1188
[main] CMD: bwa mem -t 30 -R @RG\tID:HM\tPL:illumina\tSM:HM\tLB:HM /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD001/HM_1.fq.gz /data2/sharedData/Jiao/ABD001/HM_2.fq.gz
[main] Real time: 57750.036 sec; CPU: 1706862.997 sec
** bwa mapping done **
[2]+  完成                  bwa mem -t 30 -R '@RG\tID:HM\tPL:illumina\tSM:HM\tLB:HM' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD001/HM_1.fq.gz /data2/sharedData/Jiao/ABD001/HM_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HM.pe.bam && echo "** bwa mapping done **"

样品	存储大小	 CPU数	real time(sec)	CPUtime 实际CPU
HM_1.fq.gz 38G 30	57750.036(16h) 474.129 hour	29.5

###  正式运行命令BWA  ###
###  正式运行命令BWA  ###
###  正式运行命令BWA  ###

#1 建立列表
ls *_1.fq.gz |sed 's/_1.fq.gz//g' > SM.txt 
#2 生成命令在一个脚本里
cat SM.txt | while read a;do echo "bwa mem -t 100 -R '@RG\tID:${a}\tPL:illumina\tSM:${a}\tLB:${a}' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD/${a}_1.fq.gz /data2/sharedData/Jiao/ABD/${a}_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** bwa mapping done **' " ;done > bwa.sh
#3 生成命令 一个样本一个脚本(只是测试用)
cat /Users/Aoyue/Documents/SM.txt | while read a;do echo "bwa mem -t 20 a.fq.gz ${a}_1.fq.gz ${a}_2.fq.gz && echo 'Well done'" >${a}.bwa.sh;done

##20181113 任务停掉   重新开始运行
[aoyue@lulab1:/data2/aoyue]sh bwa.sh &
[1] 16023

##20181116 任务停掉   重新开始运行
[aoyue@lulab1:/data2/aoyue]sh bwa.sh &
[1] 223638

########################### 排序、去重复和建立索引 ##########################################
########################### 排序、去重复和建立索引 ##########################################
########################### 排序、去重复和建立索引 ##########################################

### 20181122 周四 测试CPU时间 ###
### 20181122 周四 测试CPU时间 ###
### 20181122 周四 测试CPU时间 ###
HRV-L1.pe.bam 138G数据量 samtools sort ||| 5:29pm开始 6:54pm结束 85分钟10个线程  90G sorted 结果
index ||| 21:07开始  22:16结束  69分钟
#命令如下：
samtools sort -m 8G -o /data2/aoyue/output/bamfile/HRV-L1.pe.sorted.bam -O bam -@ 10 /data2/aoyue/output/bamfile/HRV-L1.pe.bam &
samtools index /data2/aoyue/output/bamfile/HRV-L1.pe.sorted.bam && echo "** index done **" &

[aoyue@lulab1:/data2/aoyue/output/bamfile]time samtools sort -m 8G -o /data2/aoyue/output/bamfile/HRV-L1.pe.sorted.bam -O bam -@ 10 /data2/aoyue/output/bamfile/HRV-L1.pe.bam &
[1] 3612
[aoyue@lulab1:/data2/aoyue/output/bamfile]samtools index /data2/aoyue/output/bamfile/HRV-L1.pe.sorted.bam && echo "** index done **" &
[1] 19245

### 20181127 周二 测试线程数加到50 ###
### 20181127 周二 测试线程数加到50 ###
### 20181127 周二 测试线程数加到50 ###
GRC-L1.pe.bam 133G数据量 samtools sort ||| 3:39pm开始 5:02pm结束 83分钟50个线程  87G sorted 结果 
samtools sort -m 8G -o /data2/aoyue/output/bamfile/GRC-L1.pe.sorted.bam -O bam -@ 50 /data2/aoyue/output/bamfile/GRC-L1.pe.bam && echo '** samtools sort done **' &
再测试一遍
HNSH.pe.bam 132G数据量 samtools sort ||| 5:24pm开始 6:48pm结束 84分钟50个线程  85G sorted 结果 
samtools sort -m 8G -o /data2/aoyue/output/bamfile/HNSH.pe.sorted.bam -O bam -@ 50 /data2/aoyue/output/bamfile/HNSH.pe.bam && echo '** samtools sort done **' &
说明排序的时候，与线程数无关，可能与内存有关系。到底这个排序和哪个有关，待增加内存测试。
再测试一遍
HUN-L1.pe.bam 135G数据量 samtools sort ||| 19:18开始 20:54结束 96分钟10个线程  90G sorted 结果 
samtools sort -m 30G -o /data2/aoyue/output/bamfile/HUN-L1.pe.sorted.bam -O bam -@ 10 /data2/aoyue/output/bamfile/HUN-L1.pe.bam && echo '** samtools sort done **' &

所以，最终决定，用 8G memory, 10 thread.

### 20181127 周二 测试去重复时间   15:39 开始 ###
### 20181127 周二 测试去重复时间   15:39 开始 ###
### 20181127 周二 测试去重复时间   15:39 开始 ###
已失败！samtools markdup /data2/aoyue/output/bamfile/HRV-L1.pe.sorted.bam /data2/aoyue/output/bamfile/HRV-L1.pe.sorted.markdup.bam && echo "** samtools markdup done **"
方案：先抽样，再进行markup命令测试。
### 20181128 周三 抽样测试去重复流程    开始 ####
### 20181128 周三 抽样测试去重复流程    开始 ####
### 20181128 周三 抽样测试去重复流程    开始 ####
#1 抽样 fq.gz  几秒钟完成
[aoyue@lulab1 ~]$ java -jar sample.jar &
[1] 245639 
[aoyue@lulab1 ~]$ 100000 reads are sampled from HRV-L1
这是BWAJiao出口
[1]+  完成                  java -jar sample.jar
#2 比对  3:14pm开始 3:15pm结束 
#3 排序  
1.按名字顺序排序
samtools sort -n -m 1G -o /data2/aoyue/output/testbamfile/HRV-L1.pe.sorted.bam -O bam -@ 10 /data2/aoyue/output/testbamfile/HRV-L1.pe.bam &
2.用fixmate为文件添加ms和MC标签
samtools fixmate -m /data2/aoyue/output/testbamfile/HRV-L1.pe.sorted.bam /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.bam &
3.按照自然顺序进行排序
samtools sort -m 1G -o /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.pos.bam -O bam -@ 10 /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.bam &
4.标记重复
samtools markdup /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.pos.bam /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.pos.mkdup.bam &
5.如果不标记重复，直接去重复，看结果。
samtools markdup -r /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.pos.bam /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.pos.rmdup.bam &
#4 建立索引
samtools index HRV-L1.pe.fixmate.pos.bam
samtools index HRV-L1.pe.fixmate.pos.mkdup.bam
samtools index HRV-L1.pe.fixmate.pos.rmdup.bam

pe.bam
namesort.bam
fixmate.bam
1. rm -f /data2/aoyue/output/bamfile/${a}.pe.bam
fixmate.pos.bam
2. rm -f /data2/aoyue/output/bamfile/${a}.namesort.bam
rmdup.bam
3. rm -f /data2/aoyue/output/bamfile/${a}.fixmate.bam

#小的分析
最终生成6个文件， 名字顺序排序，加标签，自然顺序排序，标记重复。
-rw-rw-r--. 1 aoyue aoyue 27M 11月 28 15:14 HRV-L1.pe.bam
-rw-rw-r--. 1 aoyue aoyue 27M 11月 28 15:22 HRV-L1.pe.namesort.bam
-rw-rw-r--. 1 aoyue aoyue 27M 11月 28 15:27 HRV-L1.pe.fixmate.bam
-rw-rw-r--. 1 aoyue aoyue 27M 11月 28 15:30 HRV-L1.pe.fixmate.pos.bam
-rw-rw-r--. 1 aoyue aoyue 26M 11月 28 15:36 HRV-L1.pe.fixmate.pos.rmdup.bam

#查看各个文件的文件头不同之处：
samtools view -H HRV-L1.pe.bam     
samtools view -H HRV-L1.pe.fixmate.bam  结果是不一样的，多加了一个@HD VN:1.5	SO:queryname信息
samtools view -H HRV-L1.pe.fixmate.pos.bam  结果与上不一直，@HD	VN:1.5	SO:coordinate
#查看大文件的排序顺序 HRV-L1.pe.sorted.bam  @HD	VN:1.5	SO:coordinate 说明是按照坐标排序的
samtools view -h HRV-L1.pe.sorted.bam 
samtools view -h HRV-L1.pe.fixmate.pos.bam | head -n 100
samtools view -h HRV-L1.pe.fixmate.bam | head -n 100
samtools view -h HRV-L1.pe.fixmate.pos.mkdup.bam |  head -n 100
#提取文件比对到1号染色体上的bam文件
samtools view -h HRV-L1.pe.fixmate.pos.rmdup.bam 1  -o HRV-L1.chr1.bam
samtools view -h HRV-L1.chr1.bam > HRV-L1.chr1.sam
samtools view -h HRV-L1.bam 7:335800000-335827000  -o HRV-L1.chr7.bam
#提取43 44号染色体，看看比对到叶绿体和线粒体上的reads有多大 20181203 周一
nohup samtools view -h GRC-L1.rmdup.bam 43 -o GRC-L1.rmdup.chr43.bam &
nohup samtools view -h GRC-L1.rmdup.bam 44 -o GRC-L1.rmdup.chr44.bam &
nohup samtools view -h HM.rmdup.bam 43 -o HM.rmdup.chr43.bam &
nohup samtools view -h HM.rmdup.bam 44 -o HM.rmdup.chr44.bam &
nohup samtools view -h HNSH.rmdup.bam 43 -o HNSH.rmdup.chr43.bam &
nohup samtools view -h HNSH.rmdup.bam 44 -o HNSH.rmdup.chr44.bam &
nohup samtools view -h mergeWheat24SM.bam 43 -o mergeWheat24SM.chr43.bam &
nohup samtools view -h mergeWheat24SM.bam 44 -o mergeWheat24SM.chr44.bam &
#将bam文件转换成 sam文件，比较 HRV-L1.pe.fixmate.pos.bam 和 HRV-L1.pe.fixmate.pos.mkdup.bam 和 HRV-L1.pe.fixmate.pos.rmdup.bam 的差别
samtools view HRV-L1.pe.fixmate.pos.bam > HRV-L1.pe.fixmate.pos.sam
samtools view HRV-L1.pe.fixmate.pos.mkdup.bam > HRV-L1.pe.fixmate.pos.mkdup.sam
samtools view HRV-L1.pe.fixmate.pos.rmdup.bam > HRV-L1.pe.fixmate.pos.rmdup.sam
#将 HRV-L1.pe.fixmate.pos.bam 和 HRV-L1.pe.fixmate.pos.mkdup.bam 和 HRV-L1.pe.fixmate.pos.rmdup.bam 建立index，看看index有什么区别。结果：3个文件的Md5值不相同。
#数据质控(plot-bamstats路径查找； gnuplot下载，yum install gnuplot)
samtools stats -r /data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa HRV-L1.pe.fixmate.pos.mkdup.bam > HRV-L1.bam.bc
perl /data1/programs/samtools-1.8/misc/plot-bamstats -p HRV-L1.bam.bc HRV-L1.bam.bc
samtools stats -r /data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa HRV-L1.pe.fixmate.pos.rmdup.bam > HRV-L1.rmdup.bam.bc
perl /data1/programs/samtools-1.8/misc/plot-bamstats -p HRV-L1.rmdup.bam.bc HRV-L1.rmdup.bam.bc
#打印和检索索引内容
samtools idxstats *.bam 
#将bam文件转化为fastq文件测试
samtools fastq -N -1 GRC-L1.chr001:1-100000_1.fastq.gz GRC-L1.chr001:1-100000.bam
samtools fastq -N -2 GRC-L1.chr001:1-100000_2.fastq.gz GRC-L1.chr001:1-100000.bam
samtools fasta -N -1 GRC-L1.chr001:1-100000_1.fasta.gz GRC-L1.chr001:1-100000.bam
samtools fasta -N -2 GRC-L1.chr001:1-100000_2.fasta.gz GRC-L1.chr001:1-100000.bam



### 正式排序、去重复和建立索引 11.29 周四 脚本书写  ###
### 正式排序、去重复和建立索引 11.29 周四 脚本书写  ###
### 正式排序、去重复和建立索引 11.29 周四 脚本书写  ###

cd /Users/Aoyue/Documents/Jiao/002_script
cat SM.txt | while read a;do echo "samtools sort -n -m 4G -o /data2/aoyue/output/bamfile/${a}.namesort.bam -O bam -@ 1 /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** samtools namesort done **' && samtools fixmate -m /data2/aoyue/output/bamfile/${a}.namesort.bam /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.pe.bam && samtools sort -m 4G -o /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam -O bam -@ 1 /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.namesort.bam && samtools markdup -r /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam /data2/aoyue/output/bamfile/${a}.rmdup.bam && rm -f /data2/aoyue/output/bamfile/${a}.fixmate.bam && samtools index /data2/aoyue/output/bamfile/${a}.rmdup.bam &" ;done > namesort.sh
cat First-1-60SM.txt | while read a;do echo "nohup samtools sort -n -m 20G -o /data2/aoyue/output/bamfile/${a}.namesort.bam -O bam -@ 1 /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** samtools namesort done **' && samtools fixmate -m /data2/aoyue/output/bamfile/${a}.namesort.bam /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.pe.bam && samtools sort -m 20G -o /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam -O bam -@ 1 /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.namesort.bam && samtools markdup -r /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam /data2/aoyue/output/bamfile/${a}.rmdup.bam && rm -f /data2/aoyue/output/bamfile/${a}.fixmate.bam && samtools index /data2/aoyue/output/bamfile/${a}.rmdup.bam" ;done > namesort-onebyone-1thread.sh
先跑前30个样品，一气呵成最终结果。
## 20181130 开始run后30个样品  失败！！！！！！！决定一个一个样开始运行。
nohup sh First31-60rmdup.sh &
[1] 440303
## 20181201下午3点07跑上后30个样本的数据
## 20181204 周二 跑后面几个，在2号集群上
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60
cat First-1-60SM.txt | while read a;do echo "samtools sort -n -m 4G -o /data2/aoyue/bamfile/${a}.namesort.bam -O bam -@ 1 /data2/aoyue/bamfile/${a}.pe.bam && echo '** samtools namesort done **' && samtools fixmate -m /data2/aoyue/bamfile/${a}.namesort.bam /data2/aoyue/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/bamfile/${a}.pe.bam && samtools sort -m 4G -o /data2/aoyue/bamfile/${a}.fixmate.pos.bam -O bam -@ 1 /data2/aoyue/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/bamfile/${a}.namesort.bam && samtools markdup -r /data2/aoyue/bamfile/${a}.fixmate.pos.bam /data2/aoyue/bamfile/${a}.rmdup.bam && rm -f /data2/aoyue/bamfile/${a}.fixmate.bam && samtools index /data2/aoyue/bamfile/${a}.rmdup.bam &" ;done > First-TMM-ZYJWLQ.sh
nohup sh First-TMM-ZYJWLQ.sh &
[1] 191240
[aoyue@lulab2:/data2/aoyue]nohup: 忽略输入并把输出追加到"nohup.out"
[1]+  完成                  nohup sh First-TMM-ZYJWLQ.sh
#时间核查
-rw-rw-r--. 1 aoyue aoyue 133G 11月 29 23:54 GRC-L1.namesort.bam
-rw-rw-r--. 1 aoyue aoyue 138G 11月 30 07:17 GRC-L1.fixmate.bam

## 20181212 周三
#不能同时move很多个文件，必须
[aoyue@lulab1:/data2/aoyue/output/bamfile]mv *fiamate.pos.bam ../bamsorted/
mv: 无法获取"*fiamate.pos.bam" 的文件状态(stat): 没有那个文件或目录
[aoyue@lulab1:/data2/aoyue/output/bamfile]

cd /Users/Aoyue/Documents
cat /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.txt | while read a; do echo "mv /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam /data2/aoyue/output/bamsorted"; done > First-mv-1-60fixmate.pos.bam.sh

#20181215 周六 建立 fixmate.pos.bam的Index文件
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/
cat First-1-60SM.txt | while read a; do echo "samtools index /data2/aoyue/output/bamsorted/${a}.fixmate.pos.bam &";done > First-index-1-60fixmate.pos.bam.sh
nohup sh First-index-1-60fixmate.pos.bam.sh &
[1] 339369

################################ 20181211 周二 拆分bam文件测试 ##############################
################################ 20181211 周二 拆分bam文件测试 ##############################
################################ 20181211 周二 拆分bam文件测试 ##############################
#用java写好脚本，用linux写好全部运行的脚本
samtools view -h mergeWheat24SM.bam 44 -o mergeWheat24SM.chr44.bam &

cd /Users/Aoyue/Documents/
cat /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.txt | while read a; do echo "nohup sh /data2/aoyue/splitbam/${a}_spilt.sh &";done > spilt.sh

sh spilt.sh &
[2] 300644
[2]+  完成                  sh spilt.sh
## 20181211 周二 拆分命令全部运行上

## 20181217 周一 对拆分的bam文件进行index建立
cd /Users/Aoyue/Documents/
cat /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.txt | while read a; do echo "nohup sh /data2/aoyue/splitbam.index/${a}_index.sh &";done > splitbam.index.sh

samtools index /data2/aoyue/output/splitBamfile/040/GRC-L1.chr040.bam
cd /data2/aoyue/output/splitBamfile/040 && samtools index GRC-L1.chr040.bam

### 20181211 周二 查看拆分bam文件是否能进行 mpileup ###
### 20181211 周二 查看拆分bam文件是否能进行 mpileup ###
### 20181211 周二 查看拆分bam文件是否能进行 mpileup ###


#截取1号染色体的一段进行测试
/data2/aoyue/output/spiltBamfile/001/GRC-L1.chr001.bam
samtools view -h GRC-L1.chr001.bam 1:1-100000 -o GRC-L1.chr001:1-100000.bam
#截取失败，因为只适用于有index的bam文件，所以看看是否单条染色体可以建立index----可以
@RG	ID:GRC-L1	PL:illumina	SM:GRC-L1	LB:GRC-L1
[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.
samtools index GRC-L1.chr001.bam
samtools idxstats GRC-L1.chr001.bam 
  结果：0	480980714	0	0
   1	471304005	31998099	33980
samtools view -h GRC-L1.chr001.bam 1:1-100000 -o GRC-L1.chr001:1-100000.bam
#进行mpileup测试
samtools mpileup 
samtools mpileup -A -B -q 30 -Q 10 -f /data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa GRC-L1.chr001.bam -r 1:1-100000 -o GRC-L1.chr001:1-100000.pileup.txt
samtools mpileup -A -B -q 30 -Q 10 -f /data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa HM.chr001.bam -o HM.chr001.pileup.txt

## 20181215 周六
#对000-044文件中的每个文件进行index建立，一共 60*45=2700个文件
写JAVA程序


####################################  第二批数据的测试及运行  ################################
####################################  第二批数据的测试及运行  ################################
####################################  第二批数据的测试及运行  ################################

## 20181127 周二
[root@lulab1 release_1]# nohup cp -Rf ./* /data2/sharedData/Jiao/ABD2/ &
[1] 89180  合计5.8T数据量

## 数据转移 
scp /data2/sharedData/Jiao/ABD2/* aoyue@159.226.116.204:/data2/aoyue/ABD2  

### 20181128 周三 计划：将第二批数据的前30个样品，在集群2上跑，用80个CPU。后30个样品，在集群1上跑。###
####################################  BWA  ##############################################
####################################  BWA  ##############################################
####################################  BWA  ##############################################
#1 建立列表
ls *_1.fq.gz |sed 's/_1.fq.gz//g' > SM.txt 
#2 生成命令在一个脚本里
cat /Users/Aoyue/Documents/Jiao/002_script/Second-1-30list.txt | while read a;do echo "bwa mem -t 80 -R '@RG\tID:${a}\tPL:illumina\tSM:${a}\tLB:${a}' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/aoyue/ABD2/${a}_1.fq.gz /data2/aoyue/ABD2/${a}_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** bwa mapping done **' " ;done > bwa.sh

## 20181209  将第二批的后30个样品继续在2号集群上跑
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60
cat Second-31-60list.txt | while read a;do echo "bwa mem -t 80 -R '@RG\tID:${a}\tPL:illumina\tSM:${a}\tLB:${a}' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/aoyue/ABD2/${a}_1.fq.gz /data2/aoyue/ABD2/${a}_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** bwa mapping done **' " ;done > Second-bwa-31-60.sh
## 20181211 周二 分批跑
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/
cat /Users/Aoyue/Documents/Second-bwa-49-54.txt | while read a;do echo "bwa mem -t 50 -R '@RG\tID:${a}\tPL:illumina\tSM:${a}\tLB:${a}' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/aoyue/ABD2/${a}_1.fq.gz /data2/aoyue/ABD2/${a}_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** bwa mapping done **' " ;done > Second-bwa-49-54.sh
cat /Users/Aoyue/Documents/Second-bwa-55-60.txt | while read a;do echo "bwa mem -t 50 -R '@RG\tID:${a}\tPL:illumina\tSM:${a}\tLB:${a}' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD2/${a}_1.fq.gz /data2/sharedData/Jiao/ABD2/${a}_2.fq.gz | samtools view -S -b - > /data2/aoyue/output2/bamfile/${a}.pe.bam && echo '** bwa mapping done **' " ;done > Second-bwa-55-60.sh


nohup sh Second-bwa-49-54.sh &  在2号集群
[1] 77078
nohup sh Second-bwa-55-60.sh &   在1号集群
[1] 284917

#1号集群上单独运行
bwa mem -t 50 -R '@RG\tID:NXWH\tPL:illumina\tSM:NXWH\tLB:NXWH' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD2/NXWH_1.fq.gz /data2/sharedData/Jiao/ABD2/NXWH_2.fq.gz | samtools view -S -b - > /data2/aoyue/output2/bamfile/NXWH.pe.bam && echo '** bwa mapping done **' &
bwa mem -t 50 -R '@RG\tID:MXM\tPL:illumina\tSM:MXM\tLB:MXM' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD2/MXM_1.fq.gz /data2/sharedData/Jiao/ABD2/MXM_2.fq.gz | samtools view -S -b - > /data2/aoyue/output2/bamfile/MXM.pe.bam && echo '** bwa mapping done **' &
bwa mem -t 50 -R '@RG\tID:MKD-L1\tPL:illumina\tSM:MKD-L1\tLB:MKD-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD2/MKD-L1_1.fq.gz /data2/sharedData/Jiao/ABD2/MKD-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output2/bamfile/MKD-L1.pe.bam && echo '** bwa mapping done **' &
[root@lulab1:/data2/aoyue]nohup sh Secnodbwa-MXM-MKDL1.sh &
[2] 210972

一个样品30被其他任务挤掉，故在203上单独重跑，再转移至204上。失败！
nohup bwa mem -t 30 -R '@RG\tID:GEO-L2\tPL:illumina\tSM:GEO-L2\tLB:GEO-L2' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD2/GEO-L2_1.fq.gz /data2/sharedData/Jiao/ABD2/GEO-L2_2.fq.gz | samtools view -S -b - > /data2/aoyue/GEO-L2.pe.bam && echo "** bwa mapping done **" &
nohup bwa mem -t 30 -R '@RG\tID:GEO-L2\tPL:illumina\tSM:GEO-L2\tLB:GEO-L2' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD2/GEO-L2_1.fq.gz /data2/sharedData/Jiao/ABD2/GEO-L2_2.fq.gz | samtools view -S -b - > /data2/aoyue/GEO-L2.pe.bam && echo "** bwa mapping done **" &
[1] 228063
#将1号集群上的pe.bam文件转移到2号集群上，样品 55-60
scp XYLH.pe.bam YASY.pe.bam YMBF.pe.bam YNSW.pe.bam YZM.pe.bam ZJH.pe.bam aoyue@159.226.116.204:/data2/aoyue/output/bamfile/
scp MXM.pe.bam NXWH.pe.bam aoyue@159.226.116.204:/data2/aoyue/output/bamfile/
#20181217 DKM出问题，重新就行比对
bwa mem -t 100 -R '@RG\tID:DKM\tPL:illumina\tSM:DKM\tLB:DKM' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD2/DKM_1.fq.gz /data2/sharedData/Jiao/ABD2/DKM_2.fq.gz | samtools view -S -b - > /data2/aoyue/output2/bamfile/DKM.pe.bam && echo '** bwa mapping done **' &
[2] 392206
nohup samtools sort -n -m 20G -o /data2/aoyue/output2/bamfile/DKM.namesort.bam -O bam -@ 40 /data2/aoyue/output2/bamfile/DKM.pe.bam &
nohup samtools fixmate -m /data2/aoyue/output2/bamfile/DKM.namesort.bam /data2/aoyue/output2/bamfile/DKM.fixmate.bam &
[1] 438493
nohup samtools sort -m 20G -o /data2/aoyue/output2/bamfile/DKM.fixmate.pos.bam -O bam -@ 40 /data2/aoyue/output2/bamfile/DKM.fixmate.bam &
nohup samtools markdup -r /data2/aoyue/output2/bamfile/DKM.fixmate.pos.bam /data2/aoyue/output2/bamfile/DKM.rmdup.bam && samtools index /data2/aoyue/output2/bamfile/DKM.rmdup.bam &
[1] 457633
nohup samtools index /data2/aoyue/output2/bamfile/DKM.fixmate.pos.bam &
[1] 457463
拆分bam
[1] 105722

########################### 第二批数据·排序、去重复和建立索引 ##########################################
########################### 第二批数据·排序、去重复和建立索引 ##########################################
########################### 第二批数据·排序、去重复和建立索引 ##########################################

## 20181206 前25个可以进行去重复
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60
cat Second-1-30list.txt | while read a;do echo "samtools sort -n -m 4G -o /data2/aoyue/output/bamfile/${a}.namesort.bam -O bam -@ 1 /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** samtools namesort done **' && samtools fixmate -m /data2/aoyue/output/bamfile/${a}.namesort.bam /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.pe.bam && samtools sort -m 4G -o /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam -O bam -@ 1 /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.namesort.bam && samtools markdup -r /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam /data2/aoyue/output/bamfile/${a}.rmdup.bam && rm -f /data2/aoyue/output/bamfile/${a}.fixmate.bam && samtools index /data2/aoyue/output/bamfile/${a}.rmdup.bam &" ;done > Second-namesort-1-25.sh
#26-29可以去重复
cat Second-1-30list.txt | while read a;do echo "samtools sort -n -m 4G -o /data2/aoyue/output/bamfile/${a}.namesort.bam -O bam -@ 1 /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** samtools namesort done **' && samtools fixmate -m /data2/aoyue/output/bamfile/${a}.namesort.bam /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.pe.bam && samtools sort -m 4G -o /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam -O bam -@ 1 /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.namesort.bam && samtools markdup -r /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam /data2/aoyue/output/bamfile/${a}.rmdup.bam && rm -f /data2/aoyue/output/bamfile/${a}.fixmate.bam && samtools index /data2/aoyue/output/bamfile/${a}.rmdup.bam &" ;done > Second-namesort-26-29.sh
#20181215 31-60可以去重复
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/
cat Second-31-60list.txt | while read a;do echo "samtools sort -n -m 4G -o /data2/aoyue/output/bamfile/${a}.namesort.bam -O bam -@ 4 /data2/aoyue/output/bamfile/${a}.pe.bam && echo '** samtools namesort done **' && samtools fixmate -m /data2/aoyue/output/bamfile/${a}.namesort.bam /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.pe.bam && samtools sort -m 4G -o /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam -O bam -@ 4 /data2/aoyue/output/bamfile/${a}.fixmate.bam && rm -f /data2/aoyue/output/bamfile/${a}.namesort.bam && samtools markdup -r /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam /data2/aoyue/output/bamfile/${a}.rmdup.bam && rm -f /data2/aoyue/output/bamfile/${a}.fixmate.bam && samtools index /data2/aoyue/output/bamfile/${a}.rmdup.bam &" ;done > Second-namesort-31-60.sh

nohup sh Second-namesort-31-60.sh &
[1] 408709
#20181217 周一 fixmate.pos.bam 数据转移文件夹
cd /Users/Aoyue/Documents
cat Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.txt | while read a; do echo "mv /data2/aoyue/output/bamfile/${a}.fixmate.pos.bam /data2/aoyue/output/bamsorted/ ";done > Second-mv-fixmate.pos.bam-1-60.sh
nohup sh Second-mv-fixmate.pos.bam-1-60.sh &

#20181217 周一 建立 fixmate.pos.bam的Index文件
cd /Users/Aoyue/Documents
cat Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.txt | while read a; do echo "samtools index /data2/aoyue/output/bamsorted/${a}.fixmate.pos.bam &";done > Second-index-1-60fixmate.pos.bam.sh
nohup sh Second-index-1-60fixmate.pos.bam.sh &
[1] 20046

########################### 第二批数据·拆分bam文件 ##########################################
########################### 第二批数据·拆分bam文件 ##########################################
########################### 第二批数据·拆分bam文件 ##########################################

nohup sh Second-mv-fixmate.pos.bam-1-60.sh &

cd /Users/Aoyue/Documents
cat Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.txt | while read a; do echo "nohup sh /data2/aoyue/splitbam/${a}_split.sh &";done > Second-split.sh

nohup sh Second-split.sh &
[1] 19827


########################### 第二批数据·计算文件大小 ##########################################
########################### 第二批数据·计算文件大小 ##########################################
########################### 第二批数据·计算文件大小 ##########################################
#原始数据Fq 排序的bam和bai 去重复的bam和bai 
java -jar calFileSize.jar /data2/aoyue/ABD2/ /data2/aoyue/calSecondFq.txt .fq.gz GB
java -jar calFileSize.jar /data2/aoyue/output/bamsorted/ /data2/aoyue/calSecondFixmate.pos.bam.txt .fixmate.pos.bam GB
java -jar calFileSize.jar /data2/aoyue/output/bamsorted/ /data2/aoyue/calSecondFixmate.pos.bam.bai.txt .fixmate.pos.bam.bai MB
java -jar calFileSize.jar /data2/aoyue/output/bamfile /data2/aoyue/calSecondRmdup.bam.txt .rmdup.bam GB
java -jar calFileSize.jar /data2/aoyue/output/bamfile /data2/aoyue/calSecondRmdup.bam.bai.txt .rmdup.bam.bai MB
#拆分bam数据
/Users/Aoyue/Documents/待运行的命令/Java-jar-script.txt


##########################################################################################
##################       20181210 周一 计划要做的事情：拷贝数据到移动硬盘      ##################
##########################################################################################
1.将raw数据 bam数据分类整理-----正在进行时ing 
2.没有拷贝到移动硬盘上存储的，进行数据拷贝。

#1 数据合并
[aoyue@lulab2:/data2/aoyue/bamfile]scp ./* aoyue@159.226.116.203:/data2/aoyue/output/bamfile
#2 数据拷贝 -- 将第二批原始数据拷贝到移动硬盘上(备份2次)
cd /Users/Aoyue/Documents
#第一个硬盘 LuLab3T_50
cat /Users/Aoyue/Documents/Second1-50list.txt | while read a; do echo "cp -Rf /data2/sharedData/Jiao/ABD2/${a} /mnt/usb/ABD004_1-25/";done > Second-cp1-25.sh
#第二个硬盘 LuLab3T_51
cat /Users/Aoyue/Documents/Second51-100list.txt | while read a; do echo "cp -Rf /data2/sharedData/Jiao/ABD2/${a} /mnt/usb/ABD004_26-50/";done > Second-cp26-50.sh
#第三个硬盘 LuLab3T_13? ABD003上剩余的2T
cat /Users/Aoyue/Documents/Secnod101-120list.txt | while read a; do echo "cp -Rf /data2/sharedData/Jiao/ABD2/${a} /mnt/usb/ABD004_51-60/";done > Second-cp51-60.sh
#第二批数据文件转移 20181220 周四 4：57pm进行
scp -r /data2/aoyue/output/splitBamfile/* aoyue@159.226.116.203:/data2/aoyue/splitBamfile_Second/

#所有数据的md5生成 
203集群上
nohup md5sum /data2/aoyue/output/bamsorted/* > /data2/aoyue/output/fixmate.pos.bam.md5.txt &
nohup md5sum /data2/aoyue/output/bamfile/* > /data2/aoyue/output/ABD001_rmdup.bam.md5.txt &

[aoyue@lulab1 output]$ nohup md5sum /data2/aoyue/output/bamsorted/* > /data2/aoyue/output/fixmate.pos.bam.md5.txt &
[1] 417082
[aoyue@lulab1 output]$ nohup: 忽略输入重定向错误到标准输出端

[aoyue@lulab1 output]$ nohup md5sum /data2/aoyue/output/bamfile/* > /data2/aoyue/output/ABD001_rmdup.bam.md5.txt &
[2] 418458

204集群上
nohup md5sum /data2/aoyue/output/bamsorted/* > /data2/aoyue/output/ABD002_fixmate.pos.bam.md5.txt &
nohup md5sum /data2/aoyue/output/bamfile/* > /data2/aoyue/output/ABD002_rmdup.bam.md5.txt &


fixmate.pos.bam
rmdup.bam
split.bam

#
/data2/aoyue/splitBamfile_Second/
/data2/sharedData/Jiao/splitBamfile/

cd /Users/Aoyue/Documents
cat /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.txt | while read a; do echo "sh mv_secondSplitbam/${a}-mv_secondSplitbam.sh"; done > mv_secondSplitbam.sh

### 第一批数据的拷贝 ###
### 第一批数据的拷贝 ###
### 第一批数据的拷贝 ###

## 20181215 周六 fixmate.pos.bam数据拷贝
#第一个硬盘 LuLab3T_52   ABD001_fixmate.pos.bam_1-28/  (已完成)
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/cp/
cat First-1-25SM.txt | while read a; do echo "cp -Rf /data2/aoyue/output/bamsorted/${a}.fixmate.pos.bam /mnt/usb/ABD001_fixmate.pos.bam_1-25/";done > First-cp_fimxmate.pos.bam_1-25.sh
[root@lulab1:/data2/aoyue]nohup sh First-cp_fimxmate.pos.bam_1-25.sh &
[1] 335755
由于硬盘有剩余，故再拷贝3个bam文件。
cp -Rf /data2/aoyue/output/bamsorted/NLD-C1.fixmate.pos.bam /mnt/usb/ABD001_fixmate.pos.bam_1-25/ && cp -Rf /data2/aoyue/output/bamsorted/LZT.fixmate.pos.bam /mnt/usb/ABD001_fixmate.pos.bam_1-25/ && cp -Rf /data2/aoyue/output/bamsorted/MEX-L1.fixmate.pos.bam /mnt/usb/ABD001_fixmate.pos.bam_1-25/ &
[2] 377326
#拷贝1-28的fixmate.pos.bam.bai文件
        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
        String outfileS = doc + "First-cp_fixmate.pos.bam_1-28.sh"; //写出的脚本文件路径 (需要改)
        String infileDirS = "/data2/aoyue/output/bamsorted/"; //脚本中的命令，查找数据的地方 (需要改)
        String output = "ABD001_fixmate.pos.bam_1-28/"; //脚本中的命令，目标文件夹 （需要改）
        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
        String suffix = ".fixmate.pos.bam.bai"; //脚本中的命令，list中的后缀名

nohup sh First-cp_fixmate.pos.bam_1-28.sh &
[2] 391039
        
#第二个硬盘 LuLab3T_53   ABD001_fixmate.pos.bam_29-55/   (20181217 周一 9：20pm)
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/
cat First-1-60SM.txt | while read a; do echo "cp -Rf /data2/aoyue/output/bamsorted/${a}.fixmate.pos.bam /mnt/usb/ABD001_fixmate.pos.bam_29-55/ && cp -Rf /data2/aoyue/output/bamsorted/${a}.fixmate.pos.bam.bai /mnt/usb/ABD001_fixmate.pos.bam_29-55/";done > First-cp_fixmate.pos.bam_29-55.sh

        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
        String outfileS = doc + "First-cp_fixmate.pos.bam_29-55.sh"; //写出的脚本文件路径 (需要改)
        String infileDirS = "/data2/aoyue/output/bamsorted/"; //脚本中的命令，查找数据的地方 (需要改)
        String output = "ABD001_fixmate.pos.bam_29-55/"; //脚本中的命令，目标文件夹 （需要改）
        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
        String suffix = ".fixmate.pos.bam"; //脚本中的命令，list中的后缀名
        
mkdir /mnt/usb/ABD001_fixmate.pos.bam_29-55
nohup sh First-cp_fixmate.pos.bam_29-55.sh &
[2] 391941
#第三个硬盘 LuLab3T_54   ABD001_fixmate.pos.bam_56-60/
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/
cat First-1-60SM.txt | while read a; do echo "cp -Rf /data2/aoyue/output/bamsorted/${a}.fixmate.pos.bam /mnt/usb/ABD001_fixmate.pos.bam_56-60/ && cp -Rf /data2/aoyue/output/bamsorted/${a}.fixmate.pos.bam.bai /mnt/usb/ABD001_fixmate.pos.bam_56-60/";done > First-cp_fixmate.pos.bam_56-60.sh
mkdir /mnt/usb/ABD001_fixmate.pos.bam_56-60
nohup sh First-cp_fixmate.pos.bam_56-60.sh &
[1] 457392
## 20181217 周一 rmdup.bam数据拷贝
#第一个硬盘 LuLab3T_54   ABD001_rmdup.bam_1-25/
cd /Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/
cat First-1-60SM.txt | while read a; do echo "cp -Rf /data2/aoyue/output/bamfile/${a}.rmdup.bam /mnt/usb/ABD001_rmdup.bam_1-25/ && cp -Rf /data2/aoyue/output/bamfile/${a}.rmdup.bam.bai /mnt/usb/ABD001_rmdup.bam_1-25/";done > First-cp_rmdup.bam_1-25.sh

mkdir /mnt/usb/ABD001_rmdup.bam_1-25/
nohup sh First-cp_rmdup.bam_1-25.sh &
[1] 1412
#第二个硬盘 LuLab3T_55   ABD001_rmdup.bam_26-55/ （已做完）
 public void cp(){
        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
        String outfileS = doc + "First-cp_rmdup.bam_26-55.sh"; //写出的脚本文件路径 (需要改)
        String infileDirS = "/data2/aoyue/output/bamfile/"; //脚本中的命令，查找数据的地方
        String output = "ABD001_rmdup.bam_26-55/"; //脚本中的命令，目标文件夹 （需要改）
        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
        String suffix = ".rmdup.bam"; //脚本中的命令，list中的后缀名
mkdir /mnt/usb/ABD001_rmdup.bam_26-55/
nohup sh First-cp_rmdup.bam_26-55.sh &

#第三个硬盘 LuLab3T_56   ABD001_rmdup.bam_56-60/ （该拷贝这个了）
 public void cp(){
        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
        String outfileS = doc + "First-cp_rmdup.bam_56-60.sh"; //写出的脚本文件路径 (需要改)
        String infileDirS = "/data2/aoyue/output/bamfile/"; //脚本中的命令，查找数据的地方
        String output = "ABD001_rmdup.bam_56-60/"; //脚本中的命令，目标文件夹 （需要改）
        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
        String suffix = ".rmdup.bam"; //脚本中的命令，list中的后缀名
mkdir /mnt/usb/ABD001_rmdup.bam_56-60/
nohup sh First-cp_rmdup.bam_56-60.sh &

### 第二批数据的拷贝 ###
### 第二批数据的拷贝 ###
### 第二批数据的拷贝 ###

## 20190103 周四 fixmate.pos.bam数据拷贝
#第一个硬盘 LuLab3T_56   ABD002_fixmate.pos.bam_1-24/  
#拷贝1-24的fixmate.pos.bam文件
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "Second-cp_fixmate.pos.bam_1-24.sh"; //写出的脚本文件路径 (需要改)
//        String infileDirS = "/data2/aoyue/output/bamsorted/"; //脚本中的命令，查找数据的地方 (需要改)
//        String output = "ABD002_fixmate.pos.bam_1-24/"; //脚本中的命令，目标文件夹 （需要改）
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".fixmate.pos.bam"; //脚本中的命令，list中的后缀名
mkdir /mnt/usb/ABD002_fixmate.pos.bam_1-24/
nohup sh Second-cp_fixmate.pos.bam_1-24.sh &

#第二个硬盘 LuLab3T_57   ABD002_fixmate.pos.bam_25-52/  
#拷贝25-52的fixmate.pos.bam文件
     
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "Second-cp_fixmate.pos.bam_25-52.sh"; //写出的脚本文件路径 (需要改)
//        String infileDirS = "/data2/aoyue/output/bamsorted/"; //脚本中的命令，查找数据的地方 (需要改)
//        String output = "ABD002_fixmate.pos.bam_25-52/"; //脚本中的命令，目标文件夹 （需要改）
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".fixmate.pos.bam"; //脚本中的命令，list中的后缀名

mkdir /mnt/usb/ABD002_fixmate.pos.bam_25-52/
nohup sh Second-cp_fixmate.pos.bam_25-52.sh &  [1] 353322 2019.1.4

#第三个硬盘 LuLab3T_58   ABD002_fixmate.pos.bam_53-59/  
#拷贝53-59的fixmate.pos.bam文件
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "Second-cp_fixmate.pos.bam_53-59.sh"; //写出的脚本文件路径 (需要改)
//        String infileDirS = "/data2/aoyue/output/bamsorted/"; //脚本中的命令，查找数据的地方 (需要改)
//        String output = "ABD002_fixmate.pos.bam_53-59/"; //脚本中的命令，目标文件夹 （需要改）
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".fixmate.pos.bam"; //脚本中的命令，list中的后缀名

mkdir /mnt/usb/ABD002_fixmate.pos.bam_53-59/
nohup sh Second-cp_fixmate.pos.bam_53-59.sh &   [1] 358293

## 2019010？ 周？ rmdup.bam数据拷贝
#第一个硬盘 LuLab3T_58   ABD002_rmdup.bam_1-24/
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "Second-cp_rmdup.bam_1-24.sh"; //写出的脚本文件路径
//        String infileDirS = "/data2/aoyue/output/bamfile/"; //脚本中的命令，查找数据的地方
//        String output = "ABD002_rmdup.bam_1-24/"; //脚本中的命令，目标文件夹
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".rmdup.bam"; //脚本中的命令，list中的后缀名
mkdir /mnt/usb/ABD002_rmdup.bam_1-24/
nohup sh Second-cp_rmdup.bam_1-24.sh &    [1] 359407

#第二个硬盘 LuLab3T_59   ABD002_rmdup.bam_25-56
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "Second-cp_rmdup.bam_25-56.sh"; //写出的脚本文件路径
//        String infileDirS = "/data2/aoyue/output/bamfile/"; //脚本中的命令，查找数据的地方
//        String output = "ABD002_rmdup.bam_25-56/"; //脚本中的命令，目标文件夹
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".rmdup.bam"; //脚本中的命令，list中的后缀名
mkdir /mnt/usb/ABD002_rmdup.bam_25-56/
nohup sh Second-cp_rmdup.bam_25-56.sh &    [1] 20764

#第三个硬盘 LuLab3T_60   ABD002_rmdup.bam_57-59
        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.t.txt"; //读取的list文件路径
        String outfileS = doc + "Second-cp_rmdup.bam_57-59.sh"; //写出的脚本文件路径
        String infileDirS = "/data2/aoyue/output/bamfile/"; //脚本中的命令，查找数据的地方
        String output = "ABD002_rmdup.bam_57-59/"; //脚本中的命令，目标文件夹
        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
        String suffix = ".rmdup.bam"; //脚本中的命令，list中的后缀名
mkdir /mnt/usb/ABD002_rmdup.bam_57-59/
nohup sh Second-cp_rmdup.bam_57-59.sh &     [1] 82020

rm -f GRC-L1.chr001:1-100000.fastq.gz GRC-L1.chr001:1-100000_2.fasta.gz GRC-L1.chr001:1-100000.bam GRC-L1.chr001:1-100000.pileup.txt
rm -f GRC-L1.chr001:1-100000_1.fastq.gz GRC-L1.chr001:1-100000_2.fastq.gz
ls -l | grep 'bam' | wc -l   238

## 20190103 周四 全部split.bam数据拷贝
#第一个硬盘 LuLab3T_60
mkdir JiaoABD_split.bam_md5/
cp -Rf /data2/sharedData/Jiao/splitBamfile/md5 /mnt/usb/JiaoABD_split.bam_md5/ 
mkdir JiaoABD_split.bam_chr001-008/
nohup sh FirstSecond-cp_split.bam_chr001-008.sh &     160215
#第二个硬盘 LuLab3T_61
mkdir JiaoABD_split.bam_md5/
cp -Rf /data2/sharedData/Jiao/splitBamfile/md5 /mnt/usb/JiaoABD_split.bam_md5/ 
mkdir JiaoABD_split.bam_chr009-018/
nohup sh FirstSecond-cp_split.bam_chr009-018.sh &       [1] 166183
#第三个硬盘 LuLab3T_62
mkdir JiaoABD_split.bam_md5/
cp -Rf /data2/sharedData/Jiao/splitBamfile/md5 /mnt/usb/JiaoABD_split.bam_md5/ 
mkdir JiaoABD_split.bam_chr019-029/
nohup sh FirstSecond-cp_split.bam_chr019-029.sh &    [1] 182370
#第四个硬盘 LuLab3T_63
mkdir JiaoABD_split.bam_md5/
cp -Rf /data2/sharedData/Jiao/splitBamfile/md5 /mnt/usb/JiaoABD_split.bam_md5/ 
mkdir JiaoABD_split.bam_chr030-041/
nohup sh FirstSecond-cp_split.bam_chr030-041.sh &    200854
#第二个硬盘 LuLab3T_61
mkdir JiaoABD_split.bam_md5/
cp -Rf /data2/sharedData/Jiao/splitBamfile/md5 /mnt/usb/JiaoABD_split.bam_md5/ 
mkdir JiaoABD_split.bam_chr042-044/
nohup sh FirstSecond-cp_split.bam_chr042-044.sh &      [2] 166219












########################################弃用##############################################
########################################弃用##############################################
########################################弃用##############################################
########################################弃用##############################################

## 20181217 周一 chr000-044.bam数据拷贝 （弃用）
#第一个硬盘 LuLab3T_56   ABD001_split.bam_chr000-016/
cd /mnt/usb
mkdir ABD001_split.bam_chr000-016/
cp Rf /data2/aoyue/output/splitBamfile/016/ /mnt/usb/ABD001_split.bam_chr000-016/
nohup sh First-cp_split.bam_chr000-016.sh &

#第二个硬盘 LuLab3T_57   ABD001_split.bam_chr017-038/
cd /mnt/usb
mkdir ABD001_split.bam_chr017-038/

nohup sh First-cp_split.bam_chr017-038.sh &

#第三个硬盘 LuLab3T_50   ABD001_split.bam_chr039-040/
cd /mnt/usb
mkdir ABD001_split.bam_chr039-040/

nohup sh First-cp_split.bam_chr039-040.sh &

#第四个硬盘 LuLab3T_51   ABD001_split.bam_chr041-044/
cd /mnt/usb
mkdir ABD001_split.bam_chr041-044/

nohup sh First-cp_split.bam_chr041-044.sh &

scp /data2/aoyue/output/splitBamfile/001/* aoyue@159.226.116.203:/data2/aoyue/splitBamfile_Second/001

cat /Users/Aoyue/Documents/chr.txt | while read a; do echo "scp /data2/aoyue/output/splitBamfile/${a}/* aoyue@159.226.116.203:/data2/aoyue/splitBamfile_Second/${a}"; done > /Users/Aoyue/Documents/scp.sh

