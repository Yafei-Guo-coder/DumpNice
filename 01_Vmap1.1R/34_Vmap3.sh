#Vmap3
#yafei@204:/data2/yafei/Vmap3/sort
#204: /data3/wgs/fastq/wildEmmer_S35
zcat /data3/wgs/fastq/wildEmmer_S35/Rawdata/SFS1_1.fq.gz | paste - - - - | sort -k1,1 -t " " | tr '\t' '\n' | gzip -c > /data2/yafei/Vmap3/sort/SFS1_1.sorted.fq.gz
zcat /data3/wgs/fastq/wildEmmer_S35/Rawdata/SFS1_2.fq.gz | paste - - - - | sort -k1,1 -t " " | tr '\t' '\n' | gzip -c > /data2/yafei/Vmap3/sort/SFS1_2.sorted.fq.gz

zcat /data1/home/liping/gasdata/cleandata/CRR061687_f1.filtered.fq.gz |grep "@" |sort |uniq -c > /data2/liping/Vmap3/CRR061687_f1.id
zcat /data1/home/liping/gasdata/cleandata/CRR061687_r2.filtered.fq.gz |grep "@" |sort |uniq -c > /data2/liping/Vmap3/CRR061687_r2.id
zcat /data1/home/liping/gasdata/cleandata/CRR061717_f1.filtered.fq.gz |grep "@" |sort |uniq -c > /data2/liping/Vmap3/CRR061717_f1.id
zcat /data1/home/liping/gasdata/cleandata/CRR061717_r2.filtered.fq.gz |grep "@" |sort |uniq -c > /data2/liping/Vmap3/CRR061717_r2.id






