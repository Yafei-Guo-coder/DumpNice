#安装mumer，biopython packages
nucdiff Chinese_Spring_chr1A.fa Jagger_chr1A.fa ./Out my_prefix CS_Jagger
#无结果
#切断query序列，逢N就断开，再比对
nohup cat jagger_chr1A.fa | seqkit fx2tab | cut -f 2 | sed -r 's/n+/\n/gi' | cat -n | seqkit tab2fx | seqkit replace -p "(.+)" -r "Contig{nr}" > jagger_chr1A_contig.fa &
nucdiff Chinese_Spring_chr1A.fa jagger_chr1A_contig.fa ./Out_contig my_prefix CS_Jagger_contig
