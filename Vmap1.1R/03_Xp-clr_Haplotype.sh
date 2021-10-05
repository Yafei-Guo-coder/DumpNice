#working directory
#203:yafei:/data2/yafei/003_Project3/Vmap1.1/E6/Xp-CLR_V2

#extract gene_region VCF from 225 all vcf files
grep -w -f 87Gene_id.txt gene_v1.1_Lulab.gff3 | awk -F";" '{print $1}' |awk '{print $1"\t"$4-5000"\t"$5+5000"\t"$9}' | awk -F"\tID=" '{print $2"\t"$1}' > 87gene_5k.txt
#bash getVcf.sh 87gene_5k.txt
cat $1 |while read gene chr from to
do
    #echo $chr $from $to
    #if echo $2 |grep -q '.*.vcf.gz$';then
    #    vcftools --gzvcf $2 --chr $chr --from-bp $from --to-bp $to  --recode --recode-INFO-all --out $gene.$chr.$from-$to 
    #elif echo $2 |grep -q '.*.vcf$';then
        vcftools --gzvcf /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225/chr$chr.E6_Landrace_locate.vcf.gz --chr $chr --from-bp $from --to-bp $to  --recode --recode-INFO-all --out 5k.$gene.$chr.$from-$to
    #fi
done

#传到本地用Java 将vcf转换成单倍型txt格式，使用07_VCF_Haplotype_Visual.r进行可视化。

