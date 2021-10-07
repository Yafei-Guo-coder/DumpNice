#working directory
#yafei@204:/data1/home/yafei/009_GWAS/Rht/WEGA_out
#extract VCF for 349_WEGA samples from Vmap3.0 dataset
for chr in {001..009}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/Fastcall2/04_Vmap3.0vcf/V.3/chr${chr}.vcf.gz --keep 349_WEGA.txt --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --stdout | bgzip -c > chr${chr}.WEGA.vcf.gz &
done
sleep 2h;
for chr in {010..042}
do
vcftools --gzvcf /data2/yafei/004_Vmap3/Fastcall2/04_Vmap3.0vcf/V.3/chr${chr}.vcf.gz --keep 349_WEGA.txt --min-alleles 2 --max-alleles 2 --maf 0.05 --recode  --recode-INFO-all --stdout | bgzip -c > chr${chr}.WEGA.vcf.gz &
done

#filter the VCF according LD:10 5 0.3
for chr in {001..009}
do
plink --vcf chr${chr}.WEGA.vcf.gz --indep-pairwise 10 5 0.3 --out chr${chr}.WEGA_LD --allow-extra-chr --double-id --autosome-num 42
plink --vcf chr${chr}.WEGA.vcf.gz --make-bed --extract chr${chr}.WEGA_LD.prune.in --out chr${chr}.WEGA_filter.prune.in --double-id --autosome-num 42
plink --bfile chr${chr}.WEGA_filter.prune.in --recode vcf-iid --out chr${chr}.WEGA_LD --autosome-num 42
rm *prune*
rm *nosex
bgzip -c chr${chr}.WEGA_LD.vcf > chr${chr}.WEGA_LD.vcf.gz
done
for chr in {010..042}
do
plink --vcf chr${chr}.WEGA.vcf.gz --indep-pairwise 10 5 0.3 --out chr${chr}.WEGA_LD --allow-extra-chr --double-id --autosome-num 42
plink --vcf chr${chr}.WEGA.vcf.gz --make-bed --extract chr${chr}.WEGA_LD.prune.in --out chr${chr}.WEGA_filter.prune.in --double-id --autosome-num 42
plink --bfile chr${chr}.WEGA_filter.prune.in --recode vcf-iid --out chr${chr}.WEGA_LD --autosome-num 42
rm *prune*
rm *nosex
bgzip -c chr${chr}.WEGA_LD.vcf > chr${chr}.WEGA_LD.vcf.gz
done

#Working directory
#yafei@66:/data1/home/yafei/009_GWAS/WEGA_out
#concat file from 42 to 21
vcf-concat VCF/chr001.WEGA_LD.vcf.gz VCF/chr002.WEGA_LD.vcf.gz | bgzip -c > 1A.WEGA.vcf.gz
vcf-concat VCF/chr003.WEGA_LD.vcf.gz VCF/chr004.WEGA_LD.vcf.gz | bgzip -c > 1B.WEGA.vcf.gz
vcf-concat VCF/chr005.WEGA_LD.vcf.gz VCF/chr006.WEGA_LD.vcf.gz | bgzip -c > 1D.WEGA.vcf.gz
vcf-concat VCF/chr007.WEGA_LD.vcf.gz VCF/chr008.WEGA_LD.vcf.gz | bgzip -c > 2A.WEGA.vcf.gz
vcf-concat VCF/chr009.WEGA_LD.vcf.gz VCF/chr010.WEGA_LD.vcf.gz | bgzip -c > 2B.WEGA.vcf.gz
vcf-concat VCF/chr011.WEGA_LD.vcf.gz VCF/chr012.WEGA_LD.vcf.gz | bgzip -c > 2D.WEGA.vcf.gz
vcf-concat VCF/chr013.WEGA_LD.vcf.gz VCF/chr014.WEGA_LD.vcf.gz | bgzip -c > 3A.WEGA.vcf.gz
vcf-concat VCF/chr015.WEGA_LD.vcf.gz VCF/chr016.WEGA_LD.vcf.gz | bgzip -c > 3B.WEGA.vcf.gz
vcf-concat VCF/chr017.WEGA_LD.vcf.gz VCF/chr018.WEGA_LD.vcf.gz | bgzip -c > 3D.WEGA.vcf.gz
vcf-concat VCF/chr019.WEGA_LD.vcf.gz VCF/chr020.WEGA_LD.vcf.gz | bgzip -c > 4A.WEGA.vcf.gz
vcf-concat VCF/chr021.WEGA_LD.vcf.gz VCF/chr022.WEGA_LD.vcf.gz | bgzip -c > 4B.WEGA.vcf.gz
vcf-concat VCF/chr023.WEGA_LD.vcf.gz VCF/chr024.WEGA_LD.vcf.gz | bgzip -c > 4D.WEGA.vcf.gz
vcf-concat VCF/chr025.WEGA_LD.vcf.gz VCF/chr026.WEGA_LD.vcf.gz | bgzip -c > 5A.WEGA.vcf.gz
vcf-concat VCF/chr027.WEGA_LD.vcf.gz VCF/chr028.WEGA_LD.vcf.gz | bgzip -c > 5B.WEGA.vcf.gz
vcf-concat VCF/chr029.WEGA_LD.vcf.gz VCF/chr030.WEGA_LD.vcf.gz | bgzip -c > 5D.WEGA.vcf.gz
vcf-concat VCF/chr031.WEGA_LD.vcf.gz VCF/chr032.WEGA_LD.vcf.gz | bgzip -c > 6A.WEGA.vcf.gz
vcf-concat VCF/chr033.WEGA_LD.vcf.gz VCF/chr034.WEGA_LD.vcf.gz | bgzip -c > 6B.WEGA.vcf.gz
vcf-concat VCF/chr035.WEGA_LD.vcf.gz VCF/chr036.WEGA_LD.vcf.gz | bgzip -c > 6D.WEGA.vcf.gz
vcf-concat VCF/chr037.WEGA_LD.vcf.gz VCF/chr038.WEGA_LD.vcf.gz | bgzip -c > 7A.WEGA.vcf.gz
vcf-concat VCF/chr039.WEGA_LD.vcf.gz VCF/chr040.WEGA_LD.vcf.gz | bgzip -c > 7B.WEGA.vcf.gz
vcf-concat VCF/chr041.WEGA_LD.vcf.gz VCF/chr042.WEGA_LD.vcf.gz | bgzip -c > 7D.WEGA.vcf.gz
#change file format
for i in {1..7}
do
for j in {A,B,D}
do
run_pipeline.pl -Xmx200g -fork1 -vcf ${i}${j}.WEGA.vcf.gz -export ${i}${j}.WEGA -exportType Hapmap -runfork1
done
done
#GWAS of plant height and stoma for each chromosome
for i in {1..7}
do
for j in {A,B,D}
do
run_pipeline.pl -Xmx200g -fork1 -h ${i}${j}.WEGA.hmp.txt -fork2 -r height.txt -fork3 -q matrix_10.Q -excludeLastTrait -fork4 -k all_kinship.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export height/${i}${j}_height_mlm -runfork1 -runfork2 -runfork3
run_pipeline.pl -Xmx200g -fork1 -h ${i}${j}.WEGA.hmp.txt -fork2 -r stoma.txt -fork3 -q matrix_10.Q -excludeLastTrait -fork4 -k all_kinship.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export stoma/${i}${j}_stoma_mlm -runfork1 -runfork2 -runfork3
done
done

#change result file format
#QQ plot file prepare
for i in `ls *all.txt`; do sed -i '1d' $i; done 
for i in `ls *all.txt`; do wc -l $i ; done | awk '{print "shuf -n "int($1*0.005),$2" | sort -k2,2n -k3,3n > "$2".sub.txt"}'

#Manhattan plot file prepare
ls *mlm2* | xargs -n1 > txt.names
for i in `cat txt.names`
do
awk '{print $2"\t"$3"\t"$4"\t"$7"\t"(-log($7)/log(10))}' $i | awk '{if($5>2 && $4!="NaN") print $0}' | awk '{print $1"\t"$2"\t"$3"\t"$4}' > ${i::-15}.man.txt
done





