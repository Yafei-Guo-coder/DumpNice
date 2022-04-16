#selscan
#参考资料
#hapFLK对具有对瓶颈效应或迁移的群体有着更强大的检测选择信号的能力。
#H12对于部分扫荡选择性清除(soft selective sweep)和全扫荡选择性清除(hard selective sweep)都有比较好的检测效率。
#LASSI-Plus: A program written in C++ that implements the likelihood ratio Λ statistic of DeGiorgio and Szpiech (2021) for detecting selective sweeps and inferring their softness using the spatial distribution of distortions in the haplotype frequency spectrum.
#SweepFinder2: A program written in C that can perform genomic scans for recent selective sweeps selection while controlling for background selection and mutation rate variation (DeGiorgio et al. 2016).
#ihs
selscan --ihs --threads 30  --gap-scale 200000 --max-gap 2000000 --maf 0.01 --vcf chr001.phased.imputed.pos.vcf.gz --map chr001.map4 --out chr001 --keep-low-freq
#selscan --ihs --threads 30 --gap-scale 200000 --max-gap 2000000 --maf 0.01  --vcf chr001.phased.imputed.pos.maf005.vcf.gz --map chr001.maf005.map4 --out chr001.maf005 --keep-low-freq

data <- read.table("/Users/guoyafei/公共资源/01_数据资源/iwgsc_refseqv1.0_mapping_data_42chr.txt",header=T,stringsAsFactors = F)
don <- data %>% 
  # Compute chromosome size
  group_by(chromosome) %>% 
  lines(lowess(data$physicalPosition,data$geneticPosition)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

