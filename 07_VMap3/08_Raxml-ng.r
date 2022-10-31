#Raxml-ng
#yafei@204:/data2/yafei/004_Vmap3/Tree
#run_pipeline.pl -Xms512m -Xmx5g -vcf lineageD.vcf.gz -export lineageD -exportType Phylip_Inter
#raxml-ng --all --msa lineageAB.phy --seed 12346 --model GTR+G --bs-trees 100 --threads 80 

setwd("/Users/guoyafei/Desktop/")
data <- read.table("/Users/guoyafei/Desktop/207_data.txt",header=F,stringsAsFactors = F)
data <- read.table("/Users/guoyafei/Desktop/source_data_107.txt",header=F,stringsAsFactors = F)
library(ggplot2)

ggplot(data, aes(x = V1)) +
  geom_histogram(binwidth = 5, fill = "lightblue", colour = "black")+
  xlim(1899,2010)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 15))+
  ylab("Count")+
  xlab("Source Date (USDA)")


