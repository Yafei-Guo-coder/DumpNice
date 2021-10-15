library(ggplot2)
library(RColorBrewer)
library(aplot)
display.brewer.all()
col <- brewer.pal(n = 8, name = "Spectral")

setwd("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/16_分群/02_gene_region/")
eigvec <- "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/16_分群/02_gene_region/Alineage.gcta.eigenvector.xls"
eigval <- "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/16_分群/02_gene_region/Alineage.gcta.eigenvalue.xls"
popinfo <- "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/16_分群/popinfo.txt"
key <- "Alineage_PCA"
od <- "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/16_分群/02_gene_region/"

poptable <- read.table(popinfo, header = T, comment.char = "",sep="\t")
pop <- unique(poptable[1:1355,])
colnames(pop) <- c("order","exp_id","vcf_id","group","color","pch","group2")
print(pop)
source("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/16_分群/pca.plot2d.r")
pca_plot(eigenvector = eigvec, eigenvalue = eigval,
         group = popinfo, key = key, outdir = od,
         shape = T, shapes = pop$group2, border = T, border_size = 2.5,
         line0 = T, line0_size = 1)

