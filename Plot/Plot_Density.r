install.packages('RIdeogram')

require(RIdeogram)
setwd("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/06_Test/")
gene_density <- read.table("SNP_density.txt", header=T, stringsAsFactors = F)
wheat_karyotype <- read.table("Centrome.txt", header=T, stringsAsFactors = F)
ideogram(karyotype = wheat_karyotype, overlaid = gene_density)
convertSVG("chromosome.svg", device = "png")
