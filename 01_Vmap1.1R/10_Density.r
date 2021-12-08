#install.packages('RIdeogram')
require(RIdeogram)

setwd("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/13_Plots/01_Density")
#gene_density <- read.table("ArinaLrFor_LTR_1.txt", header=T,stringsAsFactors = F)
gene_density2 <- read.table("gene_density.txt", header=T, stringsAsFactors = F)
gene_density <- read.table("VMap3_SnpDensity.txt", header=T, stringsAsFactors = F)
wheat_karyotype <- read.table("wheat_karyotype.txt", header=T, stringsAsFactors = F)
ideogram(karyotype = wheat_karyotype, overlaid = gene_density)
convertSVG("chromosome.svg", device = "png")

library(CMplot)
setwd("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据")
mydata<-read.table("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据/fastcall2_001_pos.txt",header=TRUE,sep="\t")
head(mydata)
# snp         chr       pos
# snp1_1    1        2041
# snp1_2    1        2062
# snp1_3    1        2190
CMplot(mydata,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300) 
mydata<-read.table("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据/fastcall_001_pos.txt",header=TRUE,sep="\t")
head(mydata)
# snp         chr       pos
# snp1_1    1        2041
# snp1_2    1        2062
# snp1_3    1        2190
CMplot(mydata,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300) 
