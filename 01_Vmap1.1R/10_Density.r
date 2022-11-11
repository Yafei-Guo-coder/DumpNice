#install.packages('RIdeogram')
require(RIdeogram)

setwd("/Users/guoyafei/Documents/01_Migration/01_BasicStatistic/13_Plots/01_Density")
#gene_density <- read.table("ArinaLrFor_LTR_1.txt", header=T,stringsAsFactors = F)
gene_density2 <- read.table("gene_density.txt", header=T, stringsAsFactors = F)
gene_density <- read.table("21-1M_VMap3_SnpDensity.txt", header=T, stringsAsFactors = F)
wheat_karyotype <- read.table("wheat_karyotype.txt", header=T, stringsAsFactors = F)
ideogram(karyotype = wheat_karyotype, overlaid = gene_density2)
convertSVG("chromosome.svg", device = "pdf")

library(CMplot)
setwd("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据")
mydata<-read.table("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据/fastcall2_001_pos.txt",header=TRUE,sep="\t")
head(mydata)
# snp chr pos
# snp1_1  1 2041
# snp1_2  1 2062
# snp1_3  1 2190
CMplot(mydata,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300) 
mydata<-read.table("/Users/guoyafei/Documents/01_个人项目/02_VmapIII/03_Fastcall2/测试数据/fastcall_001_pos.txt",header=TRUE,sep="\t")
head(mydata)
# snp         chr       pos
# snp1_1    1        2041
# snp1_2    1        2062
# snp1_3    1        2190
CMplot(mydata,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300) 


#lineage density
A <- read.table("Alineage.txt",header=F,stringsAsFactors = F)
ggplot(A, mapping = aes(x = V4)) +
  geom_density( alpha = 0.5, color = "#999999",fill="#999999") +
  #geom_density(fill = "blue",alpha = 0.3,color="blue")+
  theme_bw()
  

B <- read.table("Blineage.txt",header=F,stringsAsFactors = F)
ggplot(B, mapping = aes(x = V4)) +
  geom_density( alpha = 0.5, color = "#377EB8",fill="#377EB8") +
  #geom_density(fill = "blue",alpha = 0.3,color="blue")+
  theme_bw()

D <- read.table("Dlineage.txt",header=F,stringsAsFactors = F)
ggplot(D, mapping = aes(x = V4)) +
  geom_density( alpha = 0.5, color = "#FF7F00",fill="#FF7F00") +
  #geom_density(fill = "blue",alpha = 0.3,color="blue")+
  theme_bw()

ggplot(D, mapping = aes(x = V4),color = "blue",fill="blue") +
  geom_line(stat = "density")



