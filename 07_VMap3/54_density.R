setwd("/Users/guoyafei/Documents/Vmap3/density")
require(RIdeogram)


gene_density <- read.table("all.count.txt", header=T, stringsAsFactors = F)

wheat_karyotype <- read.table("wheat_karyotype.txt", header=T, stringsAsFactors = F)
ideogram(karyotype = wheat_karyotype, overlaid = gene_density)
convertSVG("chromosome.svg", device = "pdf")

#wheat_karyotype
#Chr Start       End  CE_start    CE_end
#1     0 248956422 122026459 124932724
#2     0 242193529  92188145  94090557

#gene_density
#Chr   Start     End Value
#1       1 1000000    65
#1 1000001 2000000    76
library(ggExtra)
setwd("/Users/guoyafei/Desktop/毕业论文/fig/depth")

data <- read.table("v3.0.depth.DD.txt", header=T, stringsAsFactors = F)
p <- ggplot() +
  geom_point(data = data,aes(x=mean,y=sd),alpha=0.1) +
  xlim(0,15)+
  ylim(0,10)+
  theme_bw()
p1 <- ggMarginal(p, type = "density", margins = "both", size = 5, fill = "lightblue", color = "blue")

data <- read.table("v3.2.depth.DD.txt", header=T, stringsAsFactors = F)
p <- ggplot() +
  geom_point(data = data,aes(x=mean,y=sd),alpha=0.1) +
  xlim(6,14)+
  ylim(3,8)+
  theme_bw()
p2 <- ggMarginal(p, type = "density", margins = "both", size = 5, fill = "lightblue", color = "blue")

data <- read.table("v3.0.depth.AB.txt", header=T, stringsAsFactors = F)
p <- ggplot() +
  geom_point(data = data,aes(x=mean,y=sd),alpha=0.1) +
  xlim(0,15)+
  ylim(0,10)+
  theme_bw()
p3 <- ggMarginal(p, type = "density", margins = "both", size = 5, fill = "lightblue", color = "blue")

data <- read.table("v3.2.depth.AB.txt", header=T, stringsAsFactors = F)
p <- ggplot() +
  geom_point(data = data,aes(x=mean,y=sd),alpha=0.1) +
  xlim(4,10)+
  ylim(3,8)+
  theme_bw()
p4 <- ggMarginal(p, type = "density", margins = "both", size = 5, fill = "lightblue", color = "blue")


data <- read.table("v3.0.depth.ABD.txt", header=T, stringsAsFactors = F)
p <- ggplot() +
  geom_point(data = data,aes(x=mean,y=sd),alpha=0.1) +
  xlim(0,15)+
  ylim(0,10)+
  theme_bw()
p5 <- ggMarginal(p, type = "density", margins = "both", size = 5, fill = "lightblue", color = "blue")

data <- read.table("v3.2.depth.ABD.txt", header=T, stringsAsFactors = F)
p <- ggplot() +
  geom_point(data = data,aes(x=mean,y=sd),alpha=0.1) +
  xlim(5,11)+
  ylim(3,7)+
  theme_bw()
p6 <- ggMarginal(p, type = "density", margins = "both", size = 5, fill = "lightblue", color = "blue")

pdf("ori.D.pdf", width=4,height=4)
print(p1)
dev.off()

pdf("filter.D.pdf", width=4,height=4)
print(p2)
dev.off()

pdf("ori.AB.pdf", width=4,height=4)
print(p3)
dev.off()

pdf("filter.AB.pdf", width=4,height=4)
print(p4)
dev.off()

pdf("ori.ABD.pdf", width=4,height=4)
print(p5)
dev.off()

pdf("filter.ABD.pdf", width=4,height=4)
print(p6)
dev.off()

