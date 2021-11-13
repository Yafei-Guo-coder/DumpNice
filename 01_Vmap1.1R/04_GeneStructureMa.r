#Ppd-1和VRN基因随纬度变化的structure分布
#Working diirectory: 203:yafei:/data1/home/yafei/008_Software/snpEff/Xp-clr_6VIP
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)
require(gridExtra)
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure")
#load Qmatrix files
#Ppd-1
sfiles <- list.files(path="/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/Ppd-1", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))
#VRN2-2
sfiles <- list.files(path="/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Structure/VRN2-2", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))









