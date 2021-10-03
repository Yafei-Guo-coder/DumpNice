library(ggmap)#Load libraries
library(ggplot2)
hpars<-read.table("https://sites.google.com/
site/arunsethuraman1/teaching/hpars.dat?revision=1")


par(mfrow = c(2, 2))
pie(rep(1,4),labels=c("AF","AM","AS","EU"), col = c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"), main = "Branch-Region")
pie(rep(1,7),labels=c("Cultivar","Domesticated_emmer","Free_threshing_tetraploids","Landrace","OtherHexaploid","OtherTetraploid","Wild_emmer"), col = c("#836FFF","#87CEFF","#8B636C","#8DD3C7","#9B30FF","#B3DE69","#B452CD"), main = "Label-Subspecies")
pie(rep(1,5),labels=c("Kong","LuLab","WEGA","ZKS","ZN"), col = c("#8DD3C7","#FB8072","#B3DE69","#8B636C","#80B1D3"), main = "Branch-Source")
pie(rep(1,4),labels=c("ExclusionHexaploid","ExclusionTetraploid","Hexaploid","Tetraploid"), col = c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"), main = "Label-In-Exclu")

