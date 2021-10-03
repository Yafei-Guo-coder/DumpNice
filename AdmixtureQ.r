setwd("/Users/guoyafei/Desktop/Qplot")
library(pophelper)
library(ggplot2)
require(gridExtra)
sfiles <- list.files(path="/Users/guoyafei/Desktop/Qplot/", full.names=T)
sfile <- "/Users/guoyafei/Desktop/Qplot/chrA_LD.3.Q"
slist <- readQ(files=sfiles)
twolabset <- read.delim("/Users/guoyafei/Desktop/labels", header=T,stringsAsFactors=F)
p1 <- plotQ(slist[1],returnplot=T,exportplot=F,quiet=T,basesize=11,
            grplab=onelabset,grplabsize=4,linesize=0.8,pointsize=4)
#slist <- readQ(files=sfiles)
#tabulateQ(qlist=readQ(sfiles))
#summariseQ(tabulateQ(qlist=readQ(sfiles)))

p1 <- plotQ(alignK(sortQ(slist)[1:11]),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,sortind="all",sharedindlab=F,showindlab=T,showyaxis=T,showticks=T)
p2 <- plotQ(alignK(sortQ(slist)[1:11]),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,grplab=twolabset,grplabsize=1.8,linesize=0.8,pointsize=3,sharedindlab=F,sortind="all",ordergrp=TRUE)
pdf("structure2.pdf")
grid.arrange(p2$plot[[1]])
dev.off()
onelabset <- twolabset[,1,drop=FALSE]
p1 <- plotQ(slist[1],returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="all",showyaxis=T,showticks=T)

p2 <- plotQ(sortQ(slist[c(1,2)]),imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="all",sharedindlab=F,showindlab=T,showyaxis=T,showticks=T)

grid.arrange(p1$plot[[1]],p2$plot[[1]],ncol=2)
