#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
library(pophelper)
library(ggplot2)
require(gridExtra)
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/02_Structure_map/Select/")
#load Qmatrix files
sfiles <- list.files(path="/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/02_Structure_map/Select/QMatrix", full.names=T)
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))


pdf("Q.pdf",width = 24,height = 8)
p1 <- plotQ(slist[1:7],returnplot=T,exportplot=F,quiet=T,basesize=11,
            sortind="all",showindlab=F,showyaxis=T,showticks=T,sharedindlab=T,clustercol=brewer.pal(7, "Set2"))
grid.arrange(p1$plot[[1]],p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],nrow=5)
dev.off()
