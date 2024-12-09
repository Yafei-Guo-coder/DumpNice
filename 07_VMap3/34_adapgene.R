#adaptative genes
setwd("~/Desktop/adapgene/")
#画图：UpSetR
library(UpSetR)
Name <- c("elevation.txt","prec1.txt","prec2.txt","prec3.txt","prec4.txt","prec5.txt","prec6.txt","prec7.txt","prec8.txt","soil11.txt","soil12.txt","soil13.txt","soil2.txt","soil4.txt","soil5.txt","soil6.txt","soil7.txt","soil8.txt","soil9.txt","solar1.txt","solar2.txt","solar3.txt","temp10.txt","temp11.txt","temp1.txt","temp2.txt","temp3.txt","temp4.txt","temp5.txt","temp6.txt","temp7.txt","temp8.txt","temp9.txt")
p <- list()
for(i in c(1:length(Name))){
  data <- read.table(Name[i],header=F,stringsAsFactors = F)
  p[[i]] <- data
}
listInput <- list(SF2 = p[[1]][,1], H12 = p[[2]][,1], TajimaD = p[[3]][,1])
upset(fromList(listInput),  order.by = "freq",text.scale = c(2),point.size = 3.5, line.size = 2)

data <- read.table("/Users/guoyafei/Desktop/adapgene/rep.jaccard.pos.crosstab", header=T, row.names = 1,stringsAsFactors = F)
corrgram(data, order = F,lower.panel = panel.shade,upper.panel = panel.pie,text.panel = panel.txt,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")

data <- read.table("/Users/guoyafei/Desktop/adapgene/rep.num.pos.crosstab", header=T, row.names = 1,stringsAsFactors = F)
corrgram(data, order = F,lower.panel = panel.shade,upper.panel = panel.pie,text.panel = panel.txt,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")

data <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)
sub <- data[,c(3,6:40)]
sub2 <- sub[,c(1,4,8,11,13,24,34,36)]

corrgram(sub2, order = F,lower.panel = panel.shade,upper.panel = panel.pie,text.panel = panel.txt,gap = 0.1,
         main="Correlogram of environment variables intercorrelations")

library(usdm)
vif(sub)
v1 <- vifcor(sub2,th=0.9)

v1
















