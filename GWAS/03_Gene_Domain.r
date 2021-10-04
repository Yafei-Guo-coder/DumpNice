###画这个基因的各个部分的分类图
a = 449254000
b = 449260000
plot(NULL, xlim=c(a, b), ylim=c(0, 5), main="",axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
polygon(c(a, a, b,b), c(2.7,3,3,2.7),col="grey", border="grey")
rect(449254955,2.7,449255062,3,col="#00868B",border = NA)
rect(449255063,2.7,449255198,3,col="#00FFFF",border = NA)
rect(449255307,2.7,449255709,3,col="#00FFFF",border = NA)
rect(449256551,2.7,449257487,3,col="#00FFFF",border = NA)
rect(449259163,2.7,449259241,3,col="#00868B",border = NA)

