#01.map.r
library(readxl)
a <- read_excel("/Users/guoyafei/RstudioProjects/GitHub/R_Code/07_VMap3/Lulab_germplasm_Info.xlsx", sheet=1, na='NA')
b <- a[,c(4,7,8,10,14)]
