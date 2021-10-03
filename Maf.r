#MAF-frequence
setwd("//Users/guoyafei/Downloads/Results/Impute/")
input <- read.table("freq_stat_input.frq", header=TRUE,stringsAsFactors = F)
output <- read.table("freq_stat.frq",header=TRUE,stringsAsFactors = F)
xlab <- c("0.05-0.1","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50")
a <- c(40081,26663,39421,45852,45121,34408,21774,8301,4469)
b <- a/266090

> length(which(output[,5] >= 0.05 &output[,5] < 0.1))
[1] 40081
> length(which(output[,5] >= 0.1 &output[,5] < 0.15))
[1] 26663
> length(which(output[,5] >= 0.15 &output[,5] < 0.2))
[1] 39421
> length(which(output[,5] >= 0.2 &output[,5] < 0.25))
[1] 45852
> length(which(output[,5] >= 0.25 &output[,5] < 0.3))
[1] 45121
> length(which(output[,5] >= 0.3 &output[,5] < 0.35))
[1] 34408
> length(which(output[,5] >= 0.35 &output[,5] < 0.4))
[1] 21774
> length(which(output[,5] >= 0.4 &output[,5] < 0.45))
8301
> length(which(output[,5] >= 0.45 &output[,5] < 0.5))
4469
plot(b,ylim=c(0,0.5),xlim=c(1,9),cex=3)

xlab2 <- c("0.00-0.01","0.01-0.02","0.02-0.03","0.03-0.04","0.04-0.05")
length(which(output[,5] > 0.0 &output[,5] < 0.01))
1631338
length(which(output[,5] >= 0.01 &output[,5] < 0.02))
292504
length(which(output[,5] >= 0.02 &output[,5] < 0.03))
119416
length(which(output[,5] >= 0.03 &output[,5] < 0.04))
45907
length(which(output[,5] >= 0.04 &output[,5] < 0.05))
21405
a <- c(1631338,292504,119416,45907,21405)
b <- a/2110570
plot(b,ylim=c(0,1),xlim=c(1,5),cex=3)




