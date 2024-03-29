#Rscript --no-save --no-restore
library(vegan)
library(psych)
library(data.table)
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}
setwd("/usr_storage/syp/data3/RDA")
gen=fread("input/pkdata3.RDA.geno",header=F)
POP=read.table("input/popID.txt",header=F)
LOCI=read.table("input/ID.txt",header=F)
rownames(gen)=as.character(POP$V1)
colnames(gen)=as.character(LOCI$V1)
env <- read.csv("now.csv",head=T)
env$IND <- as.character(env$IND)
######pred=env[,5:23]  1 3 5 13 15 19
pred=env[,c(5,7,9,17,19,23)]
pk.rda <- rda(gen ~ ., data=pred, scale=T)
load.rda <- scores(pk.rda, choices=c(1:6), display="species") 
write.table(load.rda,file="result/loaduniq_rda.txt", quote=F)
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
# if you needed to be very conservative and only identify those loci under very strong selection (i.e., minimize false positive rates), you could increase the number of standard deviations to 3.5.
cand1 <- outliers(load.rda[,1],3)  
cand2 <- outliers(load.rda[,2],3) 
cand3 <- outliers(load.rda[,3],3) 
cand4 <- outliers(load.rda[,4],3) 
cand5 <- outliers(load.rda[,5],3) 
cand6 <- outliers(load.rda[,6],3) 

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))
cand6 <- cbind.data.frame(rep(6,times=length(cand6)), names(cand6), unname(cand6))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4)<- colnames(cand5)<-colnames(cand6) <- c("axis","snp","loading")
write.table(cand1,file="result/cand1.txt",quote=F)
write.table(cand2,file="result/cand2.txt",quote=F)
write.table(cand3,file="result/cand3.txt",quote=F)
write.table(cand4,file="result/cand4.txt",quote=F)
write.table(cand5,file="result/cand5.txt",quote=F)
write.table(cand6,file="result/cand6.txt",quote=F)
ncand <- length(cand1) + length(cand2) + length(cand3)+length(cand4) + length(cand5) + length(cand6)


cand <- rbind(cand1, cand2, cand3,cand4, cand5, cand6)
cand$snp <- as.character(cand$snp)
foo <- matrix(nrow=(ncand), ncol=6) 

colnames(foo) <- c("BIO1","BIO3","BIO5","BIO13","BIO15","BIO19")
gen <- as.data.frame(gen)

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo) 
cand <- cand[!duplicated(cand$snp),]

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,23] <- names(which.max(abs(bar[4:22]))) # gives the variable
  cand[i,24] <- max(abs(bar[4:22]))              # gives the correlation
}

colnames(cand)[23] <- "predictor"
colnames(cand)[24] <- "correlation"
write.table(cand,file="result/cand.txt",quote=F)

pdf("RDA.pdf")
opar<-par(no.readonly=TRUE)
par(mfrow=c(3,2))
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
hist(load.rda[,4], main="Loadings on RDA4")
hist(load.rda[,5], main="Loadings on RDA5")
hist(load.rda[,6], main="Loadings on RDA6")
par(opar)
dev.off()

#### Next, take the intersection of the results obtained by RDA and lfmm
