rs#     alleles chrom   pos     strand  assembly#       center  protLSID        assayLSID       panelLSID       QCcode
data <- read.table("output.txt",header=T,stringsAsFactors = F)
pos <- read.table("WorldFilterd.pos",header=T,stringsAsFactors = F)
colnames(data)[4] <- "rs#"
colnames(data)[5] <- "strand"
data[,5] <- NA
colnames(data)[6] <- "chrom"
colnames(data)[7] <- "pos"
colnames(data)[1] <- "assembly#"
colnames(data)[2] <- "center"
colnames(data)[3] <- "protLSID"
colnames(data)[8] <- "assayLSID"

alleles <- paste(pos[,2],pos[,3],sep="/")
data$alleles <- alleles
data$panelLSID <- NA
data$QCcode <- NA
data[,1] <- NA
data[,2] <- NA
data[,3] <- NA
data[,8] <- NA
all <- data[,c(4,755,6,7,5,1,2,3,8,756,757,9:754)]
all2 <- all[which(all$chrom != "NA"),]
write.table(all2,"World.hmp.txt",quote=F,row.names = F,sep="\t")
