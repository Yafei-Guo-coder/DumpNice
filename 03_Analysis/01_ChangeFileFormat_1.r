#------------------------------------------------------------------转换WW原始文件.R-------------------------------------------------------------------------------------------------------------
#读取文件：203服务器
setwd("/data1/home/yafei/SNP/WW")
all <- read.table("Balfourier_et_al_Wheat_Phylogeography_DataS2.tab",header=T,stringsAsFactors=F)
names <- read.table("741names.txt",header=F,stringsAsFactors=F)
convers <- read.table("tabw280kAlleleConversionTable.tab",header=T,stringsAsFactors=F)
geno <- all[,which(colnames(all) %in% names[,1])]
header <- all[,1:5]
data <- rbind(header,geno)
#write.table(alldata,"741.tab",row.names=F,quote=F,sep="\t")
#data <- read.table("741.tab",header=T,stringsAsFactors=F)
#修改染色体号和snp位点
chr <- paste(rep(paste("chr",1:7,sep=""),each=3),rep(c("A","B","D"),times=7), sep = "")
pos <- c(471304005,438720154,452179604, 462376173, 453218924,462216879,454103970,448155269,476235359,452555092,451014251,451004620,453230519,451372872,451901030,452440856,452077197,450509124,450046986,453822637,453812268)
for(i in 1:21){
  data[which(data$chromosome == chr[i]),3]  <- i*2-1
  data[which(data$chromosome == i*2-1 & data$snpPosition > pos[i]),3] <- i*2
  data[which(data$chromosome == i*2 & data$snpPosition > pos[i]),4] <- data[which(data$chromosome == i*2),4] - pos[i]
  orderd <- data[order(data[,3], data[,4]), ]
}
sub <- merge(data[,c(1:4)], convers,by.x="psId",by.y="probesetId",all.x=TRUE)
sub_orderd <- sub[order(sub[,3], sub[,4]), ]
#写42个(.convers)文件
varName <- paste("chr",1:42,".convers",sep="")
for(i in 1:42){
  write.table(sub_orderd[which(sub_orderd$chromosome == i),], varName[i],row.names=F,quote=F,sep="\t")
}

#写42个(.tab)文件
varName <- paste("chr",1:42,".tab",sep="")
for(i in 1:42){
  write.table(orderd[which(orderd$chromosome == i),], varName[i],row.names=F,quote=F,sep="\t")
}
#写42个(.bed)文件
varName <- paste("chr",1:42,".bed",sep="")
for(i in 1:42){
  var <- cbind(orderd[which(orderd$chromosome == i),3],orderd[which(orderd$chromosome == i),4],orderd[which(orderd$chromosome == i),4]+1)
  write.table(var, varName[i],row.names=F,col.names=F,quote=F,sep="\t")
}
