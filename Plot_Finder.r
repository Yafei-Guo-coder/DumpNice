#Top5% finder result
#yafei@203:/data2/yafei/Project3/Finder/top5

#AABBDD
chr <- c("1A","1B", "1D", "2A", "2B","2D", "3A", "3B", "3D", "4A","4B", "4D", "5A", "5B", "5D","6A", "6B", "6D", "7A", "7B","7D")
subchr <- c("1","1", "1", "2", "2","2", "3", "3", "3", "4","4", "4", "5", "5", "5","6", "6", "6", "7", "7","7")
genome <- c("A","B", "D", "A", "B","D", "A", "B", "D", "A","B", "D", "A", "B", "D","A", "B", "D", "A", "B","D")

for(i in c("Spelt","Macha","Club_wheat","Indian_dwarf","Xinjiang_wheat","Yunan_wheat","Tibetan_semi_wild","Vavilovii","Landrace","EA","EU","SCA","WA")){
  path <- paste("/data2/yafei/Project3/Finder/",i,"/xbsj/merge",sep="")
  fileNames <- dir(path) 
  filePath <- sapply(fileNames, function(x){
    paste(path,x,sep='/')})   
  data <- lapply(filePath, function(x){
    read.table(x, header=T,stringsAsFactors=F)})
  for(j in 1:length(data)){
    data[[j]]$chr <- chr[j]
    data[[j]]$subchr <- subchr[j]
    data[[j]]$genome <- genome[j]
  }
  all <- data.frame()
  for(m in 1:length(data)){
    all <- rbind(all, data[[m]])
  }
  num <- dim(all)[1]*0.05
  top5 <- all[order(-all$LR),][1:num,]
  topdata <- as.data.frame(top5[,5])
  #topdata$start <- NA
  #for (n in seq(dim(top5)[1])) { 
  #	topdata[n,2] <- ifelse(top5[n,1]>1000000, top5[n,1] - 1000000, top5[n,1])
  #}
  topdata$start <- top5$location
  topdata$stop <- top5$location + 1000000
  outFile=paste("top5/",i,"_top5",sep="")
  write.table(topdata,outFile,row.names=F,col.names=F,quote=F,sep="\t")
}

#AABB
chr <- c("1A","1B","2A", "2B","3A", "3B","4A","4B","5A", "5B","6A", "6B","7A", "7B")
subchr <- c("1","1","2", "2","3", "3","4","4","5", "5","6", "6","7", "7")
genome <- c("A","B","A","B","A","B","A","B","A","B","A", "B","A", "B")
for(i in c("Wild_emmer","Domesticated_emmer","Persian_wheat","Polish_wheat","Rivet_wheat","Khorasan_wheat","Durum")){
  path <- paste("/data2/yafei/Project3/Finder/",i,"/xbsj/merge",sep="")
  fileNames <- dir(path) 
  filePath <- sapply(fileNames, function(x){
    paste(path,x,sep='/')})   
  data <- lapply(filePath, function(x){
    read.table(x, header=T,stringsAsFactors=F)})
  for(j in 1:length(data)){
    data[[j]]$chr <- chr[j]
    data[[j]]$subchr <- subchr[j]
    data[[j]]$genome <- genome[j]
  }
  all <- data.frame()
  for(m in 1:length(data)){
    all <- rbind(all, data[[m]])
  }
  num <- dim(all)[1]*0.05
  top5 <- all[order(-all$LR),][1:num,]
  topdata <- as.data.frame(top5[,5])
  #topdata$start <- NA
  #for (n in seq(dim(top5)[1])) { 
  #	topdata[n,2] <- ifelse(top5[n,1]>1000000, top5[n,1] - 1000000, top5[n,1])
  #}
  topdata$start <- top5$location
  topdata$stop <- top5$location + 1000000
  outFile=paste("top5/",i,"_top5",sep="")
  write.table(topdata,outFile,row.names=F,col.names=F,quote=F,sep="\t")
}

#DD
chr <- c("1D","2D","3D","4D","5D","6D","7D")
subchr <- c("1","2","3","4","5","6","7")
genome <- c("D","D","D","D","D","D","D")
for(i in c("Meyeri","Anathera","Strangulata")){
  path <- paste("/data2/yafei/Project3/Finder/",i,"/xbsj/merge",sep="")
  fileNames <- dir(path) 
  filePath <- sapply(fileNames, function(x){
    paste(path,x,sep='/')})   
  data <- lapply(filePath, function(x){
    read.table(x, header=T,stringsAsFactors=F)})
  for(j in 1:length(data)){
    data[[j]]$chr <- chr[j]
    data[[j]]$subchr <- subchr[j]
    data[[j]]$genome <- genome[j]
  }
  all <- data.frame()
  for(m in 1:length(data)){
    all <- rbind(all, data[[m]])
  }
  num <- dim(all)[1]*0.05
  top5 <- all[order(-all$LR),][1:num,]
  topdata <- as.data.frame(top5[,5])
  #topdata$start <- NA
  #for (n in seq(dim(top5)[1])) { 
  #	topdata[n,2] <- ifelse(top5[n,1]>1000000, top5[n,1] - 1000000, top5[n,1])
  #}
  topdata$start <- top5$location
  topdata$stop <- top5$location + 1000000
  outFile=paste("top5/",i,"_top5",sep="")
  write.table(topdata,outFile,row.names=F,col.names=F,quote=F,sep="\t")
}

#AA
chr <- c("1A","2A","3A","4A","5A","6A","7A")
subchr <- c("1","2","3","4","5","6","7")
genome <- c("A","A","A","A","A","A","A")
for(i in c("Wild_Einkorn","Domesticated_einkorn","Urartu")){
  path <- paste("/data2/yafei/Project3/Finder/",i,"/xbsj/merge",sep="")
  fileNames <- dir(path) 
  filePath <- sapply(fileNames, function(x){
    paste(path,x,sep='/')})   
  data <- lapply(filePath, function(x){
    read.table(x, header=T,stringsAsFactors=F)})
  for(j in 1:length(data)){
    data[[j]]$chr <- chr[j]
    data[[j]]$subchr <- subchr[j]
    data[[j]]$genome <- genome[j]
  }
  all <- data.frame()
  for(m in 1:length(data)){
    all <- rbind(all, data[[m]])
  }
  num <- dim(all)[1]*0.05
  top5 <- all[order(-all$LR),][1:num,]
  topdata <- as.data.frame(top5[,5])
  #topdata$start <- NA
  #for (n in seq(dim(top5)[1])) { 
  #	topdata[n,2] <- ifelse(top5[n,1]>1000000, top5[n,1] - 1000000, top5[n,1])
  #}
  topdata$start <- top5$location
  topdata$stop <- top5$location + 1000000
  outFile=paste("top5/",i,"_top5",sep="")
  write.table(topdata,outFile,row.names=F,col.names=F,quote=F,sep="\t")
}

#SS
chr <- c("1B","2B","3B","4B","5B","6B","7B")
subchr <- c("1","2","3","4","5","6","7")
genome <- c("B","B","B","B","B","B","B")
for(i in c("Speltoides")){
  path <- paste("/data2/yafei/Project3/Finder/",i,"/xbsj/merge",sep="")
  fileNames <- dir(path) 
  filePath <- sapply(fileNames, function(x){
    paste(path,x,sep='/')})   
  data <- lapply(filePath, function(x){
    read.table(x, header=T,stringsAsFactors=F)})
  for(j in 1:length(data)){
    data[[j]]$chr <- chr[j]
    data[[j]]$subchr <- subchr[j]
    data[[j]]$genome <- genome[j]
  }
  all <- data.frame()
  for(m in 1:length(data)){
    all <- rbind(all, data[[m]])
  }
  num <- dim(all)[1]*0.05
  top5 <- all[order(-all$LR),][1:num,]
  topdata <- as.data.frame(top5[,5])
  #topdata$start <- NA
  #for (n in seq(dim(top5)[1])) { 
  #	topdata[n,2] <- ifelse(top5[n,1]>1000000, top5[n,1] - 1000000, top5[n,1])
  #}
  topdata$start <- top5$location
  topdata$stop <- top5$location + 1000000
  outFile=paste("top5/",i,"_top5",sep="")
  write.table(topdata,outFile,row.names=F,col.names=F,quote=F,sep="\t")
}

#make_bed.sh
cat $1 | awk 'BEGIN{FS="\t";OFS="\t"}{
	if($1=="1A" && $2 < 471304005)	{print "1",$2,$3}
		else if($1=="1A" && $2 > 471304005){print "2",$2-471304005,$3-471304005}
        else if($1=="2A" && $2 < 462376173){print "7",$2,$3}
        else if($1=="2A" && $2 > 462376173){print "8",$2-462376173,$3-462376173}
        else if($1=="3A" && $2 < 454103970){print "13",$2,$3}
        else if($1=="3A" && $2 > 454103970){print "14",$2-454103970,$3-454103970}
        else if($1=="4A" && $2 < 452555092){print "19",$2,$3}
        else if($1=="4A" && $2 > 452555092){print "20",$2-452555092,$3-452555092}
        else if($1=="5A" && $2 < 453230519){print "25",$2,$3}
        else if($1=="5A" && $2 > 453230519){print "26",$2-453230519,$3-453230519}
        else if($1=="6A" && $2 < 452440856){print "31",$2,$3}
        else if($1=="6A" && $2 > 452440856){print "32",$2-452440856,$3-452440856}
        else if($1=="7A" && $2 < 450046986){print "37",$2,$3}
        else if($1=="7A" && $2 > 450046986){print "38",$2-450046986,$3-450046986}
        else if($1=="1B" && $2 < 438720154){print "3",$2,$3}
		else if($1=="1B" && $2 > 438720154){print "4",$2-438720154,$3-438720154}
		else if($1=="2B" && $2 < 453218924){print "9",$2,$3}
		else if($1=="2B" && $2 > 453218924){print "10",$2-453218924,$3-453218924}
		else if($1=="3B" && $2 < 448155269){print "15",$2,$3}
		else if($1=="3B" && $2 > 448155269){print "16",$2-448155269,$3-448155269}
		else if($1=="4B" && $2 < 451014251){print "21",$2,$3}
		else if($1=="4B" && $2 > 451014251){print "22",$2-451014251,$3-451014251}
		else if($1=="5B" && $2 < 451372872){print "27",$2,$3}
		else if($1=="5B" && $2 > 451372872){print "28",$2-451372872,$3-451372872}
		else if($1=="6B" && $2 < 452077197){print "33",$2,$3}
		else if($1=="6B" && $2 > 452077197){print "34",$2-452077197,$3-452077197}
		else if($1=="7B" && $2 < 453822637){print "39",$2,$3}
		else if($1=="7B" && $2 > 453822637){print "40",$2-453822637,$3-453822637}
		else if($1=="1D" && $2 < 452179604){print "5",$2,$3}
		else if($1=="1D" && $2 > 452179604){print "6",$2-452179604,$3-452179604}
		else if($1=="2D" && $2 < 462216879){print "11",$2,$3}
		else if($1=="2D" && $2 > 462216879){print "12",$2-462216879,$3-462216879}
		else if($1=="3D" && $2 < 476235359){print "17",$2,$3}
		else if($1=="3D" && $2 > 476235359){print "18",$2-476235359,$3-476235359}
		else if($1=="4D" && $2 < 451004620){print "23",$2,$3}
		else if($1=="4D" && $2 > 451004620){print "24",$2-451004620,$3-451004620}
		else if($1=="5D" && $2 < 451901030){print "29",$2,$3}
		else if($1=="5D" && $2 > 451901030){print "30",$2-451901030,$3-451901030}
		else if($1=="6D" && $2 < 450509124){print "35",$2,$3}
		else if($1=="6D" && $2 > 450509124){print "36",$2-450509124,$3-450509124}
		else if($1=="7D" && $2 < 453812268){print "41",$2,$3}
		else if($1=="7D" && $2 > 453812268){print "42",$2-453812268,$3-453812268}
    else{}
    }'

#for i in `ls *top5`; do bash make_bed.sh $i | sort -k1,1n -k2,2n >$i.bed; done

#提取gff基因feature
#for i in `ls *bed`
#do
#bedtools intersect -a $i -b /data1/publicData/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3 -wb | awk '{if($6 == "gene") print $1";"$12}' \
#|awk -F '=|;' '{print $1"\t"$3}' |awk '$1 == 1||$1== 2 ||$1== 7||$1== 8|| $1 == 13||$1== 14 ||$1== 19||$1== 20||$1 == 25||$1== 26 \
#	||$1== 31||$1== 32 || $1 == 37||$1== 38  {print $0}'|awk '{print $2}' > /data2/yafei/Project3/Finder/top5Urartu/anno_subgenome/A/${i::-3}txt

#bedtools intersect -a $i -b /data1/publicData/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3 -wb | awk '{if($6 == "gene") print $1";"$12}' \
#|awk -F '=|;' '{print $1"\t"$3}' |awk '$1 == 3||$1== 4 ||$1== 9||$1== 10|| $1 == 15||$1== 16 ||$1== 21||$1== 22||$1 == 27||$1== 28 \
#	||$1== 33||$1== 34 || $1 == 39||$1== 40  {print $0}'|awk '{print $2}' > /data2/yafei/Project3/Finder/top5Urartu/anno_subgenome/B/${i::-3}txt

#bedtools intersect -a $i -b /data1/publicData/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3 -wb | awk '{if($6 == "gene") print $1";"$12}' \
#|awk -F '=|;' '{print $1"\t"$3}' |awk '$1 == 5||$1== 6 ||$1== 11||$1== 12|| $1 == 17||$1== 18 ||$1== 23||$1== 24||$1 == 29||$1== 30 \
#	||$1== 35||$1== 36 || $1 == 41||$1== 42  {print $0}'|awk '{print $2}'> /data2/yafei/Project3/Finder/top5Urartu/anno_subgenome/D/${i::-3}txt
#done

#for i in `ls *bed`
#do
#bedtools intersect -a $i -b /data1/publicData/wheat/gene/v1.1/wheat_v1.1_Lulab.gff3 -wb | awk '{if($6 == "gene") print $1";"$12}' \
#|awk -F '=|;' '{print $1"\t"$3}' | awk '{print $2}' > /data2/yafei/Project3/Finder/top5/anno_subspecies/${i::-3}txt
#done
