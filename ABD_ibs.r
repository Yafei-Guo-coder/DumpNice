#服务器工作目录：yafei@203:/data2/yafei/003_project3/Project3/CS_Vmap
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/01_IBS/02_ABD_IBS")
#读取group名的文件
path <- "/Users/guoyafei/Documents/01_个人项目/02_Migration/03_基本统计/14_Subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 
length(data)

#28个
#nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_Durum","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_Khorasan_wheat","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_Polish_wheat","AABB_Rivet_wheat","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")
#25个:7,11,16,17是freethreshing
#nameABD <- c("B_Speltoides","DD_Anathera","AABBDD_Club_wheat","AABBDD_Cultivar","A_Domesticated_einkorn","AABB_Domesticated_emmer","AABB_freethreshing","AABB_Georgian_wheat","AABBDD_Indian_dwarf_wheat","AABB_Ispahanicum","AABB_freethreshing","AABBDD_Landrace","AABBDD_Macha","DD_Meyeri", "AABB_Persian_wheat","AABB_freethreshing","AABB_freethreshing","AABBDD_Spelt","B_Speltoides","DD_Strangulata", "AABBDD_Synthetic","AABBDD_Tibetan_semi_wild","A_Urartu","AABBDD_Vavilovii","A_Wild_einkorn","AABB_Wild_emmer","AABBDD_Xinjiang_wheat","AABBDD_Yunan_wheat")

A <- read.table("Alineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
rownames(A) <- A$Dxy
Alineage <- A[which(A$Dxy %in% c(data[[7]][,1],data[[11]][,1],data[[16]][,1],data[[17]][,1])),which(colnames(A)%in%data[[12]][,1])]
B <- read.table("Blineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
rownames(B) <- B$Dxy
Blineage <- B[which(B$Dxy %in% c(data[[7]][,1],data[[11]][,1],data[[16]][,1],data[[17]][,1])),which(colnames(B)%in%data[[12]][,1])]
D <- read.table("Dlineage_withCS.all.ibs.txt",header=T,stringsAsFactors=F)
rownames(D) <- D$Dxy
Dlineage <- D[which(D$Dxy %in% c(data[[20]][,1])),which(colnames(D)%in%data[[12]][,1])]
Alineage$mean <- apply(Alineage,1,mean)
Blineage$mean <- apply(Blineage,1,mean)
Dlineage$mean <- apply(Dlineage,1,mean)
newD <- Dlineage[-33,]
