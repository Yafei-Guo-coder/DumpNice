#Required-----
library(RColorBrewer)
#读取Subspecies文件，共12个。
path <- "/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/05_VcfCheck/04_TaxaGroup/03_Subspecies"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors=F)}) 
length(data)
#A_IBS----
#读取ibs文件,并且计算每个subspecies与CS的IBS平均值。
spName <- c("CS", "Cultivar", "Domesticated_emmer", "Free_threshing_tetraploid", "Landrace", "Other_hexaploid", "Other_tetraploid", "Stranglata", "Tibetan_semi_wild", "Wild_emmer")
A <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/IBS/1A.IBS.txt",header=T,stringsAsFactors=F)
#dim(A)
#1355 1356
A$Lineage <- NA
rownames(A) <- A$Dxy
#给每一行赋subspecies的值

for(i in 1:length(data)){
  A[which(A$Dxy %in% data[[i]][,1]),1357] <- spName[i]
}

#提取CS的两列
subA <- cbind(A[,which(colnames(A) %in% data[[1]][,1])],A[,1357])
subA$lineage <- "A"
colnames(subA) <- c("IBS1","IBS2","subspecies","lineage")
#B_IBS----
#读取ibs文件,并且计算每个subspecies与CS的IBS平均值。
spName <- c("CS", "Cultivar", "Domesticated_emmer", "Free_threshing_tetraploid", "Landrace", "Other_hexaploid", "Other_tetraploid", "Stranglata", "Tibetan_semi_wild", "Wild_emmer")
B <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/IBS/1B.IBS.txt",header=T,stringsAsFactors=F)
#dim(B)
#1355 1356
B$Lineage <- NA
rownames(B) <- B$Dxy
#给每一行赋subspecies的值
for(i in 1:length(data)){
  B[which(B$Dxy %in% data[[i]][,1]),1357] <- spName[i]
}
#提取CS的两列
subB <- cbind(B[,which(colnames(B) %in% data[[1]][,1])],B[,1357])
subB$lineage <- "B"
colnames(subB) <- c("IBS1","IBS2","subspecies","lineage")
#D_IBS----
#读取ibs文件,并且计算每个subspecies与CS的IBS平均值。
spName <- c("CS", "Cultivar", "Domesticated_emmer", "Free_threshing_tetraploid", "Landrace", "Other_hexaploid", "Other_tetraploid", "Stranglata", "Tibetan_semi_wild", "Wild_emmer")
D <- read.table("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/12_Vmap3Test/07_Test0720/depth_statics/IBS/1D.IBS.txt",header=T,stringsAsFactors=F)
#dim(D)
#1334 1336
D$Lineage <- NA
rownames(D) <- D$Dxy
#给每一行赋subspecies的值
for(i in 1:length(data)){
  D[which(D$Dxy %in% data[[i]][,1]),1337] <- spName[i]
}
#提取CS的两列
subD <- cbind(D[,which(colnames(D) %in% data[[1]][,1])],D[,1337])
subD$lineage <- "D"
colnames(subD) <- c("IBS1","IBS2","subspecies","lineage")

#画每组比较的boxplot图
all <- rbind(subA,subB,subD)
p <- ggplot(all, aes(subspecies, IBS1)) + 
  geom_boxplot(outlier.shape = NA, outlier.colour = NA, aes(fill=factor(lineage))) +
  scale_fill_brewer(palette = "Accent") +
  guides(fill=guide_legend(titlec = NULL)) +
  scale_x_discrete(limits = c("Wild_emmer","Domesticated_emmer", "Free_threshing_tetraploid", "Other_tetraploid", "Landrace", "Cultivar",  "Other_hexaploid", "Tibetan_semi_wild",  "Stranglata"))+
  ylab("PI") +
  xlab("Lineages") +
  scale_y_continuous(limits=c(0,0.25)) + 
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle = 30, hjust = 0.5, vjust = 0.5), axis.text.y = element_text(size=15))



