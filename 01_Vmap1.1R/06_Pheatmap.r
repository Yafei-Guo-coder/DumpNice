library(cluster)
library(factoextra)
data <- read.table("/Users/guoyafei/Documents/Lulab/傅老师/latlon.txt",header=F,stringsAsFactors = F, fill=TRUE)
colnames(data) <- c("id","Latitude", "Longitude")
rownames(data) <- data$id
data2 <- data[!is.na(data$Latitude),-1]
data3 <- data[is.na(data$Latitude),-1]
#data2 <- data[,-1]
#desc_stats = data.frame( Min=apply(data2, 2, min),#minimum
#                         Med=apply(data2, 2, median),#median
#                         Mean=apply(data2, 2, mean),#mean
#                         SD=apply(data2, 2, sd),#Standard deviation
#                         Max=apply(data2, 2, max))
#desc_stats = round(desc_stats, 1)

df = scale(data2)

#先求样本之间两两相似性
result <- dist(df, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")

data2$type <- cutree(result_hc, k=50)
lat_mean <- tapply(data2[,1],data2$type,mean,na.rm = TRUE)
lon_mean <- tapply(data2[,2],data2$type,mean,na.rm = TRUE)
data2$cluster1 <- NA
data2$cluster2 <- NA
for(i in 1:305) {
  for(j in 1:50){
    if(data2[i,3] == j ){
      data2[i,4] <- as.numeric(lat_mean[j])
      data2[i,5] <- as.numeric(lon_mean[j])
    } 
  }
}
write.table(data2,"data2.txt",quote=F,row.names = T)


#批量画图
#统一设置：注释内容及颜色
library(pheatmap)
annotation_col <- read.table("/Users/guoyafei/Documents/Lulab/傅老师/Annotation_col.txt",header=T,stringsAsFactors = T,fill = TRUE)
rownames(annotation_col) = c(1:349)
labels_col = rep(c(""), 349)
ann_colors = list(
  Year = c(plus = "#1B9E77",minus = "yellow",zero = "firebrick"),
  Growing_Habit = c(Facultative = "yellow", Spring="firebrick",Winter="blue"),
  Region = c(AM = "#1B9E77",EA = "#D95F02",SCA="#7570B3",WA="#E7298A")
)

#批量读取文件夹下面的文件
path <- "/Users/guoyafei/Documents/Lulab/傅老师/Haplo_TXT"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F,stringsAsFactors=F)}) 

#批量修改列名,并储存列的顺序
all <- c(1:349)
for(i in 1:length(data)){
  colnames(data[[i]]) <- c(1:349)
  out <- pheatmap(data[[i]], show_colnames=FALSE, cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors)
  names <- colnames(data[[i]][,out$tree_col[["order"]]])
  all <- cbind(all, names)
  
}

require(reshape)
require (rworldmap)
require(rworldxtra)

for(i in 1:length(data)){
  colnames(data[[i]]) <- c(1:349)
}
pdf("Rht.pdf")
dev.off()
#30
latlon <- read.table("/Users/guoyafei/Documents/Lulab/傅老师/latlon.txt",header=F,stringsAsFactors = F, fill=TRUE)

data<- read.table("/Users/guoyafei/Documents/Lulab/傅老师/Haplo_TXT/30_haplo.txt",header=F,stringsAsFactors = F)
data<- data[1,]
colnames(data) <- c(1:349)
out <- pheatmap(data, show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="Rht")

names <- as.numeric(colnames(data[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==4,)
which(wheat1$id==2,)
wheat1$sub[1:52] <- 1
wheat1$sub[53:288] <- 2
wheat1$sub[289:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="29")
#31
colnames(data[[31]]) <- c(1:349)
out <- pheatmap(data[[31]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="Rht-A1")
names <- as.numeric(colnames(data[[31]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==274,)
which(wheat1$id==251,)
wheat1$sub[1:9] <- 1
wheat1$sub[10:25] <- 2
wheat1$sub[26:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="30")
#32
colnames(data[[32]]) <- c(1:349)
out <- pheatmap(data[[32]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="Rht-B1")
names <- as.numeric(colnames(data[[32]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==21,)
which(wheat1$id==36,)
wheat1$sub[1:5] <- 1
wheat1$sub[6:22] <- 2
wheat1$sub[23:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")

dev.off()
pdf("plot.pdf")
#1
out <- pheatmap(data[[1]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GID-A1")
names <- as.numeric(colnames(data[[1]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==154,)
which(wheat1$id==328,)
wheat1$sub[1:4] <- 1
wheat1$sub[5:10] <- 2
wheat1$sub[11:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="1")
#2
out <- pheatmap(data[[2]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GID-B1")
names <- as.numeric(colnames(data[[2]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==260,)
which(wheat1$id==100,)
wheat1$sub[1:56] <- 1
wheat1$sub[57:81] <- 2
wheat1$sub[82:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="2")
#3
out <- pheatmap(data[[3]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GID-D1")
names <- as.numeric(colnames(data[[3]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==194,)
which(wheat1$id==279,)
wheat1$sub[1:341] <- 1
wheat1$sub[342] <- 2
wheat1$sub[343:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="3")
#4
out <- pheatmap(data[[4]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GID2-A1")
names <- as.numeric(colnames(data[[4]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==306,)
which(wheat1$id==329,)
wheat1$sub[1:74] <- 1
wheat1$sub[75:83] <- 2
wheat1$sub[84:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="4")
#5
out <- pheatmap(data[[5]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GID2-B1")
names <- as.numeric(colnames(data[[5]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==101,)
which(wheat1$id==286,)
wheat1$sub[1:197] <- 1
wheat1$sub[198:335] <- 2
wheat1$sub[336:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="5")

#6
out <- pheatmap(data[[6]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GID2-D1")
names <- as.numeric(colnames(data[[6]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==339,)
which(wheat1$id==136,)
wheat1$sub[1:344] <- 1
wheat1$sub[345] <- 2
wheat1$sub[346:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="6")
#7
out <- pheatmap(data[[7]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GRF4-A1")
names <- as.numeric(colnames(data[[7]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==95,)
which(wheat1$id==51,)
wheat1$sub[1:342] <- 1
wheat1$sub[343] <- 2
wheat1$sub[344:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="7")
#8
out <- pheatmap(data[[8]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GRF4-B1")
names <- as.numeric(colnames(data[[8]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==235,)
which(wheat1$id==251,)
wheat1$sub[1:336] <- 1
wheat1$sub[337:339] <- 2
wheat1$sub[340:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="8")
#9
out <- pheatmap(data[[9]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GRF4-D1")
names <- as.numeric(colnames(data[[9]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==305,)
which(wheat1$id==239,)
wheat1$sub[1:336] <- 1
wheat1$sub[337:340] <- 2
wheat1$sub[341:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="9")

#10
out <- pheatmap(data[[10]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA20ox-1A")
names <- as.numeric(colnames(data[[10]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==248,)
which(wheat1$id==2,)
wheat1$sub[1:5] <- 1
wheat1$sub[6:344] <- 2
wheat1$sub[345:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="10")

#11
colnames(data[[11]]) <- c(1:349)
out <- pheatmap(data[[11]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA20ox-1B")
names <- as.numeric(colnames(data[[11]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==210,)
which(wheat1$id==321,)
wheat1$sub[1:16] <- 1
wheat1$sub[17:128] <- 2
wheat1$sub[129:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="11")

#12
colnames(data[[12]]) <- c(1:349)
out <- pheatmap(data[[12]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA20ox-1D")
names <- as.numeric(colnames(data[[12]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==276,)
which(wheat1$id==226,)
wheat1$sub[1:339] <- 1
wheat1$sub[340:342] <- 2
wheat1$sub[343:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="11")

#13

data13 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/txt文件/13_haplo.txt",header=F,stringsAsFactors = F)
colnames(data13) <- c(1:349)
out <- pheatmap(data13, show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA20ox-D3_1")

names <- as.numeric(colnames(data13[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==144,)
which(wheat1$id==2,)
wheat1$sub[1] <- 1
wheat1$sub[2:346] <- 2
wheat1$sub[347:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="12")

#14
colnames(data[[14]]) <- c(1:349)
out <- pheatmap(data[[14]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA20ox-D3_2")
names <- as.numeric(colnames(data[[14]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==31,)
which(wheat1$id==340,)
wheat1$sub[1:13] <- 1
wheat1$sub[14:16] <- 2
wheat1$sub[17:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
write.table(wheat1,"test.txt",quote=F,row.names = F)
wheat1 <- read.table("test.txt",header=T,stringsAsFactors = F)
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="13")
#15
colnames(data[[15]]) <- c(1:349)
out <- pheatmap(data[[15]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA20ox-D3_3")
names <- as.numeric(colnames(data[[15]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==231,)
which(wheat1$id==289,)
wheat1$sub[1:342] <- 1
wheat1$sub[343] <- 2
wheat1$sub[344:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="14")
#16
colnames(data[[16]]) <- c(1:349)
out <- pheatmap(data[[16]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA20ox-D3_4")
names <- as.numeric(colnames(data[[16]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==339,)
which(wheat1$id==172,)
wheat1$sub[1:346] <- 1
wheat1$sub[347] <- 2
wheat1$sub[348:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="15")
#17
colnames(data[[17]]) <- c(1:349)
out <- pheatmap(data[[17]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox3-A1")
names <- as.numeric(colnames(data[[17]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==291,)
which(wheat1$id==61,)
wheat1$sub[1:256] <- 1
wheat1$sub[257:270] <- 2
wheat1$sub[271:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="16")

#18
colnames(data[[18]]) <- c(1:349)
out <- pheatmap(data[[18]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox3-B1")
names <- as.numeric(colnames(data[[18]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==240,)
which(wheat1$id==278,)
wheat1$sub[1:2] <- 1
wheat1$sub[3] <- 2
wheat1$sub[4:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="17")
#19
colnames(data[[19]]) <- c(1:349)
out <- pheatmap(data[[19]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox3-D1")
names <- as.numeric(colnames(data[[19]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==324,)
which(wheat1$id==115,)
wheat1$sub[1:2] <- 1
wheat1$sub[3:10] <- 2
wheat1$sub[11:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="18")
#20
#colnames(data[[20]]) <- c(1:349)
#out <- pheatmap(data[[20]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="")
#names <- as.numeric(colnames(data[[20]][,out$tree_col[["order"]]]))
#wheat1 <- latlon[names,]
#wheat1$sub <- NA
#colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
#which(wheat1$id==267,)
#which(wheat1$id==145,)
#wheat1$sub[1:7] <- 1
#wheat1$sub[8:79] <- 2
#wheat1$sub[80:349] <- 3
#colnames(wheat1)[4] <- "Sub.population.3"
#wheat1$value <- 1
#wheat1 <- wheat1[,c(2,3,1,4,5)]
#wheat2 = wheat1[!is.na(wheat1$Latitude),]
#wheat = wheat2[!is.na(wheat2$Sub.population.3),]
#wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
#wheat_reshape2 <- as.data.frame(wheat_reshape)
#mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
#        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="19")

#21
#colnames(data[[21]]) <- c(1:349)
#out <- pheatmap(data[[21]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="")
#names <- as.numeric(colnames(data[[21]][,out$tree_col[["order"]]]))
#wheat1 <- latlon[names,]
#wheat1$sub <- NA
#colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
#which(wheat1$id==202,)
#which(wheat1$id==11,)
#wheat1$sub[1:3] <- 1
#wheat1$sub[4:43] <- 2
#wheat1$sub[44:349] <- 3
#colnames(wheat1)[4] <- "Sub.population.3"
#wheat1$value <- 1
#wheat1 <- wheat1[,c(2,3,1,4,5)]
#wheat2 = wheat1[!is.na(wheat1$Latitude),]
#wheat = wheat2[!is.na(wheat2$Sub.population.3),]
#wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
#wheat_reshape2 <- as.data.frame(wheat_reshape)
#mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
#        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="20")

#22
colnames(data[[22]]) <- c(1:349)
out <- pheatmap(data[[22]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-A10")
names <- as.numeric(colnames(data[[22]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==20,)
which(wheat1$id==102,)
wheat1$sub[1:41] <- 1
wheat1$sub[42:51] <- 2
wheat1$sub[52:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="21")
#23
colnames(data[[23]]) <- c(1:349)
out <- pheatmap(data[[23]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS1B02G145600")
names <- as.numeric(colnames(data[[23]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==4,)
which(wheat1$id==22,)
wheat1$sub[1:57] <- 1
wheat1$sub[58:310] <- 2
wheat1$sub[311:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="22")
#24
colnames(data[[24]]) <- c(1:349)
out <- pheatmap(data[[24]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-D10")
names <- as.numeric(colnames(data[[24]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==258,)
which(wheat1$id==152,)
wheat1$sub[1:321] <- 1
wheat1$sub[322:327] <- 2
wheat1$sub[328:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="23")
#25
colnames(data[[25]]) <- c(1:349)
out <- pheatmap(data[[25]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-D7")
names <- as.numeric(colnames(data[[25]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==202,)
which(wheat1$id==11,)
wheat1$sub[1:3] <- 1
wheat1$sub[4:43] <- 2
wheat1$sub[44:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="24")

#26
colnames(data[[26]]) <- c(1:349)
out <- pheatmap(data[[26]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox5-2")
names <- as.numeric(colnames(data[[26]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==267,)
which(wheat1$id==145,)
wheat1$sub[1:7] <- 1
wheat1$sub[8:79] <- 2
wheat1$sub[80:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="25")
#27
colnames(data[[27]]) <- c(1:349)
out <- pheatmap(data[[27]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS7A02G450000")
names <- as.numeric(colnames(data[[27]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==295,)
which(wheat1$id==301,)
wheat1$sub[1:338] <- 1
wheat1$sub[339:340] <- 2
wheat1$sub[341:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="26")
#28
colnames(data[[28]]) <- c(1:349)
out <- pheatmap(data[[28]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-B2")
names <- as.numeric(colnames(data[[28]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==341,)
which(wheat1$id==94,)
wheat1$sub[1:29] <- 1
wheat1$sub[30:31] <- 2
wheat1$sub[32:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="27")
#29
colnames(data[[29]]) <- c(1:349)
out <- pheatmap(data[[29]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-D2")
names <- as.numeric(colnames(data[[29]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==306,)
which(wheat1$id==225,)
wheat1$sub[1:13] <- 1
wheat1$sub[14:263] <- 2
wheat1$sub[264:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="28")
#33
colnames(data[[33]]) <- c(1:349)
out <- pheatmap(data[[33]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS1A02G334400")
names <- as.numeric(colnames(data[[33]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==111,)
which(wheat1$id==69,)
wheat1$sub[1:182] <- 1
wheat1$sub[183:227] <- 2
wheat1$sub[228:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#34
colnames(data[[34]]) <- c(1:349)
out <- pheatmap(data[[34]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-B4")
names <- as.numeric(colnames(data[[34]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==20,)
which(wheat1$id==244,)
wheat1$sub[1:246] <- 1
wheat1$sub[247:268] <- 2
wheat1$sub[269:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#35
colnames(data[[35]]) <- c(1:349)
out <- pheatmap(data[[35]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-D4_5BL")
names <- as.numeric(colnames(data[[35]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==31,)
which(wheat1$id==23,)
wheat1$sub[1:344] <- 1
wheat1$sub[345:346] <- 2
wheat1$sub[347:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#36
colnames(data[[36]]) <- c(1:349)
out <- pheatmap(data[[36]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-A8")
names <- as.numeric(colnames(data[[36]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==31,)
which(wheat1$id==23,)
wheat1$sub[1:165] <- 1
wheat1$sub[166:287] <- 2
wheat1$sub[288:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#37
colnames(data[[37]]) <- c(1:349)
out <- pheatmap(data[[37]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS1B02G420800")
names <- as.numeric(colnames(data[[37]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==31,)
which(wheat1$id==23,)
wheat1$sub[1:257] <- 1
wheat1$sub[258:261] <- 2
wheat1$sub[262:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#38
data38 <- read.table("/Users/guoyafei/Documents/Lulab/付老师/txt文件/38_haplo.txt",header=F,stringsAsFactors = F)
colnames(data38) <- c(1:349)
out <- pheatmap(data38, show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS1D02G400700")
names <- as.numeric(colnames(data38[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==13,)
which(wheat1$id==190,)
wheat1$sub[1:340] <- 1
wheat1$sub[341] <- 2
wheat1$sub[342:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")

#39
colnames(data[[39]]) <- c(1:349)
out <- pheatmap(data[[39]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-A6")
names <- as.numeric(colnames(data[[39]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==281,)
which(wheat1$id==94,)
wheat1$sub[1:297] <- 1
wheat1$sub[298:344] <- 2
wheat1$sub[345:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")

#40
colnames(data[[40]]) <- c(1:349)
out <- pheatmap(data[[40]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS2B02G396000")
names <- as.numeric(colnames(data[[40]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==97,)
which(wheat1$id==33,)
wheat1$sub[1:64] <- 1
wheat1$sub[65:74] <- 2
wheat1$sub[75:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#41
colnames(data[[41]]) <- c(1:349)
out <- pheatmap(data[[41]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS2D02G375300")
names <- as.numeric(colnames(data[[41]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==307,)
which(wheat1$id==127,)
wheat1$sub[1] <- 1
wheat1$sub[2:345] <- 2
wheat1$sub[346:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#42
colnames(data[[42]]) <- c(1:349)
out <- pheatmap(data[[42]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-A11")
names <- as.numeric(colnames(data[[42]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==4,)
which(wheat1$id==51,)
wheat1$sub[1:196] <- 1
wheat1$sub[197:204] <- 2
wheat1$sub[205:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")

#43
colnames(data[[43]]) <- c(1:349)
out <- pheatmap(data[[43]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-B11")
names <- as.numeric(colnames(data[[43]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==8,)
which(wheat1$id==281,)
wheat1$sub[1:113] <- 1
wheat1$sub[114:346] <- 2
wheat1$sub[347:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#44
colnames(data[[44]]) <- c(1:349)
out <- pheatmap(data[[44]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS4D02G271300")
names <- as.numeric(colnames(data[[44]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==11,)
which(wheat1$id==277,)
wheat1$sub[1] <- 1
wheat1$sub[2:3] <- 2
wheat1$sub[4:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#45
colnames(data[[45]]) <- c(1:349)
out <- pheatmap(data[[45]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="TraesCS5A02G543100")
names <- as.numeric(colnames(data[[45]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==18,)
which(wheat1$id==129,)
wheat1$sub[1:8] <- 1
wheat1$sub[9:16] <- 2
wheat1$sub[17:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#46
colnames(data[[46]]) <- c(1:349)
out <- pheatmap(data[[46]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="GA2ox-B13")
names <- as.numeric(colnames(data[[46]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==20,)
which(wheat1$id==99,)
wheat1$sub[1:86] <- 1
wheat1$sub[87] <- 2
wheat1$sub[88:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#47
colnames(data[[47]]) <- c(1:349)
out <- pheatmap(data[[47]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="NGR5_1")
names <- as.numeric(colnames(data[[47]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==237,)
which(wheat1$id==256,)
wheat1$sub[1:3] <- 1
wheat1$sub[4:5] <- 2
wheat1$sub[6:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#48
colnames(data[[48]]) <- c(1:349)
out <- pheatmap(data[[48]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="NGR5_2")
names <- as.numeric(colnames(data[[48]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==223,)
which(wheat1$id==147,)
wheat1$sub[1:148] <- 1
wheat1$sub[149:166] <- 2
wheat1$sub[167:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")

#49
colnames(data[[49]]) <- c(1:349)
out <- pheatmap(data[[49]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="NGR5_3")
names <- as.numeric(colnames(data[[49]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==309,)
which(wheat1$id==11,)
wheat1$sub[1:15] <- 1
wheat1$sub[16:63] <- 2
wheat1$sub[64:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#50
colnames(data[[50]]) <- c(1:349)
out <- pheatmap(data[[50]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="PIF_1")
names <- as.numeric(colnames(data[[50]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==52,)
which(wheat1$id==155,)
wheat1$sub[1:343] <- 1
wheat1$sub[344] <- 2
wheat1$sub[345:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#51
colnames(data[[51]]) <- c(1:349)
out <- pheatmap(data[[51]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="PIF_2")
names <- as.numeric(colnames(data[[51]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==89,)
which(wheat1$id==46,)
wheat1$sub[1:18] <- 1
wheat1$sub[19:32] <- 2
wheat1$sub[33:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
#52
colnames(data[[52]]) <- c(1:349)
out <- pheatmap(data[[52]], show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="PIF_3")
names <- as.numeric(colnames(data[[52]][,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==295,)
which(wheat1$id==75,)
wheat1$sub[1:287] <- 1
wheat1$sub[288:317] <- 2
wheat1$sub[318:349] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="31")
dev.off()


