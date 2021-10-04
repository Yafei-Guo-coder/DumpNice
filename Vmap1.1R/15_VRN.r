#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
setwd("/Users/guoyafei/Documents/Lulab/Project-2-Migration/VRN2/")
annotation_col <- read.table("Annotation_col_vmap.txt",header=T,stringsAsFactors = T)
rownames(annotation_col) = c(1:145)
labels_col = rep(c(""), 145)
ann_colors = list(
  Growing_Habit = c(Facultative = "yellow", Spring="red",Winter="blue"),
  Region = c(AM = "#228B22",EA = "#D95F02",SCA="#1E90FF",WA="#E7298A", EU="#000000", AF="#FFD700",EAF="#8B008B")
)
latlon <- read.table("latlon.txt",header=F,stringsAsFactors = F, fill=TRUE)
#VRN-A2-250k
VRN_A1_250 <- read.table("1_haplo.txt",header=T,stringsAsFactors = F)
colnames(VRN_A1_250) <- c(1:145)
rownames(latlon) <- c(1:145)
out <- pheatmap(VRN_A1_250, show_colnames=TRUE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=5, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="VRN-A2")
names <- as.numeric(colnames(VRN_A1_250[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==78,)
which(wheat1$id==79,)
which(wheat1$id==72,)
which(wheat1$id==52,)
which(wheat1$id==141,)
which(wheat1$id==81,)
wheat1$sub[1] <- 1
wheat1$sub[2:4] <- 2
wheat1$sub[5:6] <- 3
wheat1$sub[7:35] <- 4
wheat1$sub[36:145] <- 5
wheat1$sub[48:145] <- 6
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4','5'),symbolSize=1,
        zColours=c('red','blue','yellow','orange','green'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="29")

#VRN-D2-250k
VRN_A1_250 <- read.table("2_haplo.txt",header=T,stringsAsFactors = F)
colnames(VRN_A1_250) <- c(1:145)
rownames(latlon) <- c(1:145)
out <- pheatmap(VRN_A1_250, show_colnames=TRUE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="VRN-D2")
names <- as.numeric(colnames(VRN_A1_250[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==124,)
which(wheat1$id==72,)
wheat1$sub[1:16] <- 1
wheat1$sub[17:21] <- 2
wheat1$sub[22:145] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="29")


#VRN-A3-250k
VRN_A1_250 <- read.table("3_haplo.txt",header=T,stringsAsFactors = F)
colnames(VRN_A1_250) <- c(1:145)
rownames(latlon) <- c(1:145)
out <- pheatmap(VRN_A1_250, show_colnames=TRUE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=5, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="VRN-A1")
names <- as.numeric(colnames(VRN_A1_250[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==41,)
which(wheat1$id==63,)
which(wheat1$id==123,)
which(wheat1$id==24,)

wheat1$sub[1] <- 1
wheat1$sub[2:11] <- 2
wheat1$sub[12:140] <- 3
wheat1$sub[141:145] <- 4
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3','4'),symbolSize=1,
        zColours=c('red','blue','yellow','orange'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="29")

#VRN-B3-250k
VRN_A1_250 <- read.table("4_haplo.txt",header=T,stringsAsFactors = F)
colnames(VRN_A1_250) <- c(1:145)
rownames(latlon) <- c(1:145)
out <- pheatmap(VRN_A1_250, show_colnames=TRUE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="VRN-B3")
names <- as.numeric(colnames(VRN_A1_250[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==61,)
which(wheat1$id==103,)
wheat1$sub[1:45] <- 1
wheat1$sub[46:141] <- 2
wheat1$sub[142:145] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="29")

#VRN-D3-250k
VRN_A1_250 <- read.table("5_haplo.txt",header=T,stringsAsFactors = F)
colnames(VRN_A1_250) <- c(1:145)
rownames(latlon) <- c(1:145)
out <- pheatmap(VRN_A1_250, show_colnames=TRUE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="VRN-D3")
names <- as.numeric(colnames(VRN_A1_250[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==45,)
which(wheat1$id==5,)
wheat1$sub[1:8] <- 1
wheat1$sub[9] <- 2
wheat1$sub[10:145] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="29")

