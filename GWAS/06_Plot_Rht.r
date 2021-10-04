library(ggplot2)
library(ggseqlogo)
setwd("/Users/guoyafei/Documents/个人项目/傅老师/Rht_maf/")

substrRight <- function(x){
  num = nchar(x)-2
  gsub('[0-9.]', '', substr(x, 6, nchar(x)))
}

names <- c("GID-A1","GID-B1","GID-D1","GID2-A1","GID2-B1","GID2-D1","GRF4-A1","GRF4-B1","GRF4-D1","GA20ox-1A","GA20ox-1B","GA20ox-1D","GA20ox-D3_1","GA20ox-D3_2","GA20ox-D3_3","GA20ox-D3_4","GA2ox3-A1","GA2ox3-B1","GA2ox3-D1","NA ","NA ","GA2ox-A10","B Homologous","GA2ox-D10","GA2ox-D7","GA2ox5-2","A Homologous","GA2ox-B2","GA2ox-D2","Rht","Rht-A1","Rht-B1","A Homologous","GA2ox-B4","GA2ox-D4_5BL","GA2ox-A8","B Homologous","D Homologous","GA2ox-A6","B Homologous","D Homologous","GA2ox-A11","GA2ox-B11","D Homologous","A Homologous","GA2ox-B13","NGR5_1","NGR5_2","NGR5_3","PIF_1","PIF_2","PIF_3")
names <- c("Rht")
#pdf("plot1.pdf",width = 15, height = 8)
#for (i in c(1,2,3,4,5,6,10,11,12,13,14,15,16,17,18,19,22,23,24,25,26,33,34,35,36,37,38,39,40,41,42,43,44,47,48,49)){
pdf("plot2.pdf",width = 60, height = 8)
for (i in c("B1")){
  #for (i in 1:52){
  file1 <- paste("/Users/guoyafei/Documents/个人项目/傅老师/Rht_maf/",i,".logo.seq",sep="")
  fasta_file = read.table(file1,header=F,stringsAsFactors = F)
  fasta = fasta_file[,2]
  file2 <- paste("/Users/guoyafei/Documents/个人项目/傅老师/Rht_maf/",i,".snpEff",sep="")
  snpEff <- read.table(file2,header=F,stringsAsFactors = F,fill=TRUE,sep="\t")
  p1 <- ggseqlogo(fasta,method="prob")+theme(axis.text.x = element_blank())+labs(title = names[i])+theme(plot.title = element_text(hjust = 0.5,size = 25))
  
  if(i==13|| i==38){
    snpEff$V3 <- "T"
  }
  Ref <- as.data.frame(snpEff$V3)
  colnames(Ref) <- "letter"
  Alt <- as.data.frame(snpEff$V4)
  colnames(Alt) <- "letter"
  
  if(dim(snpEff)[2] >6){
    Old <- as.data.frame(substring(snpEff$V7,3,5)) 
    colnames(Old) <- "letter"
    Pos <- as.data.frame(gsub('[a-zA-Z.*]', '', snpEff$V7))
    colnames(Pos) <- "letter"
    New <- as.data.frame(substrRight(snpEff$V7))
    colnames(New) <- "letter"
    all <- rbind(Ref,Alt,Old,Pos,New)
    num <- dim(Ref)[1]
    aln <- data.frame(
      letter = all,
      Type=rep(c("Ref","Alt","Old","Pos","New"), each=num),
      x=rep(1:num,5)
    )
    p2 <- ggplot(aln, aes(x, Type)) +
      geom_text(aes(label=letter,)) +
      scale_y_discrete(limits=c("New","Pos","Old","Alt","Ref"))+
      scale_x_continuous(breaks=1:10, expand = c(0.105, 0)) + xlab('') +
      theme_logo() +
      theme(legend.position = 'none', axis.text.x = element_blank())
  } else{
    all <- rbind(Ref,Alt)
    num <- dim(Ref)[1]
    aln <- data.frame(
      letter = all,
      species=rep(c("Ref","Alt"), each=num),
      x=rep(1:num,2)
    )
    p2 <- ggplot(aln, aes(x, Type)) +
      geom_text(aes(label=letter,),check_overlap = TRUE) +
      scale_y_discrete(limits=c("Alt","Ref"))+
      scale_x_continuous(breaks=1:10, expand = c(0.105, 0)) + xlab('') +
      theme_logo() +
      theme(legend.position = 'none', axis.text.x = element_blank())
    
  }
  snpEff$impact <- NA
  snpEff[which(snpEff$V6=="HIGH"),8] <- 4
  snpEff[which(snpEff$V6=="MODERATE"),8] <- 3
  snpEff[which(snpEff$V6=="LOW"),8] <- 2
  snpEff[which(snpEff$V6=="MODIFILER"),8] <- 1
  
  bp_data <- data.frame(
    x=snpEff$V2,
    Impact=snpEff$impact
  )
  
  p3 <- ggplot(bp_data, aes(factor(x), Impact))+
    geom_bar(stat = "identity", fill="grey")+
    theme_logo()+
    xlab("")+
    theme(axis.text.x = element_text(angle = 90))
  
  suppressMessages(require(cowplot))
  p <- plot_grid(p1,p2,p3,ncol = 1, align = "v")
  print(p)
}
dev.off()

###################################################################################
#统一设置：注释内容及颜色
library(pheatmap)
require(reshape)
require (rworldmap)
require(rworldxtra)
setwd("/Users/guoyafei/Documents/个人项目/傅老师/Rht_maf/")

#VRN-A1-250k
annotation_col <- read.table("Annotation_col.txt",header=T,stringsAsFactors = T)
rownames(annotation_col) = c(1:349)
labels_col = rep(c(""), 349)
ann_colors = list(
  Growing_Habit = c(Facultative = "#228B22", Spring="#D95F02", Winter="#1E90FF"),
  Region = c(AM = "#228B22",EA = "#D95F02",SCA="#1E90FF",WA="#E7298A", IND="#000000", JAP="#FFD700")
)
VRN_A1_250 <- read.table("B1.txt",header=T,stringsAsFactors = F)
latlon <- read.table("data2.txt",header=F,stringsAsFactors = F, fill=TRUE)
#pdf("Rht.pdf")
colnames(VRN_A1_250) <- c(1:349)
pheatmap(VRN_A1_250, clustering_method = "ward.D", show_colnames=TRUE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="Rht-B1")



rownames(latlon) <- c(1:349)
names <- as.numeric(colnames(VRN_A1_250[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==3,)
which(wheat1$id==43,)
wheat1$sub[1:65] <- 1
wheat1$sub[66:88] <- 2
wheat1$sub[89:145] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','2','3'),symbolSize=1,
        zColours=c('red','blue','yellow'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="29")

#VRN-A1-300k
annotation_col <- read.table("Annotation_col.txt",header=T,stringsAsFactors = T)
rownames(annotation_col) = c(1:145)
labels_col = rep(c(""), 145)
ann_colors = list(
  Growing_Habit = c(Facultative = "yellow", Spring="red",Winter="blue"),
  Region = c(AM = "#228B22",EA = "#D95F02",SCA="#1E90FF",WA="#E7298A", EU="#000000", AF="#FFD700",EAF="#8B008B")
)
VRN_A1_300 <- read.table("2_haplo.txt",header=T,stringsAsFactors = F)
#latlon <- read.table("latlon.txt",header=F,stringsAsFactors = F, fill=TRUE)
#pdf("Rht.pdf")
colnames(VRN_A1_300) <- c(1:145)
rownames(latlon) <- c(1:145)
out <- pheatmap(VRN_A1_300, show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors, main="VRN-A1_300k")


names <- as.numeric(colnames(VRN_A1_300[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==96,)
which(wheat1$id==31,)
wheat1$sub[1:119] <- 1
wheat1$sub[120:141] <- 2
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

#VRN-B1
annotation_col <- read.table( "Annotation_col.txt", header=T, stringsAsFactors = T )
rownames(annotation_col) = c(1:145)
labels_col = rep(c(""), 145)
ann_colors = list(
  Growing_Habit = c(Facultative = "yellow", Spring="red",Winter="blue"),
  Region = c(AM = "#228B22",EA = "#D95F02",SCA="#1E90FF",WA="#E7298A", EU="#000000", AF="#FFD700",EAF="#8B008B")
)
VRN <- read.table("3_haplo.txt", header=T, stringsAsFactors = F)
#latlon <- read.table("latlon.txt", header=F, stringsAsFactors = F, fill=TRUE)
#pdf("Rht.pdf")
colnames(VRN) <- c(1:145)
rownames(latlon) <- c(1:145)
out <- pheatmap(VRN, show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, main="VRN-B1")
names <- as.numeric(colnames(VRN[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==139,)
which(wheat1$id==4,)
wheat1$sub[1:2] <- 1
wheat1$sub[3:137] <- 2
wheat1$sub[138:145] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2, xlim=c(-120,140), ylim=c(35,40), nameX="Longitude", nameY="Latitude", nameZs=c('2','3'), symbolSize=1, zColours=c('red','blue'), barOrient='vert', oceanCol="#D1EEEE", landCol="#FFDAB9", main="29")

#VRN-D1
annotation_col <- read.table("Annotation_col.txt",header=T,stringsAsFactors = T)
rownames(annotation_col) = c(1:145)
labels_col = rep(c(""), 145)
ann_colors = list(
  Growing_Habit = c(Facultative = "yellow", Spring="red",Winter="blue"),
  Region = c(AM = "#228B22",EA = "#D95F02",SCA="#1E90FF",WA="#E7298A", EU="#000000", AF="#FFD700",EAF="#8B008B")
)
VRN <- read.table("4_haplo.txt",header=T,stringsAsFactors = F)
latlon <- read.table("data2.txt",header=F,stringsAsFactors = F, fill=TRUE)
#pdf("Rht.pdf")
colnames(VRN) <- c(1:145)
rownames(latlon) <- c(1:145)
out <- pheatmap(VRN, show_colnames=FALSE,legend_breaks = -1:2, legend_labels = c("./.", "0/0", "0/1", "1/1"), cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, main="VRN-D1")
names <- as.numeric(colnames(VRN[,out$tree_col[["order"]]]))
wheat1 <- latlon[names,]
wheat1$sub <- NA
colnames(wheat1)[1:3] <- c("id","Latitude", "Longitude")
which(wheat1$id==7,)
which(wheat1$id==16,)
wheat1$sub[1:136] <- 1
wheat1$sub[137] <- 2
wheat1$sub[138:145] <- 3
colnames(wheat1)[4] <- "Sub.population.3"
wheat1$value <- 1
wheat1 <- wheat1[,c(2,3,1,4,5)]
wheat2 = wheat1[!is.na(wheat1$Latitude),]
wheat = wheat2[!is.na(wheat2$Sub.population.3),]
wheat_reshape <- wheat_reshape <- cast(wheat,Latitude+Longitude~Sub.population.3) 
wheat_reshape2 <- as.data.frame(wheat_reshape)
mapPies(wheat_reshape2,xlim=c(-120,140),ylim=c(35,40),nameX="Longitude",nameY="Latitude",nameZs=c('1','3'),symbolSize=1,
        zColours=c('red','blue'),barOrient='vert',oceanCol="#D1EEEE",landCol="#FFDAB9",main="29")