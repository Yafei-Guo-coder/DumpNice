library(ggplot2)
library(dplyr)
library(ggridges)
library(RColorBrewer)
display.brewer.all()
col1 = c(brewer.pal(9,'Pastel1'),brewer.pal(10,'Set3'))
col1 <- c(rep("#FB8072", times=7),rep("#80B1D3", times=6),rep("#FDB462", times=6))
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density")
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Gene/V5/Fst_density" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F,sep="\t")})
names <- read.table("../getFst.txt",header=F,stringsAsFactors = F)
all <- data[[1]]
all$Pop <- NA
colnames(all)[1:6] = c("Chr","Pos1","Pos2","Site","WEIGHTED_FST","MEAN_FST")
for(i in c(1:length(data))){
  colnames(data[[i]]) = c("Chr","Pos1","Pos2","Site","WEIGHTED_FST","MEAN_FST")
  data[[i]]$Pop <- names[i,1]
  all <- rbind(all,data[[i]])
}
all = all[which(all$MEAN_FST>= 0),]
all = all[which(all$Pop!="NA"),]
unique(all$Pop)
all$Pop = factor(all$Pop,levels=rev(unique(all$Pop)))
ggplot(all, aes(x=MEAN_FST, y=Pop,fill=Pop)) +
  scale_fill_manual(values = col1) +
  geom_density_ridges(scale=3) +
  theme_ridges(grid = FALSE) +
  theme(axis.title.y = element_blank()) +
  xlab("Fst") + ylab("Population")
values = colorspace::sequential_hcl(5, palette = "Peach")
ggplot(data) +
  geom_point(aes(wt, mpg), color = 'red') +
  geom_text(aes(wt, mpg, label = rownames(data))) +
  theme_classic(base_size = 16)
ggplot(data) +
  geom_point(aes(wt, mpg), color = 'red') +
  geom_text_repel(aes(wt, mpg, label = rownames(data))) +
  theme_classic(base_size = 16)
