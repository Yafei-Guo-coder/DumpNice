setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/ABlineage")
a = read.table("AAexontreeIN.eigenvec",head =T)
head(a[,1:4])
av = read.table("AAexontreeIN.eigenval",head =F)
explained_ratio = av[,1]/sum(av[,1])
pca_color = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/treeRangeColor_E2.txt",head =F)
dim(a)
dim(pca_color)
index = match(pca_color$V1,a$IID)
head(explained_ratio)
#par(bg="black")
#par(mfrow=c(2,1), oma=c(4.5,4, 4, 2.5), mar=c(0.5,1,0.5,1), cex=1)
#layout(matrix(c(1:2), 1, 2, byrow = F), widths=rep(1,2), heights=c(rep(1,2),0.05))
#bty="l" 只保留左和下两条边框
plot(a[index,3:4],xlab="PC1 (42.01%)",ylab="PC2 (13.70%)",cex.lab=1.5,mgp=c(2.5,1,0),cex.axis=1,col= as.character(pca_color[,3]),pch=19,bg="transparent",cex=1.5)

#####
rm(list=ls())
library(shape)
library(cluster)
library(geosphere)
library(RColorBrewer)
display.brewer.all()
brewer.pal(12,'Set3')
brewer.pal(8,'Dark2')
brewer.pal(9,'OrRd')
brewer.pal(9,'RdPu')
brewer.pal(9,'Blues')
brewer.pal(9,'GnBu')
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/pca/MDS/AB/A")
taxaID = read.table("A_IN_matrix.mdist.id")
info1 = read.delim("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/map/map_file_r.csv",head =T,sep=",")
index = match(taxaID$V1,info1$ID)
meta = info1[index,c(1,7,8,9)]
gen_dist <- read.table("A_IN_matrix.mdist")
gen_null <- sum(unlist(gen_dist))
col_pal <- matrix(rep(brewer.pal(max_k, "Paired"), max_k),nrow=max_k)
glob <- cmdscale(gen_dist, eig=TRUE,k=3)
xl  <- round(glob$eig[1]/sum(glob$eig[which(glob$eig > 0)]),3)*100
yl <- round(glob$eig[2]/sum(glob$eig[which(glob$eig > 0)]),3)*100
##先画一个空的
plot(-glob$points[,2],-glob$points[,1], xlab=yl, ylab=xl, main="Metric MDS", xaxt = 'n', yaxt = 'n', col=  NULL)
index1 = which(meta$Common_name=="Wild emmer")
#points(-glob$points[index1,2],-glob$points[index1,1], col = "#FA9FB5", cex = 1, pch =8)
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(250/255,159/255,181/255,0.8), cex = 1, pch =8)
index1 = which(meta$Common_name=="Domesticated emmer")
#points(-glob$points[index1,2],-glob$points[index1,1], col = "#BEBADA", cex = 1.5, pch =18)
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(190/255,186/255,218/255,0.7), cex = 1.7, pch =18)
index1 = which(meta$Common_name=="Ispahanicum")
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(190/255,186/255,218/255,0.9), cex = 1.2, pch =9)
index1 = which(meta$Common_name=="Georgian wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(190/255,186/255,218/255,0.9), cex = 1, pch =3)
index1 = which(meta$Common_name=="Rivet wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(179/255,222/255,105/255,0.7), cex = 1.28, pch =17)
index1 = which(meta$Common_name=="Polish wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#B3DE69", cex = 1, pch =11)
index1 = which(meta$Common_name=="Persian wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#B3DE69", cex = 1.2, pch =10)
index1 = which(meta$Common_name=="Khorasan wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#B3DE69", cex = 1, pch =5)
index1 = which(meta$Common_name=="Durum")
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(179/255,222/255,105/255,0.7), cex = 1.28, pch =15)


index1 = which(meta$Common_name=="Landrace")
#points(-glob$points[index1,2],-glob$points[index1,1], col = "#80B1D3", cex = 1, pch =6)
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(127/255,177/255,211/255,0.8), cex = 1, pch =6)
index1 = which(meta$Common_name=="Cultivar")
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(127/255,177/255,211/255,0.8), cex = 1, pch =12)
index1 = which(meta$Common_name=="Spelt")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#80B1D3", cex = 1, pch =0)
index1 = which(meta$Common_name=="Macha")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#80B1D3", cex = 1, pch =7)
index1 = which(meta$Common_name=="Indian_dwarf_wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = rgb(127/255,177/255,211/255,0.7), cex = 1, pch =16)
index1 = which(meta$Common_name=="Yunan_wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#80B1D3", cex = 1, pch =4)
index1 = which(meta$Common_name=="Tibetan_semi_wild")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#80B1D3", cex = 1, pch =1)
index1 = which(meta$Common_name=="Xinjiang_wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#80B1D3", cex = 1, pch =13)
index1 = which(meta$Common_name=="Vavilovii")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#80B1D3", cex = 1, pch =14)
index1 = which(meta$Common_name=="Club_wheat")
points(-glob$points[index1,2],-glob$points[index1,1], col = "#80B1D3", cex = 1, pch =2)

library(gdata)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
display.brewer.all()
scale_fill_brewer(palette = "RdYiBu")
setwd("/Users/guoyafei/Documents/02_VmapIII/04_Statistics/05_MDS")
anno <- read.xls("/Users/guoyafei/Documents/02_VmapIII/01_表格/V.2/Lulab_germplasm_Info.xlsx",sheet=1,na.strings=c("NA","#DIV/0!"))
data <- read.table("ABlineage.maf0.01.mds", header=T, stringsAsFactors = F)
data[1356:1408,1] <- c("ABD_1144","ABD_1145","ABD_1146","ABD_1147","ABD_1148","ABD_1149","ABD_1150","ABD_1151","ABD_1152","ABD_1153","ABD_1154","ABD_1155","ABD_1156","ABD_1157","ABD_1158","ABD_1159","ABD_1160","ABD_1161","ABD_1162","ABD_1163","ABD_1164","ABD_1165","ABD_1166","ABD_1167","ABD_1168","ABD_1169","ABD_1170","ABD_1171","ABD_1172","ABD_1173","ABD_1174","ABD_1175","ABD_1176","ABD_1177","ABD_1178","ABD_1179","ABD_1180","ABD_1181","ABD_1182","ABD_1183","ABD_1184","ABD_1185","ABD_1186","ABD_1187","ABD_1188","ABD_1189","ABD_1190","ABD_1191","ABD_1192","ABD_1193","ABD_1194","ABD_1195","ABD_1196")
test <- merge(data,anno,by.x="FID",by.y="Taxa.vmap3.")

color = scale_color_brewer(brewer.pal(6, "Set2")[c(1,2,3,4,6)])
p <- ggplot(test, aes(C1, C2, color=as.factor(PCA_color) )) +
  geom_point(size=3,shape=as.factor(test$PCA_shape)) +
  scale_color_manual(values = c("#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
  theme(panel.border = element_blank()) +
  guides(fill=guide_legend(title=NULL))

#四倍体
tetra <- test[which(test$Common_name == "Landrace" | test$Common_name == "Tibet_Landrace" |test$Common_name =="Club_wheat"  | test$Common_name == "Indian_dwarf_wheat" | test$Common_name == "Yunan_wheat" | test$Common_name == "Xinjiang_wheat" | test$Common_name == "Vavilovii"),]
tetra <- test[which(test$PCA_color ==2 &test$PCA_color ==1 & test$PCA_color ==3),]
tetra <- test[which(test$PCA_color ==2 & test$add_continent != "NA" ),]
tetra <- test[which(test$add_continent != "NA" & test$Common_name == "Landrace" |test$Common_name == "Tibet_Landrace" | test$Common_name =="Macha"  | test$Common_name == "Spelt" | test$Common_name == "Persian_wheat" | test$Common_name == "Xinjiang_wheat" | test$Common_name == "Polish_wheat" | test$Common_name == "Domesticated_emmer"| test$Common_name == "Rivet_wheat"),]

tetra <- test[which(test$add_continent != "NA" & test$Common_name == "OtherHexaploid" | test$Common_name == "Landrace"),]
tetra <- test[which(test$Common_name != "Cultivar" &test$Common_name != "Cultivar" & test$Common_name != "OtherHexaploid" & test$Common_name != "Synthetichexaploid"),]
tetra <- test[which(test$Common_name == "Landrace" | test$Common_name == "Tibet_Landrace" |test$Common_name =="Club_wheat"  | test$Common_name == "Indian_dwarf_wheat" | test$Common_name == "Yunan_wheat" | test$Common_name == "Xinjiang_wheat" | test$Common_name == "Vavilovii"  | test$Common_name == "Strangulata"),]
tetra <- test[which(test$Common_name == "Landrace" | test$Common_name == "Tibet_Landrace" |test$Common_name =="Club_wheat"  | test$Common_name == "Indian_dwarf_wheat" | test$Common_name == "Yunan_wheat" | test$Common_name == "Xinjiang_wheat" | test$Common_name == "Vavilovii" ),]
tetra <- test[which(test$Common_name == "Cultivar"),]
tetra <- test[which(test$Common_name == "Cultivar" | test$Common_name == "Landrace"| test$Common_name == "Strangulata" ),]

ggplot(tetra, aes(C1, C2, color=as.factor(tetra$add_continent),shape=as.factor(tetra$Common_name))) +
  geom_point(size=3) +
  #scale_color_manual(values = c("#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  scale_color_manual(values = c(brewer.pal(12, "Set3")[c(1:12)]),)+
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11))+
  theme_classic() +
  geom_label_repel(aes(label =tetra$FID), size = 3,  show.legend = FALSE)


ggplot(tetra, aes(C1, C2, color=as.factor(tetra$PCA_color),shape=as.factor(tetra$Common_name))) +
  geom_point(size=3) +
  #scale_color_manual(values = c("#FC8D62","#8DA0CB","#E78AC3","#FFD92F")) +
  scale_color_manual(values = brewer.pal(12, "Set3")[c(1:12)])+
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11))+
  theme_classic() +
  geom_label_repel(aes(label =tetra$FID), size = 3,  show.legend = FALSE)

strang <- read.table("strangulata.txt",header=F,stringsAsFactors = F)
library(maps)
library(ggplot2)
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70") 
mp<-ggplot()+mapworld+ylim(-50,80)
mp_13<-mp+geom_point(aes(x=strang$V4, y=strang$V3,color = as.factor(strang$V6)))+
  scale_size(range=c(1,1))+ 
  theme_classic()


