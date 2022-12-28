library(ggplot2)
#library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(cluster)
library(factoextra)
setwd("/Users/guoyafei/Documents/02_VmapIII/20_ppi")
data <- read.table("Network.line.0.5cfms.txt",header=T,stringsAsFactors = F)
colnames(data) <- c("ID","variable","value")

data2 <- data[,c(2,1,3)]
colnames(data2) <- c("ID","variable","value")

all <- rbind(data,data2)
cast <- dcast(all,ID~variable,fun.aggregate = mean,fill=0)
cast2 <- cast[,-1]
rownames(cast2) <- cast[,1]

for (i in c(1:1117)) {
  cast2[i,i] <- 1
}

df = scale(cast2,center = T,scale = T)
# kmeans聚类
mydata <- df
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:100) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
plot(1:100, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#kmeans聚类，标准化后
data2<- data2[,c(1:22)]
km <- kmeans(df, 15,iter.max = 10000) #用于画地图
fviz_cluster(km, data = df,
             #palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid",
             star.plot = T, 
             repel = T,
             ggtheme = theme_minimal()
)
cast2$type <- km$cluster

a <- cast2[,1118,drop=F]
a$ID <- rownames(a)
write.table(a,"network.node.txt", row.names = F, quote=F, sep="\t")

#hcluster
result <- dist(df, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")
#Ward: 最小方差方法旨在寻找紧凑的球形簇的完整的联动方法找到相似集群。

cast2$type <- cutree(result_hc, k=15)
a <- cast2[,1118,drop=F]
a$ID <- rownames(a)
write.table(a,"network.node_hcluster.txt", row.names = F, quote=F, sep="\t")

################################################ pan-protein ppi soft versus hard sweep #############################################
library(ggplot2)
#library(gridExtra)
library(RColorBrewer)
#library(fitdistrplus)
library(ggrepel)

setwd("/Users/guoyafei/Documents/02_VmapIII/20_ppi/h12/hap100/landrace")
name <- c("14-33","19S-proteasome","20S-proteasome","30S-50S-ribosome","30S-ribosome","40S-ribosome","50S-ribosome","60S-ribosome","90S-preribosome","AHL-dna-binding","AP1-2","ATP-synthase","Argonaute","BRISC","CCT","CHIB-OSM34","CLPP","COG","COP9","CPN10-CPN20","CSP41","CYP18-22","Carbamoyl-phosphate-synthase","Chaperonin-60","Chromatin-+-pre-mRNA","Chromatin-remodelling","Coatomer","Complex-I-III","Cytochrome-b6-f","Dynamin-related","EF1","EIF2-+-4E","EIF2B","EIF3","EIFISO4E-G","ETF","Exon-Junction-Complex","FACT","FTSH1-2","FY-CPSF","GAPDH-PRK","GARP","Glu-1-P-adenylyltransferase","Glucosidase-II","Glutamyl-tRNA-amidotransferase","HOP2-MND1","HOS1-NUP85-107-144-160","HSP70-90","HSP70","Histone-2a","KELP-KIWI","KU70-80","LSM","Lish","Lumenal-NDH","MCC","MCM-large","MCM-small","MOS4-associated","Mediator-14-16","Mediator","Mitochondrial-carriers","Mitochondrial-processing-peptidase","NACA","NADH-flavoprotein","NEDD8-activating-enzyme-E1","Nuclear-pore","OST-EMC","Oxalate-synthesis","PAF1C","PIP-NUDT3","POLD","PP2A","PSBS-CP24","PTAC-rpo","PeBOW","Phenylalanine--tRNA-ligase","Photosystem-I-reaction-center","Photosystem","Proteasome","Pyrophosphate--fructose-6-phosphate-1-phosphotransferase","Pyruvate-dehydrogenase-E1","RAD1-HUS1","RFC","RNA-polymerase-II,-IV-and-V","RNA-polymerase-III","RNA-polymerase","RPN13-UCH1-2","RZ1B-C-VRN1","RuBisCo","SAE1-2","SAGA","SDH","SEC13-31","SEC23-24","SMC","SR-pre-mRNA-splicing","SRP68-72","SWI-SNF-","SWR1","Serine-arginine-rich-splicing-factor","Serrate-cap-binding-complex","SnoRNP","Splicing-factor-3A-B","THO-part-1","THO-part-2","TIM","TRAPP","TSET","VHA","VPS29-35-26","exosome","naa10-15","nap1-npr1","prefoldin","pwp2-utp","pyruvate-dehydrogenase","small-nuclear-ribonucleoprotein-complex","tRNA-multisynthetase")
landrace_name <- paste("landrace.",name,".module.AB.h12",sep="")
cultivar_name <- paste("cultivar.",name,".module.AB.h12",sep="")
landrace_out <- paste("landrace.",name,".pdf",sep="")
cultivar_out <- paste("cultivar.",name,".pdf",sep="")
#thresh
#hap75
cultivar_t <- 0.06642
landrace_t <- 0.06306

#hap100
cultivar_t <- 0.041964165079476
landrace_t <- 0.033131204459876

#cultivar
all1 <- read.table("shuf5k.landrace.AB.h12", header=F,stringsAsFactors = F)
for (i in c(1:length(landrace_name))) {
  gene1 <- read.table(landrace_name[i],header=F,stringsAsFactors = F)
  p <- ggplot(all1, aes(V2, V3)) +
    geom_point(size=2, alpha = 0.5,colour = "grey",shape = 20) +
    geom_point(data = gene1, aes(V2, V3, size=2, alpha = 0.5,colour = "red")) +
    geom_vline(xintercept=landrace_t,color='red',linetype = "dotted")+
    geom_hline(yintercept=0.059,color='red',linetype = "dotted")+
    #geom_text_repel(data = sub_gene1,aes(V2, V3, label = name),max.overlaps = 100) +
    xlab("H12")+
    ylab("H2/H1")+
    theme_bw()+
    theme(legend.position="none",legend.title=element_blank(),axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15))
  pdf(landrace_out[i],height = 6,width = 6)
  print(p)
  dev.off()
  
}


