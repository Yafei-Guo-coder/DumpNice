#env & trait heatmap
setwd("/Users/guoyafei/Desktop/表型/trait/blue/")
cul <- read.table("/Users/guoyafei/Desktop/表型/group/cultivar.txt",header=F,stringsAsFactors = F)
land <- read.table("/Users/guoyafei/Desktop/表型/group/landrace.txt",header=F,stringsAsFactors = F)

#env.trait.cor3_land.txt
env <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)
#envTraitCor_land.txt
env <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/35_variable.txt", header=T, row.names = 1,stringsAsFactors =F)
data <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/flowering.blue.txt", header=F, row.names = 1,stringsAsFactors = F)
sub1 <- env[which(row.names(env) %in% row.names(data)),]
sub2 <- data[which(row.names(data) %in% row.names(env)),,drop=F]
sub3 <- sub1[rownames(sub2),]
#sub3 <- sub1[which(row.names(sub1) %in% land$V1),]
#sub4 <- sub2[which(row.names(sub2) %in% land$V1),,drop=F]

a <- vector()
b <- vector()

for ( i in c(1,2,3,6:40)) {
  p <- cor.test(sub3[,i],sub2$V2)
  a <- append(a,p$p.value)
  b <- append(b,p$estimate)
}
c <- cbind(a,as.numeric(b))
rownames(c) <- colnames(sub1)[c(1,2,3,6:40)]

colnames(c) <- c("p-green-yellow","cor-green-yellow")
all <- cbind(all,c)
colnames(all) <- c("p-heading","cor-heading","p-flowering","cor-flowering","p-green","cor-green","p-yellow","cor-yellow","p-heading-flowering","cor-heading-flowering","p-flowering-green","cor-flowering-green","p-green-yellow","cor-green-yellow")
write.table(all,"env.trait.cor3_land.txt", quote=F,row.names = T,sep="\t")

library(pheatmap)
out <- all[,c(2,4,6,8,10,12,14)]
out <- all[,c(2,4,6,8,12)]
colnames(out) <- c("heading","flowering","green","yellow","flowering-green")
p <- pheatmap(out, show_colnames=FALSE, cutree_cols=3, cluster_row = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors)
p <- pheatmap(out, show_colnames=T, cluster_row = T)
p

#环境相关性分析
library(corrgram)
#env.trait.cor3_land.txt
env <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)
data <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/heading.blue.txt", header=F, row.names = 1,stringsAsFactors = F)
sub1 <- env[which(row.names(env) %in% row.names(data)),]
sub2 <- sub1[rownames(sub2),]
sub <- sub2[,c(1,2,3,6:40)]
#order <- c("prec8","prec6","soil13","prec3","soil12","soil9","prec7","soil2","soil10","temp7","prec2","prec5","prec1","temp4","Elevation","soil1","soil7","Latitude","soil3","temp8","Longitude","soil8","temp6","temp1","temp11","solar2","solar3","soil6","soil4","soil5","temp9","temp10","temp5","soil11","temp3","prec4","solar1","temp2")
order <- c("prec7","solar1","soil8","temp8","solar2","soil7","soil5","temp9","prec8","temp4","soil4","prec4","temp5","temp2","temp10","Elevation","temp7","solar3","soil11","temp3","temp6","temp1","prec3","soil10","prec6","prec5","temp11","soil9","soil6","soil12","soil13","soil2","prec2","prec1","soil3","soil1")

sub3 <- sub[,order]
corrgram(sub3, order = F,lower.panel = NULL,upper.panel=panel.pie,text.panel = panel.txt, gap = 0.1,
         main="Correlogram of environment variables intercorrelations")

#env.trait.cor.txt
env <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/35_variable.txt", header=T, row.names = 1, stringsAsFactors = F)
data <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/heading.blue.txt", header=F, row.names = 1,stringsAsFactors = F)
sub1 <- env[which(row.names(env) %in% row.names(data)),]
sub2 <- sub1[rownames(sub2),]
order <- c("solar2","temp6","solar3","temp1","temp11","soil3","temp8","soil4","soil5","soil8","soil6","soil7","soil12","soil2","prec4","soil11","temp10","temp5","solar1","temp2","temp3","temp9","soil1","soil13","temp7","prec7","soil10","prec2","prec5","prec3","prec6","prec8","soil9","prec1","temp4")
sub3 <- sub2[,order]
corrgram(sub3, order = F,lower.panel = panel.pie,upper.panel=NULL,text.panel = panel.txt, gap = 0.1,
         main="Correlogram of environment variables intercorrelations")

#环境距离和表型距离的相关性
#env.trait.cor3_land.txt
library(tidyr)
library(reshape2)
env <- read.table("/Users/guoyafei/Desktop/bio-RDA/471_baypass_taxa.Info", header=T, row.names = 1, stringsAsFactors = F)
data <- read.table("/Users/guoyafei/Desktop/表型/trait/blue/green-yellow.blue.txt", header=F, row.names = 1,stringsAsFactors = F)
sub1 <- env[which(row.names(env) %in% row.names(data)),]
sub2 <- data[which(row.names(data) %in% row.names(env)),,drop=F]
sub3 <- sub1[rownames(sub2),c(1,2,3,6:40)]
#sub <- sub2[,c(1,2,3,6:40)]
env_sub <- sub3[,c("prec8","prec6","soil13","soil12","soil9","prec3","temp6","temp1","temp11","solar2","solar3")]

df = scale(env_sub,center = T,scale = T)
#df<-env_sub
colnames(df) <- colnames(env_sub)
result1 <- as.data.frame(as.matrix(dist(df, method = "euclidean")))
result1$id <- rownames(result1)
a<- melt(result1,id="id")

df2 = scale(sub2,center = T,scale = T)
#df2 <- sub2
result2 <- as.data.frame(as.matrix(dist(df2, method = "euclidean")))
result2$id <- rownames(result2)
b<- melt(result2,id="id")

geo_sub <- sub3[,c("Longitude","Latitude")]
#df3 = scale(geo_sub,center = T,scale = T)
df3 <- geo_sub 
colnames(df3) <- colnames(geo_sub)
library ("geosphere")
result3 <- as.data.frame(as.matrix(distm(geo_sub, fun = distGeo)))
result3 <- as.data.frame(as.matrix(dist(geo_sub, method = "euclidean")))
result3$id <- rownames(result3)
c<- melt(result3,id="id")

env_all <- sub3[,1:38]
df4 = scale(env_all,center = T,scale = T)
#df4<-env_all
colnames(df4) <- colnames(env_all)
result4 <- as.data.frame(as.matrix(dist(df4, method = "euclidean")))
result4$id <- rownames(result4)
d<- melt(result4,id="id")

all <- cbind(a,b,c,d)
colnames(all) <- c("id1","id2","env_sub_dist","id3","id4","trait_dist_green-yellow","id1","id2","geo_dist","id3","id4","env_all_dist")
write.table(all, "env_trait_dist_green-yellow.txt", quote=F,row.names = F )

library(ggplot2)
env_sub <- vector() 
geom <- vector()
env_all <- vector()
file <- c("env_trait_dist_heading.cut.txt","env_trait_dist_flowering.cut.txt","env_trait_dist_green.cut.txt","env_trait_dist_yellow.cut.txt","env_trait_dist_heading-flowering.cut.txt","env_trait_dist_flowering-green.cut.txt","env_trait_dist_green-yellow.cut.txt")
out1 <- c("env_trait_dist_heading.env.pdf","env_trait_dist_flowering.env.pdf","env_trait_dist_green.env.pdf","env_trait_dist_yellow.env.pdf","env_trait_dist_heading-flowering.env.pdf","env_trait_dist_flowering-green.env.pdf","env_trait_dist_green-yellow.env.pdf")
out2 <- c("env_trait_dist_heading.geo.pdf","env_trait_dist_flowering.geo.pdf","env_trait_dist_green.geo.pdf","env_trait_dist_yellow.geo.pdf","env_trait_dist_heading-flowering.geo.pdf","env_trait_dist_flowering-green.geo.pdf","env_trait_dist_green-yellow.geo.pdf")

for ( i in c(1:length(file))) {
  data <- read.table(file[i], header=T,stringsAsFactors = F)
  colnames(data) <- c("id1","id2","env_sub_dist","trait_dist","geo_dist","env_all_dist")
  a <- cor.test(data$env_sub_dist,data$trait_dist)
  b <- cor.test(data$geo_dist,data$trait_dist)
  c <- cor.test(data$env_all_dist,data$trait_dist)
  env_sub <- append(env_sub,a$p.value)
  env_sub <- append(env_sub,a$estimate)
  geom <- append(geom,b$p.value)
  geom <- append(geom,b$estimate)
  env_all <- append(env_all,c$p.value)
  env_all <- append(env_all,c$estimate)
  pdf(out1[i],width=4.04,height=3.77)
  p <- ggplot(data, aes(x=env_sub_dist, y=trait_dist)) +
    geom_point(size=0.5,alpha=0.3,color="#999999")+
    theme_bw()+
    geom_smooth(method = "lm", se=TRUE, 
                color="#0072B2", formula = y ~ x) 
  print(p)
  dev.off()
  pdf(out2[i],width=4.04,height=3.77)
  q <- ggplot(data, aes(geo_dist, trait_dist)) +
    geom_point(size=0.5,alpha=0.3,color="#999999")+
    theme_bw()+
    geom_smooth(method = "lm", se=TRUE, 
                color="#E69F00", formula = y ~ x)
  print(q)
  dev.off()
}
