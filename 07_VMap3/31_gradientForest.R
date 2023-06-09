library("gradientForest")
gfData <- read.table("/data2/yafei/polygenic/genotype/gradientforest/GF.input_35bio.txt",header = T,sep="\t",row.names = "pop")
colnames(gfData)[3:38] <- paste("bio_",1:36,sep="")
candidate <- gfData[,grep("cand",names(gfData))] 
reference <- gfData [,grep("ref",names(gfData))]
present <- gfData[,c(1,2,grep("bio",names(gfData)))] 
bioclimatic <- paste("bio_",1:36,sep = "") 
maxLevel <- log2(0.368*nrow(candidate)/2) 
gf_candidate <- gradientForest(cbind(present[,bioclimatic], candidate), predictor.vars=colnames(present[,bioclimatic]), response.vars=colnames(candidate), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
gf_reference <- gradientForest(cbind(present[,bioclimatic], reference), predictor.vars=colnames(present[,bioclimatic]), response.vars=colnames(reference), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
gf_runs <- list(gf_reference=gf_reference, gf_candidate=gf_candidate)
save(gf_runs,file = "gf_runs.R") 

load("gf_runs.R") 
present <- gfData[,c(1,2,grep("bio",names(gfData)))] 
gf_candidate <- gf_runs$gf_candidate 
gf_reference <- gf_runs$gf_reference # Once the gradient forest model has been constructed, the importance of each variables
bio_cand <- gf_candidate$overall.imp[order(gf_candidate$overall.imp,decreasing = T)] 
most_cand <- names(bio_cand[1])

pdf("import1_35bio.pdf")
barplot(bio_cand,las=2,cex.names=0.8,col=c(MaizePal::maize_pal("JimmyRed",4),rep("gr ey",15)),ylab="Weigthed importance (R-sqr)")
dev.off()

temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) 
temp_ref_overall <- cumimp(gf_reference,predictor = most_cand, type=c("Overall"),standardize = T) 
temp_ref_SNP <- cumimp(gf_reference,predictor = most_cand, type=c("Species"),standardize = T) 
#ylim <- NULL 
#for(j in 1:length(temp_cand_SNP)){ 
#	ylim <- c(ylim,max(temp_cand_SNP[[j]][[2]]))
#	}
#for(j in 1:length(temp_ref_SNP)){
#	ylim <- c(ylim,max(temp_ref_SNP[[j]][[2]]))
#	} 
#ylim <- max(ylim)
ylim <- 0.2
pdf("import2_35bio.pdf",height=4,width=8)
par(mfrow=c(1,2)) 
par(mai=c(0.9,0.8,0.4,0)) 
plot(temp_cand_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0), ylab="Cumulative importance",xlab= "bio_9")
for(j in 1:length(temp_cand_SNP)){ lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6))} 
lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
par(mai=c(0.9,0.1,0.4,0.6),tcl=-0.2) 
plot(temp_ref_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0),ylab="",xlab= "bio_9", yaxt="n") 
for(j in 1:length(temp_ref_SNP)){ 
	lines(temp_ref_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("MaizAzul")[3],alpha.f = 0.6))
	} 
lines(temp_ref_overall,col=MaizePal::maize_pal("MaizAzul")[6],lwd=4)
dev.off()

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm]) 
pdf("import3_35bio.pdf")
plot(temp_cand_overall,type="n",ylim=c(0,ylim),mgp=c(2,0.6,0), ylab="Cumulative importance",xlab= paste("Most important variable
(",most_cand,")",sep=""),main="Candidate SNPs") 
for(j in 1:length(temp_cand_SNP)){ lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 
lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4)
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 
id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 
id_w <- order(temp$bio[warm]) 
id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 
points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg=rev(id_wcol),cex=1.5) 
points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg=id_ccol,cex=1.5)
dev.off()


#########################
#$cold
#[1] "pop1"  "pop11" "pop12" "pop13" "pop14" "pop15" "pop16" "pop17" "pop18"
#[10] "pop19" "pop20" "pop22" "pop24" "pop25" "pop3"  "pop5"  "pop8"  "pop9" 

#$warm
#[1] "pop10" "pop2"  "pop21" "pop23" "pop4"  "pop6"  "pop7" 
library(readxl)
library(ggmap)
library(RColorBrewer)
library(raster)
library(rgdal)
library(rasterVis)
library(ggplot2)
library(viridis) 

setwd("/Users/guoyafei/Desktop/gradientforest/")
data <- read.table("/data2/yafei/polygenic/genotype/gradientforest/471_baypass_10PC_edited.txt", header=T, stringsAsFactors = F)
pop <- read.table("/data2/yafei/polygenic/genotype/gradientforest/25pop_35bio.env",header=T,stringsAsFactors = F)
wetpop <- c("pop2","pop4","pop6","pop7","pop10","pop21","pop23")
wet <- pop[which(pop$pop %in% wetpop),]
prec7 <- raster("/data2/yafei/polygenic/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_18.tif")
pdf("test.pdf")
plot(prec7)
points(pop$X, pop$Y,pch = 17, bg = "cyan",cex=0.5)
dev.off()

#绘图:全部的
mp <- NULL
mapworld <- borders("world",colour = "gray70",fill="white") 
mp <- ggplot() + 
  mapworld + 
  ylim(-60,90) + 
  theme_classic()
color <- brewer.pal(8, "Dark2")[c(1,2,3,4,5,6,7)]
color <- brewer.pal(8, "Dark2")[c(6)]
#AABBDD
mp+
  geom_point(aes(x=wet$X, y=wet$Y, color=wet$pop),size=0.5) +
  scale_size(range=c(1,1)) + 
  scale_color_manual(values = color) +
  #geom_point(aes(x=pop$X, y=pop$Y, color="red"), size=1) +
  #scale_size(range=c(1,1)) + 
  #scale_color_manual(values = "red") +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "null",
        panel.border = element_blank())


