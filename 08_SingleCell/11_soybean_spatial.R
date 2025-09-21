args = commandArgs(T)

if (length(args) != 6) {
	                stop ("Usage: Rscript stomics_Gmax.ZH13_cc.R <setwd> <gem> <bin> <cut> <pc> <res>\n")
}

library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(ggplot2)
library(ggsci)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(future)
library(harmony)

dir.create(args[1])
setwd(args[1]) ${sample}_bin50 
infile <- args[2] /mount/disk1/guoyafei1/01_soybean_spatial/01_stero/addUSDA110/merge/${sample}/${sample}.merged.gem.gz 
bs <- as.numeric(args[3]) 50
cut <- as.numeric(args[4]) 0
pc <- as.numeric(args[5]) 30
res <- as.numeric(args[6]) 0.5
pro <- args[1]


LD_PRELOAD=/usr/lib64/libstdc++.so.6 LD_LIBRARY_PATH=/home/guoyafei1/software/anaconda3/lib:$LD_LIBRARY_PATH 
Rscript 
/home/guoyafei1/01_soybean_spatial/stomics_${ref}.R 
${sample}_bin50 
/mount/disk1/guoyafei1/01_soybean_spatial/01_stero/addUSDA110/merge/${sample}/${sample}.merged.gem.gz 
50 
0 
30 
0.5

############################## 1. bin data  ##############################
dat <- fread(file = "cd_21.merged.gem.gz")

if(length(grep("MIDCounts|MIDCount",colnames(dat))>0)){
  colnames(dat) <- gsub("MIDCounts|MIDCount","UMICount",colnames(dat))}
out <- as.data.frame(dat)
dat$x <- trunc((dat$x - min(dat$x)) / 50 + 1)
dat$y <- trunc((dat$y - min(dat$y)) / 50 + 1)	
out <- cbind(dat$y,dat$x,out)
colnames(out)[1:2] <- c(paste0("bin",50,".y"),paste0("bin",50,".x"))
#fwrite(out,paste0(pro,"_bin",bs,"_information.txt"),col.names=T,row.names=F,sep="\t",quote=F)
dat <- dat[, sum(UMICount), by = .(geneID, x, y,Slice)]
dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
slice_info <- dat[,c("bin_ID","Slice")]
slice_info <- slice_info %>% distinct(bin_ID, .keep_all = TRUE)
slice_info$bin_ID <- as.integer(slice_info$bin_ID)
bin.coor <- dat[, sum(V1), by = .(x, y)]
geneID <- seq(length(unique(dat$geneID)))
hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
gen <- hash.G[dat$geneID, 'values']
bin_ID <- unique(dat$bin_ID)
hash.B <- data.frame(row.names = sprintf('%d', bin_ID), values = bin_ID)
bin <- hash.B[sprintf('%d', dat$bin_ID), 'values']

cnt <- dat$V1
rm(dat)
gc()
tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))
tissue_positions_list <- data.frame(row.names = paste('BIN', rownames(hash.B), sep = '.'), tissue = 1, row = bin.coor$y, col = bin.coor$x, imagerow = bin.coor$y, imagecol = bin.coor$x)
scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1, tissue_hires_scalef = 1, tissue_lowres_scalef = 1))
mat <- sparseMatrix(i = gen, j = bin, x = cnt)
rownames(mat) <- rownames(hash.G)
colnames(mat) <- paste('BIN', sprintf('%d', seq(max(hash.B[, 'values']))), sep = '.')

############################## 2. creat Spatial Object  ##############################

seurat_spatialObj <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial',min.cells=1, min.features=1)
# Add Slice information to meta.data
slice_info <- slice_info[match(colnames(seurat_spatialObj), paste('BIN', slice_info$bin_ID, sep = '.')), ]
seurat_spatialObj$Slice <- slice_info$Slice
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE)
{
  if (filter.matrix) {tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]}
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius / max(dim(image))
  return(new(Class = 'VisiumV1', image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef), coordinates = tissue.positions, spot.radius = spot.radius))
}
spatialObj <- generate_spatialObj(image = tissue_lowres_image, scale.factors = fromJSON(scalefactors_json), tissue.positions = tissue_positions_list)

spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'
seurat_spatialObj[['slice1']] <- spatialObj

rm(mat)
rm(bin.coor)
rm(hash.G)
rm(hash.B)
rm(bin)
rm(gen)
rm(cnt)

# dir.create("figures")

##############################  3. Spatial Data QC  ##############################
#QC
geneinfo <- as.data.frame(fread("/mount/disk1/guoyafei1/01_soybean_spatial/01_stero/addUSDA110/ChiheiaddUSDA110_gene_info.txt"))
geneinfo$GeneID <- gsub("_", "-", geneinfo$GeneID)

Mtgene <- as.vector(geneinfo[which(geneinfo$Chr=="ChrM"),"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr=="ChrM"),"GeneID"]) %in% rownames(seurat_spatialObj))]
cat("Mitochondrial gene:",Mtgene,"\n")
seurat_spatialObj[["percent.mt"]] <- PercentageFeatureSet(seurat_spatialObj, features = Mtgene)

Cgene <- as.vector(geneinfo[which(geneinfo$Chr=="ChrC"),"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr=="ChrC"),"GeneID"]) %in% rownames(seurat_spatialObj))]
cat("Chloroplast gene:",Cgene,"\n")
seurat_spatialObj[["percent.C"]] <- PercentageFeatureSet(seurat_spatialObj, features = Cgene)

Ugene <- as.vector(geneinfo[which(geneinfo$Chr=="NZ_CP011360.1"),"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr=="NZ_CP011360.1"),"GeneID"]) %in% rownames(seurat_spatialObj))]
cat("USDA110 gene:",Ugene,"\n")
seurat_spatialObj[["percent.U"]] <- PercentageFeatureSet(seurat_spatialObj, features = Ugene)

counts_mat <- GetAssayData(seurat_spatialObj, assay = "Spatial", slot = "counts")
usda_counts <- counts_mat[Ugene, , drop = FALSE]
nFeature_USDA110 <- Matrix::colSums(usda_counts > 0)
seurat_spatialObj$nFeature_USDA110 <- nFeature_USDA110

all_genes <- rownames(seurat_spatialObj)
soybean_genes <- setdiff(all_genes, Ugene)
soybean_counts <- counts_mat[soybean_genes, , drop = FALSE]
nFeature_soybean <- Matrix::colSums(soybean_counts > 0)
seurat_spatialObj$nFeature_soybean <- nFeature_soybean



usda_counts_total <- Matrix::colSums(counts_mat[Ugene, , drop = FALSE])
seurat_spatialObj$UMICount_USDA110 <- usda_counts_total

soybean_counts_total <- Matrix::colSums(counts_mat[soybean_genes, , drop = FALSE])
seurat_spatialObj$UMICount_soybean <- soybean_counts_total

seurat_spatialObj[["percent.S"]] <- PercentageFeatureSet(seurat_spatialObj, features = soybean_genes)

Q1 <- quantile(seurat_spatialObj$nFeature_Spatial)[2]
Q3 <- quantile(seurat_spatialObj$nFeature_Spatial)[4]
upper <- as.numeric(Q3+1.5*(Q3-Q1))
lower <- as.numeric(Q1-1.5*(Q3-Q1))

############ raw Quality
seurat_spatialObj <- subset(seurat_spatialObj,subset = nFeature_Spatial > 0 &  percent.mt < 5 & percent.C < 5)
pdf(paste0("figures/cd_21_bin50","_QC.Feature.pdf"),10,7)
VlnPlot(seurat_spatialObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt","percent.C","percent.U"), ncol=3,pt.size = 0) +theme(axis.text.x=element_text(angle=20,size=0),axis.title.x=element_text(angle=20,size=10))+labs(x=paste0("nGene: ",dim(seurat_spatialObj)[1],"; ","nBIN: ",dim(seurat_spatialObj)[2]))
plot(density(seurat_spatialObj$nFeature_Spatial))+abline(v=c(50,100,150,200,500),col="grey",lwd=2,lty=6)

SpatialFeaturePlot(seurat_spatialObj, features="nFeature_Spatial", stroke=0, pt.size.factor = 1.5) + theme(legend.position = "right") + theme_void()
SpatialFeaturePlot(seurat_spatialObj, features="nCount_Spatial", stroke=0, pt.size.factor = 1.5) + theme(legend.position = "right") + theme_void()
SpatialFeaturePlot(seurat_spatialObj, features="percent.U", stroke=0, pt.size.factor = 1.5) + theme(legend.position = "right") + theme_void()
SpatialFeaturePlot(seurat_spatialObj, features="nFeature_USDA110", stroke=0, pt.size.factor = 1.5) + theme(legend.position = "right") + theme_void()
SpatialFeaturePlot(seurat_spatialObj, features="UMICount_USDA110", stroke=0, pt.size.factor = 1.5) + theme(legend.position = "right") + theme_void()
SpatialFeaturePlot(seurat_spatialObj, features="percent.S", stroke=0, pt.size.factor = 1.5) + theme(legend.position = "right") + theme_void()
SpatialFeaturePlot(seurat_spatialObj, features="nFeature_soybean", stroke=0, pt.size.factor = 1.5) + theme(legend.position = "right") + theme_void()
SpatialFeaturePlot(seurat_spatialObj, features="UMICount_soybean", stroke=0, pt.size.factor = 1.5) + theme(legend.position = "right") + theme_void()
dev.off()

############ QC
pdf(paste0("figures/cd_21_bin50_QC.Feature.pdf"),10,7)

#seurat_spatialObj <- subset(seurat_spatialObj,subset = nFeature_Spatial > cut & percent.mt < 5 & percent.C < 5)
VlnPlot(seurat_spatialObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt","percent.C","percent.U"), ncol=3,pt.size = 0) +theme(axis.text.x=element_text(angle=20,size=0),axis.title.x=element_text(angle=20,size=10))+labs(x=paste0("nGene: ",dim(seurat_spatialObj)[1],"; ","nBIN: ",dim(seurat_spatialObj)[2]))
plot(density(seurat_spatialObj$nFeature_Spatial))+abline(v=c(50,100,150,200,500),col="grey",lwd=2,lty=6)
SpatialFeaturePlot(seurat_spatialObj, features="nFeature_Spatial", stroke=0) + theme(legend.position = "right") + theme_void()
SpatialFeaturePlot(seurat_spatialObj, features="nCount_Spatial", stroke=0) + theme(legend.position = "right") + theme_void()
dev.off()

#######################################  4.Clustering ###################################

#seurat_spatialObj <- SCTransform(seurat_spatialObj, assay = "Spatial", verbose = FALSE, variable.features.n=3000)
seurat_spatialObj <- NormalizeData(object = seurat_spatialObj, verbose = FALSE)
seurat_spatialObj <- FindVariableFeatures(object = seurat_spatialObj, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
seurat_spatialObj <- ScaleData(seurat_spatialObj)
seurat_spatialObj <- RunPCA(seurat_spatialObj,verbose = F)
#seurat_spatialObj <- RunHarmony(seurat_spatialObj, group.by.vars="Slice")
### Find cluster
#seurat_spatialObj <- FindNeighbors(seurat_spatialObj, reduction = "harmony", dims = 1:30)
#seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = FALSE,resolution = 0.5)
#seurat_spatialObj <- RunUMAP(seurat_spatialObj, reduction = "harmony", dims = 1:30)

seurat_spatialObj <- FindNeighbors(seurat_spatialObj, dims = 1:30)
seurat_spatialObj <- FindClusters(seurat_spatialObj, verbose = FALSE,resolution = 0.5)
seurat_spatialObj <- RunUMAP(seurat_spatialObj, dims = 1:30)

col <- c("#F56867","#FEB915","#59BE86","#C798EE","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785","#8DD3C7","#BEBADA","#80B1D3","#B3DE69","#FCCDE5","#BC80BD","#FFED6F","#8DA0CB","#E78AC3","#E5C494","#CCCCCC","#FB9A99","#E31A1C","#CAB2D6","#6A3D9A","#1B9E77","#FCCDE5","#80CDC1","#FEE090","#1F78B4","#7570B3","#FED9A6",'#c49c94',"#BF812D","#BC80BD","#ABD9E9","#C51B7D","#4D9221","#DECBE4","#BEAED4","#7FBC41","#FBB4AE","#FDC086","#B3CDE3","#B3DE69","#8DD3C7","#01665E","#DE77AE","#CCEBC5","#FB8072","#8C510A","#FDAE61","#7FC97F","#BEBADA","#66A61E","#276419","#A6761D","#35978F","#F56867","#FEB915","#59BE86","#C798EE","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370")

#plot2 <- ElbowPlot(seurat_spatialObj, ndims=50, reduction="pca")
pdf(paste0("figures/cd_21_bin50_UMAP.pdf"))
#print(plot2)
DimPlot(seurat_spatialObj, reduction="umap", label=TRUE,cols = col,pt.size=0.5)
SpatialDimPlot(seurat_spatialObj,stroke=0, pt.size.factor=1.2) + theme_void()+scale_fill_manual(values=col)
dev.off()

pdf(paste0("figures/cd_21_bin50_UMAP_split.pdf"), length(levels(seurat_spatialObj))*5, 6)
DimPlot(seurat_spatialObj, reduction="umap", label=TRUE, split.by = "seurat_clusters", cols = col)
dev.off()

############  Cluster Ident

pdf(paste0("figures/cd_21_bin80_Cluster_Ident.pdf"),10,7)
p <- SpatialDimPlot(seurat_spatialObj, cells.highlight=CellsByIdentities(object=seurat_spatialObj),cols.highlight = c("#DE2D26", "grey90"), facet.highlight=TRUE, stroke=0, pt.size.factor=1.2, ncol=4, combine=F)
for(i in 1:length(levels(seurat_spatialObj))) {
  p[[i]] <- p[[i]] +theme_void()
}
p
dev.off()

### save info for customizing SpatialDimPlot
# SpatialPlotInfo <- cbind(seurat_spatialObj@meta.data,seurat_spatialObj@images$slice1@coordinates)
# write.table(SpatialPlotInfo,paste0(pro,"_bin",bs,"_SpatialPlotInfo.txt"),quote=F,sep="\t",row.names=T,col.names=NA)
# info <- cbind(seurat_spatialObj@images$slice1@coordinates[c(2,3)],seurat_spatialObj@meta.data$seurat_clusters)
# colnames(info) <- c('row','col','cluster')
# write.table(info,paste0(pro,"_bin",bs,"_cluster.txt"),quote=F,sep="\t",row.names=F,col.names=T)
# write.table(col,paste0(pro,"_bin",bs,"_col.txt"),quote=F,sep="\t",row.names=F,col.names=F)

# saveRDS(seurat_spatialObj,file=paste0(pro,"_bin",bs,"_Cluster.rds"))

######  DEG of patial.assay
de_markers <- FindAllMarkers(seurat_spatialObj,only.pos = FALSE, assay="Spatial", min.pct = 0.25, logfc.threshold = 0.25)
de_markers <- de_markers[de_markers$p_val_adj<0.05,]
top10 <- de_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)

pdf(paste0("figures/cd_21_bin50_Cluster_top10Heat.pdf"),13,13)
DoHeatmap(seurat_spatialObj, assay="Spatial", features=top10$gene,slot = "data") + NoLegend()
dev.off()
# write.table(de_markers, file=paste0(pro,"_bin",bs,"_DEmarker.txt"),sep = "\t",quote = F)
# write.table(top10, file=paste0(pro,"_bin",bs,"_top10DEmarker.txt"),sep = "\t",quote = F)

q()
# celltype_col <- data.frame()
# for (i in 1:length(unique(seurat_spatialObj$celltype))){
#         n <- data.frame(celltype=levels(as.factor(seurat_spatialObj$celltype))[i],cols = col[i])
#         celltype_col <- rbind(celltype_col,n)
# }
# write.table(celltype_col,paste0(pro,"_bin",bs,"_celltype_col.txt"),quote=F,sep="\t",row.names=F,col.names=F)


######################################################################### cluster visualization ##############################################################

args <- commandArgs (T)
if (length(args) != 4) {
	stop ("Rscript plot.spatial.cluster.R <in.txt> <color.list> <bin> <plot.leg>\n");
}

library (plotrix)

rt <- read.table (args[1], header=T);
rad <- 10
row_max <- max (rt$row);
col_max <- max (rt$col);
cvs_bin <- max (row_max, col_max);
cvs_max <- cvs_bin*2*rad *4/3;
bin <- as.numeric(args[3])

pdf (paste0(args[4],"_bin",bin,"_cluster.pdf"));

w <- cvs_max;
h <- w * 1.085;
plot (c(1,h), c(1,w), cex=0, xlab='', ylab='', main='', axes=F);

row_os <- cvs_max/2 - row_max*2*rad/2;
cols <- read.table (args[2], header=F, comment.char='')$V1;
rt$row <- row_max + 1 - rt$row;
nbins <- nrow (rt)
for (i in 1:nbins) {
	x <- rt$col[i] * 2 * rad;
	y <- rt$row[i] * 2 * rad + row_os;
	draw.circle (x, y, rad, nv=100, border=NA, col=cols[rt$cluster[i]+1]);
}

bottom <- cvs_max/2 - row_max*2*rad/4;
row_os <- cvs_max/60
n_cls <- length (unique(rt$cluster));
if(n_cls > 15){
	leg_dist <- row_max*2*rad/(n_cls-1)
}
leg_dist <- row_max*2*rad/2/(n_cls-1);
for (i in 1:n_cls) {
	cls_name <- i - 1;
	cls_name <- as.character (cls_name);
	x = col_max*2*rad + cvs_max/25;
	y = cvs_max/2 - (i-n_cls/2)*leg_dist + row_os;
	draw.circle (x, y, 10000/bin/2/5, nv=100, border=NA, col=cols[i]);
	text (x+10000/bin/2, y, cls_name, cex=0.5);
	bottom = y;
}

x <- c(col_max*2*rad+cvs_max/25-2*rad, col_max*2*rad+cvs_max/25-2*rad+10000/bin*2)
y <- c(bottom-cvs_max/8, bottom-cvs_max/8)
lines (x, y, lwd=2);
text (col_max*2*rad+cvs_max/25+15*rad, bottom-cvs_max/10+8*rad, "500um", cex=0.6)

dev.off ()





/mount/disk1/guoyafei1/01_soybean_spatial/01_stero/addUSDA110/merge/cd_07/cd_07.merged.gem.gz
/mount/disk1/guoyafei1/01_soybean_spatial/01_stero/addUSDA110/merge/cd_21/cd_21.merged.gem.gz






















