library(matrixStats)
# packageVersion("matrixStats")
library(ArchR)
library(BSgenome.Athaliana.TAIR10.02022009)
library(parallel)
library(ggalt)
library(clustree)
library(dplyr)
library(patchwork)
library(gridExtra)
library(grid)
library(ggalluvial)
library(ggplot2)
library(SummarizedExperiment)
library(motifmatchr)
library(Matrix)
library(rhdf5)
library(readxl)
library(tidyverse)
library(cluster)
# library(hrbrthemes)
library(viridis)
library(pheatmap)
library(Seurat,lib = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
library(Signac)
#library(EnsDb.Hsapiens.v86)
library(cowplot)
library(parallel)
library(ComplexHeatmap)
library(S4Vectors)
library(chromVAR)
library(TFBSTools)
# library(SeuratData)
# install the dataset and load requirements
# InstallData("pbmcMultiome")
library(SeuratObject,lib= "./lib2")
#.libPaths("/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
# devtools::install_version("matrixStats", version="1.5.0", lib = "./lib")
# install.packages("universalmotif", lib = "./lib")
library(universalmotif)
library(tidydr)

#====================================================================================================================#
#                                     使用Seurat进行RNA和ATAC数据映射                                                ####
#====================================================================================================================#
# 处理RNA数据
pbmc.rna <- readRDS("root_final202502.rds")
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

# 处理atac数据
projHeme2 <- readRDS("rds/PCA_snATAC_motifs_silique.rds")
geneScoreMat <- getMatrixFromProject(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix")
geneScoreMat_sparse <- assay(geneScoreMat)
rownames(geneScoreMat_sparse) <- rowData(geneScoreMat)$name
cells <- colnames(geneScoreMat_sparse)
pbmc.atac <- CreateSeuratObject(counts = geneScoreMat_sparse, assay = "ATAC")

pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)
pbmc.atac <- RunPCA(pbmc.atac, npcs = 30)
# 使用Lsi聚类
pbmc.atac <- RunSVD(pbmc.atac, assay = "ATAC", reduction.key = "LSI_", reduction.name = "lsi")
pbmc.atac <- FindNeighbors(pbmc.atac, reduction = "lsi", dims = 1:30)
pbmc.atac <- FindClusters(pbmc.atac, resolution = 8)
pbmc.atac <- RunUMAP(pbmc.atac,  reduction = "lsi",dims = 1:30)
# 在最后画图

# 取子集计算（需要确认一下细胞的数量）
DefaultAssay(pbmc.rna) <- "RNA"
DefaultAssay(pbmc.atac) <- "ATAC"
pbmc.reference.small <- subset(pbmc.rna, cells = sample(colnames(pbmc.rna), 10000))
pbmc.reference.small <- NormalizeData(pbmc.reference.small)
pbmc.reference.small <- FindVariableFeatures(pbmc.reference.small)
pbmc.reference.small <- ScaleData(pbmc.reference.small)
pbmc.reference.small <- RunPCA(pbmc.reference.small)
pbmc.query.small <- subset(pbmc.atac, cells = sample(colnames(pbmc.atac), 8000))
pbmc.query.small <- RunTFIDF(pbmc.query.small)
pbmc.query.small <- FindTopFeatures(pbmc.query.small, min.cutoff = 'q0')
pbmc.query.small <- RunSVD(pbmc.query.small)

# anchors
transfer.anchors <- FindTransferAnchors(reference = pbmc.reference.small, query = pbmc.query.small, features = VariableFeatures(pbmc.reference.small), reference.assay = "RNA", query.assay = "ATAC", reduction = "cca")

# Transfer labels
predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = pbmc.reference.small$cell_type, weight.reduction = pbmc.query.small[['lsi']], dims = 1:30)
atac <- AddMetaData(pbmc.query.small, metadata = predicted.labels)

# 共嵌入可视化
genes.use <- VariableFeatures(pbmc.reference.small)
refdata <- GetAssayData(pbmc.reference.small, assay = "RNA", slot = "data")[genes.use, ] 
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["pca"]], dims = 2:30) 
atac[["RNA"]] <- imputation 

# 进行merge合并 
coembed <- merge(x = pbmc.reference.small, y = atac) 
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE) 
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE) 
coembed <- RunUMAP(coembed, dims = 1:30) 
coembed$celltype <- ifelse(!is.na(coembed$cell_type), coembed$cell_type, coembed$predicted.id) 
coembed$tech <- ifelse(!is.na(coembed$cell_type), "RNA", "ATAC") 

# 画图
pdf("RNA-ATAC-integration-plot2.pdf",width = 10, height = 6)
DimPlot(pbmc.rna, group.by = "cell_type", label = TRUE) + NoLegend() + ggtitle("RNA")
DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
DimPlot(coembed, group.by = "tech") 
DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE) 
DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() 
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()


