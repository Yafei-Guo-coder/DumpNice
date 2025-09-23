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

projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_globular_anchor_True.rds")
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_heart_anchor_True.rds")
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_torpedo_anchor_True.rds")
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_bent_anchor_True.rds")
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_cotyledon_anchor_True.rds")

idxPass <- which(projAll$predictedScore > 0.5)
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]

cells_back <- read.table("anchor_ED/txt/cells_bent.txt", stringsAsFactors = FALSE, comment.char = "")[,1]
idxPass <- which(!projAll$cellNames %in% cells_back)
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]

idxBack <- which(projAll$cellNames %in% cells_back)

# 把这些细胞的 predictedGroup 改成 "General embryo"
projAll$predictedGroup[idxBack] <- "General embryo"


# projAll$highlight <- ifelse(projAll$cellNames %in% cells_back, "bent_PEN", "other")
# pdf("anchor_ED/pdf/bent_PEN.pdf")
# print(p)
# dev.off()
tab <- table(projAll$Clusters, projAll$predictedGroup)

# 每行找最大值对应的细胞类型、数量和占比
result <- apply(tab, 1, function(x) {
  max_idx <- which.max(x)                     # 最大值列号
  celltype <- colnames(tab)[max_idx]          # 最大值对应的细胞类型
  count <- x[max_idx]                         # 最大值数量
  prop <- x[max_idx] / sum(x)                 # 最大值占比
  c(CellType = celltype, Count = count, Proportion = round(prop, 3))
})

# 转置成数据框
result_df <- as.data.frame(t(result))
result_df

idxPass <- which((projAll$Clusters == "C5" |  projAll$Clusters == "C6") &  projAll$predictedGroup == "Peripheral endosperm")
cellsPass_globular <- projAll$cellNames[idxPass]
idxPass <- which((projAll$Clusters == "C1" |  projAll$Clusters == "C12") &  projAll$predictedGroup == "Chalazal endosperm")
cellsPass_heart <- projAll$cellNames[idxPass]
idxPass <- which((projAll$Clusters == "C2" | projAll$Clusters == "C6" ) & ( projAll$predictedGroup == "Chalazal endosperm" | projAll$predictedGroup == "Peripheral endosperm"))
cellsPass_torpedo <- projAll$cellNames[idxPass]
# 这里考虑把另外一个的细胞也去掉
idxPass <- which((projAll$Clusters == "C5" |  projAll$Clusters == "C6") &  projAll$predictedGroup == "General endosperm")
cellsPass_bent <- projAll$cellNames[idxPass]
idxPass <- which((projAll$Clusters == "C8" ) &  projAll$predictedGroup == "Micropylar endosperm")
cellsPass_cotyledon <- projAll$cellNames[idxPass]
all_endo <- c(cellsPass_globular, cellsPass_heart, cellsPass_torpedo, cellsPass_bent, cellsPass_cotyledon)

write.table(
  all_endo,
  file = "anchor_ED/txt/all_endo.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cells_back <- read.table("anchor_ED/txt/all_endo.txt", stringsAsFactors = FALSE, comment.char = "")[,1]

# 这里读进去endo整合之后的rds文件为projAll
# 使用更新过的snRNA数据的注释进行ATAC数据细胞类型映射
# projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds_Endosperm/ArchR_EndospermD_harmony_rmTorpedo2_motif_V2.rds")
# projAll <- addImputeWeights(projAll)
# # new RNA annotation
# seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
# cell_ids <- seRNA@meta.data$cellID
# cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
# seRNA_filtered <- subset(seRNA, cells = cells_to_keep)
# projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = TRUE, force = TRUE, groupRNA = "endotype2025", nameCell = "predictedCell_endo2025", nameGroup = "predictedGroup_endo2025", nameScore = "predictedScore_endo2025", threads = 1)



saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2.rds")

projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_coNewRNA.rds")
idxPass <- which(projAll$cellNames %in% cells_back)
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]



p <- plotEmbedding(ArchRProj = projAll, colorBy = "cellColData",name = "highlight",embedding = "UMAPHarmony") + scale_color_manual(values = c("other" = "grey80", "bent_PEN" = "red"))




