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

# projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_filter_motifV2meme_UMAPTest.rds")

#====================================================================================================================#
#                                    去掉bent时期C1和C5中的以及cotyledon中聚类异常的细胞                                   #
#====================================================================================================================#
# # Version 1
# # projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_bent_anchor_True.rds")
# # idxPass <- which(projAll$predictedScore > 0.5)
# # cellsPass <- projAll$cellNames[idxPass]
# # projAll <- projAll[cellsPass, ]
# # idxPass <- which(projAll$Clusters == "C1" | projAll$Clusters == "C5")
# # cells_back <- projAll$cellNames[idxPass]
# # writeLines(cells_back, "txt/C1C5_bent.txt")
# cells_back <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/C1C5_bent.txt", stringsAsFactors = FALSE, comment.char = "")[,1]
# projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_coNewRNA.rds")
# idxPass <- which(!projAll$cellNames %in% cells_back)
# cellsPass <- projAll$cellNames[idxPass]
# projAll <- projAll[cellsPass, ]
# idxPass <- which(projAll$predictedScore_Co > 0.5)
# cellsPass <- projAll$cellNames[idxPass]
# projAll <- projAll[cellsPass, ]
# # 8608 cells

# Version 2
# projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_bent_anchor_True.rds")
# idxPass <- which(projAll$predictedScore > 0.5)
# cellsPass <- projAll$cellNames[idxPass]
# projAll <- projAll[cellsPass, ]
# idxPass <- which(projAll$Clusters == "C1" | projAll$Clusters == "C5")
# cells_back <- projAll$cellNames[idxPass]
# writeLines(cells_back, "txt/C1C5_bent.txt")
cells_back1 <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/bent_other1.txt", stringsAsFactors = FALSE, comment.char = "")[,1]
cells_back2 <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/bent_other2.txt", stringsAsFactors = FALSE, comment.char = "")[,1]
cells_back3 <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/cotyloden_other1.txt", stringsAsFactors = FALSE, comment.char = "")[,1]
cells_back4 <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/cotyloden_other2.txt", stringsAsFactors = FALSE, comment.char = "")[,1]

projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_coNewRNA.rds")
# 13189 cells
idxPass <- which(!projAll$cellNames %in% c(cells_back,cells_back2,cells_back3,cells_back4))
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]
# 8786 cells
#====================================================================================================================#
#                                                        无约束整合                                                    #
#====================================================================================================================#
projAll <- addImputeWeights(projAll)
# new RNA annotation
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
cell_ids <- seRNA@meta.data$cellID
cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
seRNA_filtered <- subset(seRNA, cells = cells_to_keep)
# new
projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA_filtered, addToArrow = TRUE, force = TRUE, groupRNA = "endotype2025", nameCell = "predictedCell_endo2025", nameGroup = "predictedGroup_endo2025", nameScore = "predictedScore_endo2025", threads = 1)
saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_filter.rds")
# 8786 cells

#====================================================================================================================#
#                                                         构建伪时序                                                   #
#====================================================================================================================#
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_filter.rds")
projAll <- addUMAP(projAll, reducedDims = "IterativeLSI",minDist = 0.6, name = "UMAPTest", nNeighbors = 30, force = TRUE)
# 8786 cells V2
idxPass <- which(projAll$predictedScore_endo2025 > 0.5)
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]
# 7835 cells V3 (使用这版的数据)

# plot 看预测结果
# projAll <- addIterativeLSI(projAll, useMatrix = "PeakMatrix", name = "IterativeLSI", iterations = 2,force = TRUE)
# #projAll <- addHarmony(projAll, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample")
# projAll <- addClusters(projAll, reducedDims = "IterativeLSI", name = "Clusters", resolution = 0.8,force = TRUE)

# 貌似有一堆异常的细胞，看一下在哪儿
# cells <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/integration_visualization/top_right_atac_ids.txt", header=F, sep="\t", stringsAsFactors = F, quote="",check.names = FALSE,comment.char = "")
# projAll <- projAll[cells[,1], ]

p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAPTest")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAPTest")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup_endo2025", embedding = "UMAPTest")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedScore_endo2025", embedding = "UMAPTest")
# p5 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
# cells_back <- read.table("cells_bent_PEN.txt", stringsAsFactors = FALSE, comment.char = "")[,1]
# projAll$highlight <- ifelse(projAll$cellNames %in% cells_back, "bent_PEN", "other")
# p6 <- plotEmbedding(projAll, colorBy = "cellColData", name = "highlight", embedding = "UMAP")
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure1_V3_UMAPtest.pdf", width = 20, height = 6)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "UMAP - Endosperm Development")
#grid.arrange(p1, p2, ncol = 2, top = "UMAP - Endosperm Development")
dev.off()

# 假设经验上知道 trajectory 顺序是以下几个 cluster
trajectory_order <- c("globular", "heart", "torpedo", "bent", "cotyledon")
# projAll <- addTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", groupBy = "Stage", trajectory = trajectory_order,  preFilterQuantile = 0.5, postFilterQuantile = 0.8, embedding = "UMAP", force = TRUE)
projAll <- addTrajectory(ArchRProj = projAll, name = "EndospermTrajectoryTest", groupBy = "Stage", trajectory = trajectory_order,  preFilterQuantile = 0.9,  embedding = "UMAPTest", force = TRUE)
# saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_coNewRNA_trajectory.rds")

pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure2_V3_UMAPtest.pdf", width = 6, height = 6)
p <- plotTrajectory(projAll, embedding = "UMAPTest", trajectory = "EndospermTrajectoryTest", colorBy = "cellColData", name = "EndospermTrajectoryTest")
print(p)
dev.off()

# 结合伪时间和真实细胞
cell_df <- as.data.frame(getCellColData(projAll))
# 假设你的数据框是 cell_df，包含 pseudotime 和 Stage 两列
# 1. 将 pseudotime 分成 50 个bin
cell_df <- cell_df %>% mutate(bin = cut(EndospermTrajectoryTest, breaks = 50, include.lowest = TRUE))
# 2. 统计每个 bin 内不同 Stage 细胞数量
bin_stage_counts <- cell_df %>% group_by(bin, Stage) %>% summarise(count = n(), .groups = 'drop')
# 3. 计算每个 bin 细胞总数
bin_totals <- bin_stage_counts %>% group_by(bin) %>% summarise(total = sum(count), .groups = 'drop')
# 4. 合并并计算比例
bin_stage_props <- bin_stage_counts %>% left_join(bin_totals, by = "bin") %>% mutate(prop = count / total)
# 5. 画堆积柱状图显示比例
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure3_V3_UMAPtest.pdf",width = 8,height = 3)
ggplot(bin_stage_props, aes(x = bin, y = prop, fill = Stage)) + geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "Stage") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

cell_df$predictedGroup_add <- paste(cell_df$predictedGroup_endo2025,cell_df$Stage,sep="_")
# cell_df$predictedGroup_add <- paste(cell_df$predictedGroup_endo,cell_df$Stage,sep="_")
# 假设伪时间列是 pseudotime，先做bin
cell_df <- cell_df %>% mutate(bin = cut(EndospermTrajectoryTest, breaks = 50, include.lowest = TRUE))
# 统计每个bin中不同predictedGroup_add的数量
bin_group_counts <- cell_df %>% group_by(bin, predictedGroup_add) %>% summarise(count = n(), .groups = 'drop')
# 计算每个bin总数
bin_totals <- bin_group_counts %>% group_by(bin) %>% summarise(total = sum(count), .groups = 'drop')
# 计算比例
bin_group_props <- bin_group_counts %>% left_join(bin_totals, by = "bin") %>% mutate(prop = count / total)
# 画堆积百分比柱状图
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure4_V3_UMAPtest.pdf",width = 8,height = 3)
ggplot(bin_group_props, aes(x = bin, y = prop, fill = predictedGroup_add)) + geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "PredictedGroup_Stage") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# 假设伪时间列是 pseudotime，先做bin
cell_df <- cell_df %>% mutate(bin = cut(EndospermTrajectoryTest, breaks = 50, include.lowest = TRUE))
# 统计每个bin中不同predictedGroup_add的数量
bin_group_counts <- cell_df %>% group_by(bin, predictedGroup_endo2025) %>% summarise(count = n(), .groups = 'drop')
# 计算每个bin总数
bin_totals <- bin_group_counts %>% group_by(bin) %>% summarise(total = sum(count), .groups = 'drop')
# 计算比例
bin_group_props <- bin_group_counts %>% left_join(bin_totals, by = "bin") %>% mutate(prop = count / total)
# 画堆积百分比柱状图
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure5_V3_UMAPtest.pdf",width = 8,height = 3)
ggplot(bin_group_props, aes(x = bin, y = prop, fill = predictedGroup_endo2025)) + geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "PredictedGroup_Stage") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

################################ plot trajectory heatmap - use new motif matrix based on pseudotime trajectory #################################
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", force=TRUE, threads = 1)
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "Stage", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addPeakMatrix(projAll, force = TRUE)
# 添加 motif
meme_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/motif/Ath_TF_binding_motifs_addV2.meme" 
motifs <- read_meme(meme_file)
# 转换为 PWMatrix 格式
corrected_pwms <- lapply(motifs, function(m) {
  matrix_data <- as.matrix(m@motif) 
  new("PWMatrix", ID = m@name, name = ifelse(length(m@altname) > 0, m@altname, m@name), profileMatrix = matrix_data, strand = "+")})
# colSums(corrected_pwms[[1]]@profileMatrix)
# profileMatrix 不是概率矩阵，需要进行归一化
corrected_pwms <- lapply(corrected_pwms, function(m) {
  m@profileMatrix <- sweep(m@profileMatrix, 2, colSums(m@profileMatrix), "/")
  return(m)})
pwmlist <- do.call(PWMatrixList, corrected_pwms)
ids <- sapply(pwmlist, function(x) x@ID)
names(pwmlist) <- ids
# peaks <- getPeakSet(projAll)
# motif_hits <- matchMotifs(pwmlist, peaks, genome = BSgenome.Athaliana.TAIR10.02022009)
# motif_matrix <- assay(motif_hits, "motifMatches")
# overlapMotifs <- motifmatchr::matchMotifs(pwmlist, peaks, genome = "BSgenome.Athaliana.TAIR10.02022009")
# colnames(overlapMotifs) <- paste0("Motif_", seq_len(ncol(overlapMotifs)))
projAll <- addMotifAnnotations(ArchRProj = projAll, force = TRUE, motifPWMs = pwmlist, name = "Motif")
# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll,force = TRUE)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)
# saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_filter_motifV2meme_UMAPTest.rds")
# projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_filter_motifV2meme_UMAPTest.rds")

#====================================================================================================================#
#                                          正向TF-调控因子鉴定                                                       ####
#====================================================================================================================#
# 第一步: 鉴定偏离TF motif
seGroupMotif <- getGroupSE(ArchRProj = projAll, useMatrix = "MotifMatrix", groupBy = "Stage")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){rowMaxs(assay(seZ) - assay(seZ)[,x])}) %>% Reduce("cbind", .) %>% rowMaxs
# 第二步: 鉴定相关的TF motif和TF基因得分/表达值
corGSM_MM <- correlateMatrices(ArchRProj = projAll, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "MotifMatrix", reducedDims = "IterativeLSI")
corGIM_MM <- correlateMatrices(ArchRProj = projAll, useMatrix1 = "GeneIntegrationMatrix", useMatrix2 = "MotifMatrix", reducedDims = "IterativeLSI")
# 第三步: 在相关性DataFrame中添加极大偏差值
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
# 第四步: 鉴定正向TF调控因子并画图
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "Unclear"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.3 & corGSM_MM$padj < 0.001  & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.9))] <- "Up"
corGSM_MM$TFRegulator[which(corGSM_MM$cor < -0.3 & corGSM_MM$padj < 0.001 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.9))] <- "Down"

# corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="Up",1])
sort(corGSM_MM[corGSM_MM$TFRegulator=="Down",1])
library(ggrepel)
p_GSM <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) + 
    geom_point() + 
    theme_ArchR() + 
    geom_vline(xintercept = 0, lty = "dashed") + 
    scale_color_manual(values = c("Unclear"="darkgrey", "Up"="firebrick3","Down"="royalblue3")) + xlab("Correlation To Gene Score") + ylab("Max TF Motif Delta") + scale_y_continuous(expand = c(0,0), limits = c(0, max(corGSM_MM$maxDelta)*1.05)) + geom_text_repel(
    data = subset(corGSM_MM, TFRegulator %in% c("Up","Down")),
    aes(label = GeneScoreMatrix_name),
    color = "black",  # Label color
    size = 3,         # Text size
    max.overlaps = Inf
  )
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "Unclear"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.3 & corGIM_MM$padj < 0.001 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.9))] <- "Up"
corGIM_MM$TFRegulator[which(corGIM_MM$cor < -0.3 & corGIM_MM$padj < 0.001 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.9))] <- "Down"

sort(corGIM_MM[corGIM_MM$TFRegulator=="Up",1])
sort(corGIM_MM[corGIM_MM$TFRegulator=="Down",1])
p_GIM <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("Unclear"="darkgrey", "Up"="firebrick3","Down"="royalblue3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(corGIM_MM$maxDelta) * 1.05)) +
  # Add text labels only for points where TFRegulator == "YES"
  geom_text_repel(
    data = subset(corGIM_MM, TFRegulator  %in% c("Up","Down")),
    aes(label = GeneIntegrationMatrix_name),
    color = "black",  # Label color
    size = 3,         # Text size
    max.overlaps = Inf
  )
  
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure6_V3_V2meme_UMAPtest.pdf",width = 7,height = 7)
print(p_GSM)
print(p_GIM)
dev.off()



#====================================================================================================================#
#                                 combine motifMatrix and geneInteration matrix                                   ####
#====================================================================================================================#
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_filter_motifV2meme_UMAPTest.rds")
# combine motifMatrix and geneInteration matrix
corGSM_MM <- correlateTrajectories(trajGSM, trajMM,corCutOff = 0.3)
idxToRemove <- grep(pattern = "deviations", x = corGSM_MM[["correlatedMappings"]]$name2)
if(length(idxToRemove > 1)){
  corGSM_MM[["correlatedMappings"]] <- corGSM_MM[["correlatedMappings"]][-idxToRemove,]
}

# corGSM_MM[["correlatedMappings"]]
trajGSM2 <- trajGSM[corGSM_MM[["correlatedMappings"]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[["correlatedMappings"]]$name2, ]

trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)

# ComplexHeatmap::draw(ht1 + ht2)
pdf("EndospermTrajectory_umap_reanchor_plot7.pdf", width = 6, height = 8)
print(ht1)
print(ht2)
dev.off()

# combine motifMatrix and geneScore matrix
corGIM_MM <- correlateTrajectories(trajGIM, trajMM,corCutOff = 0.3)
idxToRemove2 <- grep(pattern = "deviations", x = corGIM_MM[["correlatedMappings"]]$name2)
if(length(idxToRemove2 > 1)){
  corGIM_MM[["correlatedMappings"]] <- corGIM_MM[["correlatedMappings"]][-idxToRemove2,]
}
# corGIM_MM[[1]]

trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
# ComplexHeatmap::draw(ht1 + ht2)
pdf("EndospermTrajectory_umap_reanchor_plot8.pdf", width = 6, height = 8)
print(ht1)
print(ht2)


#====================================================================================================================#
#                           画关注的基因的开放分数，表达分数以及motif富集（针对转录因子）分数                                  #
#====================================================================================================================#
projAll <- addImputeWeights(projAll)
markerGenes  <- c("AT3G26744", "AT1G75080")
p1 <- plotEmbedding(
    ArchRProj = projAll, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAPTest",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projAll),
    pal = paletteContinuous(set = "solarExtra")
)
p2 <- plotEmbedding(
    ArchRProj = projAll, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    embedding = "UMAPTest",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projAll),
    pal = paletteContinuous(set = "solarExtra")
)
markerGenes  <- c("z:AT3G26744", "z:AT1G75080")
p3 <- plotEmbedding(
    ArchRProj = projAll, 
    colorBy = "MotifMatrix", 
    name = markerGenes, 
    embedding = "UMAPTest",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(projAll),
    pal = paletteContinuous(set = "solarExtra")
)
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure9_V3_UMAPtest.pdf",width = 15, height = 8)
grid.arrange(p1[[1]],p2[[1]],p3[[1]],p1[[2]],p2[[2]],p3[[2]],nrow=2)
dev.off()

#====================================================================================================================#
#                                     计算细胞类型相关性以及marker基因的表达特异性热图                                      #
#====================================================================================================================#
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
cell_ids <- seRNA@meta.data$cellID
cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
seRNA_filtered <- subset(seRNA, cells = cells_to_keep)

# 方式一，RNA高变基因以及跟ATAC overlap的基因，采用这种方法
marker <- VariableFeatures(seRNA_filtered)
# # 方式二，找RNA差异表达的基因，效果不好，不采用这种方法
# Idents(seRNA_filtered) <- seRNA_filtered$endotype2025
# seRNA_filtered.markers <- FindAllMarkers(seRNA_filtered, only.pos = TRUE)
# marker <- seRNA_filtered.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
# # 找差异开放的基因
# projAll$endotype2025_merged <- projAll$predictedGroup_endo2025
# projAll$endotype2025_merged[projAll$endotype2025_merged == "MCE-like"] <- "MCE"
# markersGS <- getMarkerFeatures(ArchRProj = projAll, useMatrix = "GeneIntegrationMatrix", groupBy = "endotype2025_merged", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
# markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# common_clusters <- intersect(unique(marker$cluster), names(markerList))
# mergedMarkers <- lapply(common_clusters, function(clust) {
#   rna_genes <- marker %>% dplyr::filter(cluster == clust) %>% dplyr::pull(gene)
#   atac_genes <- markerList[[clust]]$name
#   intersect(rna_genes, atac_genes)
# })
# names(mergedMarkers) <- common_clusters
# marker <- c(mergedMarkers[[1]],mergedMarkers[[2]],mergedMarkers[[3]],mergedMarkers[[4]],mergedMarkers[[5]])
# gene2 <- unique(c(markerList[[1]][,5],markerList[[2]][,5],markerList[[3]][,5],markerList[[4]][,5],markerList[[5]][,5]))

# 同时取ATAC基因的交集，使用1k个基因
# genes_of_interest <- intersect(marker[,7,drop=T], gim_genes)
genes_of_interest <- intersect(marker, gim_genes)
genes_of_interest <- intersect(marker_genes, gim_genes)
genes_of_interest <- head(genes_of_interest, 200)
# RNA 原始数据（可以用 scale.data 或 data）
rna_mat <- GetAssayData(seRNA_filtered, assay = "RNA", slot = "data")[genes_of_interest, ]
# 可以按细胞类型汇总平均表达
rna_avg <- aggregate(t(as.matrix(rna_mat)), by = list(celltype = seRNA_filtered$endotype2025), FUN = mean)
rownames(rna_avg) <- rna_avg$celltype
rna_avg <- t(rna_avg[,-1])

projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_filter_motifV2meme.rds")
# cells.keep <- rownames(projAll)[projAll$predictedScore_endo2025 > 0.6]
# projAll <- projAll[cells.keep, ] 
gim <- getMatrixFromProject(projAll, useMatrix = "GeneScoreMatrix")
gim_genes <- rowData(gim)$name
gene_scores <- assay(gim, "GeneScoreMatrix")
rownames(gene_scores) <- rowData(gim)$name

# 提取基因和细胞
atac_mat <- gene_scores[genes_of_interest, ]

# 按 ATAC 细胞类型平均
atac_avg <- aggregate(t(as.matrix(atac_mat)), by = list(celltype = projAll$predictedGroup_endo2025), FUN = mean)
rownames(atac_avg) <- atac_avg$celltype
atac_avg <- t(atac_avg[,-1])
cor_matrix1 <- cor(rna_avg, atac_avg, method = "pearson")  
# 行：RNA 细胞类型，列：ATAC 细胞类型

row_order <- c("SE", "MCE", "MCE-like", "ESR", "ESR_PEN", "PEN", "CZE","CZE_PEN")
col_order <- c( "ESR", "ESR_PEN", "PEN", "CZE")

cor_matrix <- cor_matrix1[row_order, col_order]
# cor_matrix_ordered2 <- cor_matrix[row_order, col_order]
# my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure7.pdf")
pheatmap(cor_matrix, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, main = "RNA vs ATAC cell type correlation")
dev.off()

symbol <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/geneID_to_symbol_rmNA.txt", header=T, sep="\t", stringsAsFactors = F, quote="",check.names = FALSE,comment.char = "")
# cosg
marker_cosg <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/EndoCells_pseudotime.DEG_cosg.tsv", header=T, stringsAsFactors=F)
top20_per_cluster <- marker_cosg %>% group_by(cluster) %>% arrange(score, .by_group = TRUE) %>% slice_head(n = 50) %>% ungroup()
# wilcox
marker_wilcox <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/EndoCells_pseudotime.DEG_wilcox.tsv", , header=T, stringsAsFactors=F)
top20_per_cluster <- marker_wilcox %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice_head(n = 50) %>% ungroup()

top20_genes <- Reduce(intersect, list(rownames(gene_scores), top20_per_cluster$gene, symbol$gene_id))
top20_per_cluster <- top20_per_cluster[which(top20_per_cluster$gene %in% top20_genes),]

atac_mat <- gene_scores[top20_genes, ]
atac_avg <- aggregate(t(as.matrix(atac_mat)), by = list(celltype = projAll$predictedGroup_endo2025), FUN = mean)
rownames(atac_avg) <- atac_avg$celltype
atac_avg <- t(atac_avg[,-1])

# 对齐注释
atac_avg_df <- atac_avg %>% as.data.frame() %>% tibble::rownames_to_column("gene")
anno <- top20_per_cluster %>% left_join(symbol, by = c("gene" = "gene_id"))
anno_aligned <- atac_avg_df %>% left_join(anno %>% select(gene, cluster, symbol), by = "gene")
# 按照cluster排序，最终的基因排序
anno_sorted <- anno_aligned %>% arrange(cluster)
# 创建行名（gene + cluster，保证唯一），这里整合了所有的元素
anno_sorted <- anno_sorted %>% mutate(gene_cluster = paste(gene, symbol,cluster, sep = "_"))

# 过滤掉没有的细胞类型以及表达模式不一致的基因
anno_sorted <- anno_sorted[which(anno_sorted$cluster %in% c("CZE","ESR","ESR_PEN","PEN")),]

anno_sorted <- anno_sorted %>% rowwise() %>% filter(get(cluster) == max(c_across(c(CZE, ESR, ESR_PEN, PEN)))) %>% ungroup()
write.table(anno_sorted, "RNA_cosg_ATAC.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(anno_sorted, "RNA_wilcox_ATAC.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# 提取表达矩阵
mat <- anno_sorted %>% column_to_rownames("gene_cluster") %>% select(-gene, -cluster, -symbol)

mat_filtered <- mat[apply(mat, 1, sd) != 0, ]
mat_z <- t(scale(t(as.matrix(mat_filtered))))
anno_sorted <- anno_sorted[, 6, drop = FALSE]
rownames(anno_sorted) <- rownames(mat)
anno_sorted <- anno_sorted[rownames(mat_z),,drop=F]

# 创建颜色向量
library(RColorBrewer)
cluster_levels <- unique(anno_sorted$cluster)
cluster_colors <- setNames(brewer.pal(length(cluster_levels), "Set2"), cluster_levels)
annotation_colors <- list(cluster = cluster_colors)

pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure8_cosg_score.pdf",width = 10,height = 15)
pheatmap(mat_z,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_row = anno_sorted,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         show_colnames = TRUE)
dev.off()

#====================================================================================================================#
#                                                   提取peak marker                                               ####
#====================================================================================================================#

markerPeaks <- getMarkerFeatures(
  ArchRProj = projAll, 
  useMatrix = "PeakMatrix", 
  groupBy = "Stage",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
for(stage in names(markerList)){
  write.table(
    as.data.frame(markerList[[stage]]),
    file = paste0(stage, "_specific_peaks.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

#====================================================================================================================#
#                                               motif足迹分析                                                     ####
#====================================================================================================================#

motifPositions <- getPositions(projAll)
motifs <- c("AT1G75080")
# 为了准确找到TF足迹，我们需要大量的reads。因此，细胞需要进行分组生成拟混池ATAC-seq谱才能用于TF足迹分析。这些拟混池谱之前在peak鉴定时就已经保存为分组覆盖文件。 如果没有在ArchRProject添加分组覆盖信息，则运行如下命令
# Stage
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage",force = TRUE)
seFoot <- getFootprints(ArchRProj = projAll, positions = motifPositions[motifs], groupBy = "Stage")
# cell type,这个画不出来
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "predictedGroup_endo2025",force = TRUE)
seFoot <- getFootprints(ArchRProj = projAll, positions = motifPositions[motifs], groupBy = "predictedGroup_endo2025")
# 当我们获取了这些足迹，我们可以使用plotFootprints()函数进行展示。该函数能够同时以多种方式对足迹进行标准化。
# Tn5偏好的足迹标准化
# 减去Tn5偏好，这里会把图画到/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/Merged_Endosperm_Filtered_V2/Plots/这个位置
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure10_subtract.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, plotName = "Footprints-Subtract-Bias", ArchRProj = projAll, normMethod = "Subtract", addDOC = FALSE, smoothWindow = 5)
dev.off()
# 这些图保存在ArchRProject的outputDirectory。如果需要绘制所有motif, 可以将其返回为ggplot2对象，需要注意这个ggplot对象会非常大
# 除以Tn5偏好
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure10_devide.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, ArchRProj = projAll, normMethod = "Divide", plotName = "Footprints-Divide-Bias", addDOC = FALSE, smoothWindow = 5)
dev.off()
# 无Tn5偏好标准化的足迹
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_figure10_normal.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, ArchRProj = projAll,normMethod = "None", plotName = "Footprints-No-Normalization", addDOC = FALSE, smoothWindow = 5)
dev.off()
# 与TSS的距离
seTSS <- getFootprints(ArchRProj = projAll, positions = GRangesList(TSS = getTSS(projAll)), groupBy = "Stage", flank = 2000)
plotFootprints(
  seFoot = seTSS,
  ArchRProj = projHeme5, 
  normMethod = "None",
  plotName = "TSS-No-Normalization",
  addDOC = FALSE,
  flank = 2000,
  flankNorm = 100
)

#====================================================================================================================#
#                                 combine motifMatrix and geneInteration matrix                                   ####
#====================================================================================================================#
gene <- c("AT1G71250","AT1G21970","AT1G08560","AT1G30690","AT1G48930","AT3G62830","AT5G13170","AT3G08900","AT1G68560","AT2G47650","AT2G25540","AT4G33600","AT3G19820","AT4G04460","AT1G62340","AT5G50750","AT3G48740","AT5G50260","AT2G19500","AT1G23320","AT2G35670","AT5G60440","AT2G35230","AT1G65330","AT1G67775","AT1G50650","AT1G62340","AT3G08900","AT2G47650","AT1G49770","AT1G68560","AT3G11520","AT1G48930","AT1G75080","AT5G07280","AT1G71250","AT3G26744","AT3G48740","AT1G11190","AT5G18270")
# plot motifMatrix
trajMM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, varCutOff = 0.9, pal = paletteContinuous(set = "solarExtra"),labelTop =150)
# plot geneScoreMatrix 
trajGSM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
# plot geneIntegrationMatrix 
trajGIM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
# trajPM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
# p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra")) 

# select gene to plot
pattern <- paste(gene, collapse = "|")
motifs <- grep(pattern, rownames(trajGSM), value = TRUE)
# motifs_z <- motifs[grep("^z:", motifs)]
trajGSM_sub <- trajGSM[motifs, ]
trajGIM_sub <- trajGIM[motifs, ]
mat <- assay(trajGSM_sub)
scaled_mat <- t(scale(t(mat)))
p2_2 <- pheatmap(scaled_mat, color = paletteContinuous(set = "horizonExtra"), cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)
mat <- assay(trajGIM_sub)
scaled_mat <- t(scale(t(mat)))
p3_2 <- pheatmap(scaled_mat, color = paletteContinuous(set = "horizonExtra"), cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)

pdf("anchor_ED/pdf/EndospermTrajectory_umap_trajectory_k50_09.pdf", width = 6, height = 8)
print(p1)
print(p2)
print(p2_2)
print(p3)
print(p3_2)
dev.off()

## combine motifMatrix and geneInteration matrix
#corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
#idxToRemove <- grep(pattern = "deviations", x = corGSM_MM[["correlatedMappings"]]$name2)
#if(length(idxToRemove > 1)){
#  corGSM_MM[["correlatedMappings"]] <- corGSM_MM[["correlatedMappings"]][-idxToRemove,]
#}
#
## corGSM_MM[["correlatedMappings"]]
#trajGSM2 <- trajGSM[corGSM_MM[["correlatedMappings"]]$name1, ]
#trajMM2 <- trajMM[corGSM_MM[["correlatedMappings"]]$name2, ]
#
#trajCombined <- trajGSM2
#assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
#
#combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
#rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
#
#ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
#ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, force = TRUE)
#
## ComplexHeatmap::draw(ht1 + ht2)
#pdf("EndospermTrajectory_umap_reanchor_plot7.pdf", width = 6, height = 8)
#print(ht1)
#print(ht2)
#dev.off()
#
## combine motifMatrix and geneScore matrix
#corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
#idxToRemove2 <- grep(pattern = "deviations", x = corGIM_MM[["correlatedMappings"]]$name2)
#if(length(idxToRemove2 > 1)){
#  corGIM_MM[["correlatedMappings"]] <- corGIM_MM[["correlatedMappings"]][-idxToRemove2,]
#}
## corGIM_MM[[1]]
#
#trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
#trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
#trajCombined <- trajGIM2
#assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
#combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
#
#rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
#ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
#ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
## ComplexHeatmap::draw(ht1 + ht2)
#pdf("EndospermTrajectory_umap_reanchor_plot8.pdf", width = 6, height = 8)
#print(ht1)
#print(ht2)
#dev.off()

#====================================================================================================================#
#                                        整合分析:peak关联gene调控                                                ####
#====================================================================================================================#
# 可以只用ATAC-seq数据进行分析，如识别peak之间的共开放性来预测调控相互作用，也可以整合scRNA-seq数据，如通过peak-基因的连锁分析预测增性子活性
# 查看peak之间的共开放相关性
projAll <- addCoAccessibility(ArchRProj = projAll, reducedDims = "IterativeLSI")
# cA <- getCoAccessibility(ArchRProj = projAll, corCutOff = 0.5, resolution = 1, returnLoops = FALSE)
# metadata(cA)[[1]]
# returnLoops的参数含义是什么？resolution的含义是什么？
# cA <- getCoAccessibility(ArchRProj = projAll, corCutOff = 0.5, resolution = 1, returnLoops = TRUE)
# cA <- getCoAccessibility(ArchRProj = projAll, corCutOff = 0.5, resolution = 1000, returnLoops = TRUE)
# cA <- getCoAccessibility(ArchRProj = projAll, corCutOff = 0.5, resolution = 10000, returnLoops = TRUE)
# 在browser track中绘制基因区域在不同细胞类型中的共开放
markerGenes  <- c("AT1G75080")
p1 <- plotBrowserTrack(ArchRProj = projAll, groupBy = "Stage", geneSymbol = markerGenes, upstream = 50000, downstream = 50000, loops = getCoAccessibility(projAll))
p2 <- plotBrowserTrack(ArchRProj = projAll, groupBy = "predictedGroup_endo", geneSymbol = markerGenes, upstream = 50000, downstream = 50000, loops = getCoAccessibility(projAll))
pdf("Endosperm-Tracks-Marker-Genes-with-CoAccessibility_plot1.pdf",width = 25,height = 10)
grid.arrange(p1[[1]],p2[[1]],nrow=1)
dev.off()

# Peak2GeneLinkages peak关联基因分析
projAll <- addPeak2GeneLinks(ArchRProj = projAll, reducedDims = "IterativeLSI")
# p2g <- getPeak2GeneLinks(ArchRProj = projAll, corCutOff = 0.45, resolution = 1, returnLoops = FALSE)
# metadata(p2g)[[1]]
# p2g <- getPeak2GeneLinks(ArchRProj = projAll, corCutOff = 0.45, resolution = 10000, returnLoops = TRUE)
# 在browser track中绘制每种细胞类型的peak-to-gene连接
markerGenes  <- c("AT1G75080")
p1 <- plotBrowserTrack(ArchRProj = projAll, groupBy = "Stage", geneSymbol = markerGenes, upstream = 50000, downstream = 50000, loops = getPeak2GeneLinks(projAll))
p2 <- plotBrowserTrack(ArchRProj = projAll, groupBy = "predictedGroup_endo", geneSymbol = markerGenes, upstream = 50000, downstream = 50000, loops = getPeak2GeneLinks(projAll))

pdf("Endosperm-Tracks-Marker-Genes-with-Peak2GeneLinks_plot2.pdf",width = 25,height = 10)
grid.arrange(p1[[1]],p2[[1]],nrow=1)
dev.off()

# 绘制Peak-to-gene连接热图
p3 <- plotPeak2GeneHeatmap(ArchRProj = projAll, groupBy = "Stage")
p4 <- plotPeak2GeneHeatmap(ArchRProj = projAll, groupBy = "predictedGroup_endo")

pdf("Endosperm-Tracks-Marker-Genes-with-Peak2Gene-heatmap.pdf",width = 25,height = 10)
#grid.arrange(p3[[1]],p4[[1]],nrow=1)
draw(p3)
draw(p4)
dev.off()

# 保存
saveRDS(projAll, "")

