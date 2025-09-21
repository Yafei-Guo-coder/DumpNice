library(matrixStats, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib")
packageVersion("matrixStats")
library(ArchR)
library(BSgenome.Athaliana.TAIR10.02022009)
library(parallel)
library(ggalt)
library(clustree)
library(dplyr)
library(patchwork)
library(gridExtra)
library(grid)
.libPaths("/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib")
# devtools::install_version("matrixStats", version="1.5.0", lib = "./lib")
# install.packages("universalmotif", lib = "./lib")
library(universalmotif)
library(tidydr)
library(ggalluvial)
library(ggplot2)
library(chromVAR)
library(TFBSTools)
library(SummarizedExperiment)
library(motifmatchr)
library(Matrix)
library(rhdf5)
library(readxl)
library(plyr)
library(tidyverse)
library(cluster)
# library(hrbrthemes)
library(viridis)
library(pheatmap)
# library(SeuratData)
# install the dataset and load requirements
# InstallData("pbmcMultiome")
library(SeuratObject, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib")
library(Seurat)
library(Signac)
#library(EnsDb.Hsapiens.v86)
library(cowplot)
library(parallel)
library(ComplexHeatmap)
library(S4Vectors)

# 从BSgenome对象中创建基因组注释
genome_annotation <- createGenomeAnnotation(genome = BSgenome.Athaliana.TAIR10.02022009)

#====================================================================================================================#
###                                             create Arrow Files                                                 ###
#====================================================================================================================#

blacklist_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/BlackList_TAIR10.bed"
bl <- read.table(blacklist_file)
bl <- bl[,1:3]
colnames(bl) <- c("chr", "start", "end")
bl$strand <- "+" #问题1:为什么要把链都设置成+
# 将数据框转换为GRanges对象，保留额外的列
blacklist <- makeGRangesFromDataFrame(bl, keep.extra.columns = T)
# 将黑名单添加到基因组注释对象中
genome_annotation$blacklist <- blacklist

# 定义工具脚本路径
gene_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.gene.txt"
trans_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.transcript.txt"
exon_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.exon.txt"

#gene_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/00_ref/Arabidopsis/Araport11/Araport11.Mar92021.gene.noCM.txt"
#trans_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/00_ref/Arabidopsis/Araport11/Araport11.Mar92021.transcript.noCM.txt"
#exon_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/00_ref/Arabidopsis/Araport11/Araport11.Mar92021.exon.noCM.txt"
archr_utils <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/archr_utils.R"
# 加载工具脚本
source(archr_utils)

# 加载基因信息，包括基因、外显子和转录本
gene_annotation <- loadGeneinfo(gene = gene_file, exon = exon_file, trans = trans_file)
# 修改外显子信息中的 "gene_id" 列，去掉 ":exon:" 后的部分
gene_annotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]] <- unlist(lapply(gene_annotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]], function(x){unlist(strsplit(x, ":exon:"))[1]}))
# 修改外显子信息中的 "symbol" 列，去掉 ":exon:" 后的部分
gene_annotation@listData[["exons"]]@elementMetadata@listData[["symbol"]] <- unlist(lapply(gene_annotation@listData[["exons"]]@elementMetadata@listData[["symbol"]], function(x){unlist(strsplit(x, ":exon:"))[1]}))
# 获取指定目录下所有以 ".R" 结尾的R脚本文件名
rscripts <- list.files(path = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/ArchR_R", pattern = ".R$", recursive = F, full.names = F)
# 循环加载每个R脚本文件
for (i in rscripts) {source(paste0("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/ArchR_R/",i))}

# 所有时期的 RDS 文件
stage_files <- list.files(path = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds", pattern = "ArchR_.*_anchor\\.rds$", full.names = TRUE)
proj_list <- list()
for (f in stage_files) {
  proj <- readRDS(f)
  stage <- sub("ArchR_(.*)_anchor\\.rds", "\\1", basename(f))
  if (!"predictedGroup" %in% names(proj@cellColData)) {
    stop(paste("No predictedGroup in:", f))
  }
  endo_cells <- which(grepl("endosperm", proj$predictedGroup) & proj$predictedScore > 0.5)
  proj_endo <- subsetArchRProject(
    ArchRProj = proj, 
    cells = proj$cellNames[endo_cells], 
    outputDirectory = paste0("Subset_", stage),
    dropCells = TRUE,
    force = TRUE
  )
  proj_endo$Stage <- stage
  proj_list[[stage]] <- proj_endo
}

#====================================================================================================================#
###                                        create archR project                                                   ####
#====================================================================================================================#

# create ArchR project from arrow file
# Get all files and folders in the current directory
all_files <- list.files(path = ".", full.names = TRUE)
# Filter out directories starting with Subset_
subset_dirs <- all_files[file.info(all_files)$isdir & grepl("^./Subset_", all_files)]
# View Results
# print(subset_dirs)
# Look for .arrow files in each Subset_*/ArrowFiles/
arrow_files_all <- unlist(lapply(subset_dirs, function(dir) {
   arrow_dir <- file.path(dir, "ArrowFiles")
   list.files(arrow_dir, pattern = "\\.arrow$", full.names = TRUE)
 }))
projAll <- ArchRProject(
   ArrowFiles = arrow_files_all,       
   geneAnnotation = gene_annotation,
   genomeAnnotation = genome_annotation,
   outputDirectory = "Merged_Endosperm",   
   showLogo = F,           
   copyArrows = F,  
   threads = 8
 )
# Assume Sample is something like "A250415004_bent1"
projAll$Stage <- sub(".*_(globular|heart|torpedo|bent|cotyledon)[0-9]*", "\\1", projAll$Sample)
# table(projAll$Stage)

# 1. 收集所有 endosperm 细胞名
endo_cells_all <- unlist(lapply(proj_list, function(p) p$cellNames))
endo_cells_all_fixed <- sapply(endo_cells_all, function(cell) {
  # 提取组织名，例如从 cell 名里或其他表里获取
  # 下面是示例，请根据你的命名规则调整
  if (grepl("bent", cell)) {
    paste0("bent.", cell)
  } else if (grepl("globular", cell)) {
    paste0("globular.", cell)
  } else if (grepl("torpedo", cell)) {
    paste0("torpedo.", cell)
  } else if (grepl("heart", cell)) {
    paste0("heart.", cell)
  } else if (grepl("cotyledon", cell)) {
    paste0("cotyledon.", cell)
  } else {
    cell  # fallback
  }
})

projAll <- addGeneScoreMatrix(input = projAll, force = TRUE)

# 2. 构建 projAll 后，筛选保留这些细胞
projAll <- subsetArchRProject(
  ArchRProj = projAll,
  cells = endo_cells_all,
  outputDirectory = "Merged_Endosperm_Filtered",
  dropCells = TRUE,
  force = TRUE)

all_colnames <- unique(unlist(lapply(proj_list, function(p) colnames(p@cellColData))))

# 合并所有 cell metadata，自动补 NA
all_metadata <- do.call(rbind, lapply(proj_list, function(p) {
  meta <- as.data.frame(p@cellColData)
  meta$cellNames <- rownames(meta)
  
  # 补齐缺失的列
  missing_cols <- setdiff(all_colnames, colnames(meta))
  for (col in missing_cols) {
    meta[[col]] <- NA
  }
  
  # 保证列顺序一致（含 cellNames）
  meta <- meta[, c(all_colnames, "cellNames")]
  return(meta)
}))

# 覆盖新项目的 metadata
projAll@cellColData <- S4Vectors::DataFrame(all_metadata)
rownames(projAll@cellColData) <- projAll@cellColData$cellNames

# Consider removing the A250410002_torpedo2 library
cells_to_keep <- projAll$cellNames[projAll$Sample != "A250410002_torpedo2"]
projAll <- projAll[cells_to_keep, ]
projAll <- addIterativeLSI(projAll, useMatrix = "MotifMatrix", name = "IterativeLSI")

projAll <- addUMAP(projAll, reducedDims = "IterativeLSI", name = "UMAP")
projAll <- addTSNE(ArchRProj = projAll, reducedDims = "IterativeLSI", name = "TSNE", perplexity = 30)

saveRDS(projAll, "ArchR_EndospermD_UMAP_rmTorpedo2.rds")

p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Clusters_merge", embedding = "UMAP")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
pdf("ArchR_EndospermD_UMAP_2.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "UMAP - Endosperm Development")
dev.off()

p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "TSNE")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "TSNE")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Clusters_merge", embedding = "TSNE")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup", embedding = "TSNE")
pdf("ArchR_EndospermD_TSNE_2.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "UMAP - Endosperm Development")
dev.off()

#====================================================================================================================#
###                                                motif                                                          ####
#====================================================================================================================#

# 单独根据每一组细胞（例如聚类）鉴定peak，这样可以分析出不同组的特异性peak
# 拟混池重复
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "predictedGroup", threads = 1)
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "predictedGroup", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addPeakMatrix(projAll)
# getAvailableMatrices(projAll)

# 保存
saveRDS(projAll, "ArchR_EndospermD_harmony_rmTorpedo2_peak.rds")
projHeme5 <- proj

# 添加motif
# projHeme5 <- readRDS("ArchR_EndospermD_harmony_rmTorpedo2_peak.rds")
meme_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/Ath_TF_binding_motifs.meme" 
motifs <- read_meme(meme_file)
# 转换为PWMatrix格式
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
projAll <- addMotifAnnotations(ArchRProj = projAll, motifPWMs = pwmlist, name = "Motif")

# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)

# 保存
saveRDS(projAll, "ArchR_EndospermD_harmony_rmTorpedo2_motif.rds")

# getAvailableMatrices(projAll)
projAll <- readRDS("ArchR_EndospermD_harmony_rmTorpedo2_motif.rds")

plotVarDev <- getVarDeviations(projAll, name = "MotifMatrix", plot = TRUE)
dev.off()
pdf("motif1-root-0316.pdf",width = 5,height = 5)
grid.arrange(plotVarDev,nrow=1)
dev.off()

# 画motif富集热图
motif_matrix <- getMatrixFromProject(ArchRProj = projHeme5, useMatrix = "MotifMatrix")
deviations <- assay(motif_matrix)
predictedGroup <- projHeme5$predictedGroup
deviation_by_predictedGroup <- sapply(unique(predictedGroup), function(ct) {
  cells <- which(predictedGroup == ct)
  rowMeans(deviations[, cells])
})
mat_z <- t(scale(t(deviation_by_predictedGroup)))
pdf("motif_heatmap1.pdf",width = 30,height = 20)
Heatmap(mat_z)
dev.off()

# 把motif聚类成家族进行富集分析
match_TF <- read.table("../motif_AT/TF_family.txt", header=T, stringsAsFactors = F)
match_TF <- match_TF[match_TF$TF %in% rownames(deviations), ]
deviation_mat_filtered <- deviations[match_TF$TF, ]
rownames(match_TF) <- match_TF$TF  
TF_family <- match_TF[rownames(deviation_mat_filtered), "Family"]
predictedGroups <- projHeme5$predictedGroup
predictedGroup_levels <- unique(predictedGroups)
dev_by_predictedGroup <- sapply(predictedGroup_levels, function(ct) {
  cells <- projHeme5$cellNames[projHeme5$predictedGroup == ct]
  rowMeans(deviation_mat_filtered[, cells, drop=FALSE])
})
rownames(dev_by_predictedGroup) <- rownames(deviation_mat_filtered)
dev_family <- as.data.frame(dev_by_predictedGroup)
dev_family$Family <- TF_family
dev_family_avg <- dev_family %>%
  group_by(Family) %>%
  summarise(across(where(is.numeric), mean))
dev_mat <- as.matrix(dev_family_avg[,-1])
rownames(dev_mat) <- dev_family_avg$Family
dev_mat_z <- t(scale(t(dev_mat)))

pdf("motif_heatmap2.pdf",width = 30,height = 20)
Heatmap(dev_mat_z, name = "TF Family\nDeviation Z-score", 
        cluster_rows = TRUE, cluster_columns = TRUE)
dev.off()

#====================================================================================================================#
####                                              trajectory                                                       ###
#====================================================================================================================#

projAll <- readRDS("ArchR_EndospermD_harmony_rmTorpedo2_motif.rds")

# anchor ##
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig4.endo_traj/Fig4a.endo_recluster/EndoCells_recluster.rds")
projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = FALSE, force = TRUE, groupRNA = "endotype0514", nameCell = "predictedCell_endo", nameGroup = "predictedGroup_endo", nameScore = "predictedScore_endo", threads = 1)

saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor.rds")

projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor.rds")
# 过滤
idxPass <- which(projAll$predictedScore_endo > 0.6)
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]

# 查看各 cluster 的时间分布，辅助选择轨迹顺序
table(projAll$Stage, projAll$predictedGroup)

# 假设经验上知道 trajectory 顺序是以下几个 cluster
trajectory_order <- c("globular", "heart", "torpedo", "bent", "cotyledon")

projAll <- addTrajectory(
  ArchRProj = projAll,
  name = "EndospermTrajectory",
  groupBy = "Stage",
  trajectory = trajectory_order,
  postFilterQuantile = 0.9,
  #reducedDims = "IterativeLSI",
  embedding = "UMAP", 
  force = TRUE
)

pdf("EndospermTrajectory_umap_anchor-zhijie-0.5.pdf", width = 6, height = 6)
p <- plotTrajectory(projAll, trajectory = "EndospermTrajectory", colorBy = "cellColData", name = "EndospermTrajectory")
print(p)
dev.off()

# 看预测结果
#projAll <- addUMAP(projAll, reducedDims = "IterativeLSI", name = "UMAP")

p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedScore_endo", embedding = "UMAP")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup_endo", embedding = "UMAP")

pdf("ArchR_EndospermD_UMAP_anchorEndo_zhijie-0.5.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "UMAP anchor - Endosperm Development")
dev.off()
dev.off()

# 结合伪时间和真实细胞
library(dplyr)
library(ggplot2)
cell_df <- as.data.frame(getCellColData(projAll))
# 假设你的数据框是 cell_df，包含 pseudotime 和 Stage 两列
# 1. 将 pseudotime 分成 50 个bin
cell_df <- cell_df %>%
  mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))

# 2. 统计每个 bin 内不同 Stage 细胞数量
bin_stage_counts <- cell_df %>%
  group_by(bin, Stage) %>%
  summarise(count = n(), .groups = 'drop')

# 3. 计算每个 bin 细胞总数
bin_totals <- bin_stage_counts %>%
  group_by(bin) %>%
  summarise(total = sum(count), .groups = 'drop')

# 4. 合并并计算比例
bin_stage_props <- bin_stage_counts %>%
  left_join(bin_totals, by = "bin") %>%
  mutate(prop = count / total)

# 5. 画堆积柱状图显示比例
pdf("Endosperm_trajectory_box1_umap-peak-zhijie-0.5.pdf")
ggplot(bin_stage_props, aes(x = bin, y = prop, fill = Stage)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "Stage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

cell_df$predictedGroup_add <- paste(cell_df$predictedGroup_endo,cell_df$Stage,sep="_")
# 假设伪时间列是 pseudotime，先做bin
cell_df <- cell_df %>%
  mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))

# 统计每个bin中不同predictedGroup_add的数量
bin_group_counts <- cell_df %>%
  group_by(bin, predictedGroup_add) %>%
  summarise(count = n(), .groups = 'drop')

# 计算每个bin总数
bin_totals <- bin_group_counts %>%
  group_by(bin) %>%
  summarise(total = sum(count), .groups = 'drop')

# 计算比例
bin_group_props <- bin_group_counts %>%
  left_join(bin_totals, by = "bin") %>%
  mutate(prop = count / total)

# 画堆积百分比柱状图

pdf("Endosperm_trajectory_box2_predictedGroup_umap-peak-zhijie-0.5.pdf")
ggplot(bin_group_props, aes(x = bin, y = prop, fill = predictedGroup_add)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "PredictedGroup_Stage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# 假设伪时间列是 pseudotime，先做bin
cell_df <- cell_df %>%
  mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))

# 统计每个bin中不同predictedGroup_add的数量
bin_group_counts <- cell_df %>%
  group_by(bin, predictedGroup_endo) %>%
  summarise(count = n(), .groups = 'drop')

# 计算每个bin总数
bin_totals <- bin_group_counts %>%
  group_by(bin) %>%
  summarise(total = sum(count), .groups = 'drop')

# 计算比例
bin_group_props <- bin_group_counts %>%
  left_join(bin_totals, by = "bin") %>%
  mutate(prop = count / total)

# 画堆积百分比柱状图

pdf("Endosperm_trajectory_box2_predictedGroup_endo_umap-peak-zhijie-0.5.pdf")
ggplot(bin_group_props, aes(x = bin, y = prop, fill = predictedGroup_endo)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "PredictedGroup_Stage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#====================================================================================================================#
####                         通过计算前2000个差异的peak然后构建空间分布然后再计算伪时序分布                                 ###
#====================================================================================================================#

projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor.rds")
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", threads = 1)
projAll <- addGeneScoreMatrix(input = projAll, force = TRUE)
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor_heatmap.rds")

# 单独根据每一组细胞（例如聚类）鉴定peak，这样可以分析出不同组的特异性peak
# 拟混池重复
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "Stage", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addPeakMatrix(projAll)
# getAvailableMatrices(projAll)
meme_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/Ath_TF_binding_motifs.meme" 
motifs <- read_meme(meme_file)
# 转换为PWMatrix格式
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
projAll <- addMotifAnnotations(ArchRProj = projAll, motifPWMs = pwmlist, name = "Motif")

# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)

saveRDS(projAll, "ArchR_EndospermD_rmTorpedo2_reanchor_motif.rds")

arrowFiles <- list.files(path = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/", pattern = "\\.arrow$", full.names = TRUE, recursive = FALSE)
projAll <- ArchRProject(
  ArrowFiles = arrowFiles,          # 输入的Arrow文件列表
  geneAnnotation = gene_annotation, # 基因注释对象
  genomeAnnotation = genome_annotation, # 基因组注释对象
  outputDirectory = "Development",       # 输出目录，用于存储项目的输出文件
  showLogo = F,                     # 是否显示ArchR的logo，设置为FALSE表示不显示
  copyArrows = F,                    # 是否复制Arrow文件，设置为FALSE表示不复制 # 推荐用T，如果你修改了Arrow文件，可以保留一个原始副本以备后用
  threads = 1
)

cells_subset <- projSub$cellNames                  
idxPass <- which(projAll$cellNames %in% cells_subset)
cellsPass <- projAll$cellNames[idxPass]
projAll2 <- projAll[cellsPass, ]
projAll2 <- addIterativeLSI(projAll2, useMatrix = "GeneScoreMatrix", name = "IterativeLSI")
projAll2 <- addUMAP(projAll2, reducedDims = "IterativeLSI", name = "UMAP")

# projAll2$EndospermTrajectory <- projSub$EndospermTrajectory
# projAll2$predictedCell_endo <- projSub$predictedCell_endo
# projAll2$predictedGroup_endo <- projSub$predictedGroup_endo
# projAll2$predictedScore_endo <- projSub$predictedScore_endo

stage_vec <- setNames(projSub$Stage, projSub$cellNames)
projAll2$Stage <- stage_vec[projAll2$cellNames]

# anchor
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig4.endo_traj/Fig4a.endo_recluster/EndoCells_recluster.rds")
projAll2 <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix2", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = FALSE, force = TRUE, groupRNA = "endotype0514", nameCell = "predictedCell_endo2", nameGroup = "predictedGroup_endo2", nameScore = "predictedScore_endo2", threads = 1)
projAll

p1 <- plotEmbedding(projAll2, colorBy = "cellColData", name = "Stage", embedding = "UMAP")
p2 <- plotEmbedding(projAll2, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p3 <- plotEmbedding(projAll2, colorBy = "cellColData", name = "predictedScore_endo", embedding = "UMAP")
p4 <- plotEmbedding(projAll2, colorBy = "cellColData", name = "predictedGroup_endo", embedding = "UMAP")

pdf("ArchR_EndospermD_UMAP_new.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "UMAP - Endosperm Development")
#grid.arrange(p1, p2, ncol = 2, top = "UMAP - Endosperm Development")
dev.off()

idxPass <- which(projAll$predictedScore_endo > 0.6)
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]

# 假设经验上知道 trajectory 顺序是以下几个cluster
trajectory_order <- c("globular", "heart", "torpedo", "bent", "cotyledon")

projSub <- addTrajectory(
  ArchRProj = projAll,
  name = "EndospermTrajectory",
  groupBy = "Stage",
  trajectory = trajectory_order,
  postFilterQuantile = 0.9,
  #reducedDims = "IterativeLSI",
  embedding = "UMAP", 
  force = TRUE
)
# 查看各 cluster 的时间分布，辅助选择轨迹顺序
table(projSub$Stage, projSub$predictedGroup)

pdf("EndospermTrajectory_umap_anchor-projAll3.pdf", width = 6, height = 6)
p <- plotTrajectory(projSub, trajectory = "EndospermTrajectory", colorBy = "cellColData", name = "EndospermTrajectory")
print(p)
dev.off()

common_cells <- intersect(cells_subset, getCellNames(projFiltered))
projFiltered@cellColData[common_cells, "Trajectory"] <- projSubset@cellColData[common_cells, "Trajectory"]
projFiltered@cellColData[common_cells, "snRNA_CellType"] <- projSubset@cellColData[common_cells, "snRNA_CellType"]

markers <- getMarkerFeatures(ArchRProj = projAll,useMatrix = "PeakMatrix",groupBy = "Stage",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"),useGroups = "cotyledon",bgdGroups = "bent")
topPeaks <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 1")
topPeakNames1 <- paste0(topPeaks[[1]]$seqnames, ":",topPeaks[[1]]$start, "-",topPeaks[[1]]$end)[1:500]

markers <- getMarkerFeatures(ArchRProj = projAll, useMatrix = "PeakMatrix", groupBy = "Stage", testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), useGroups = "bent", bgdGroups = "torpedo")
topPeaks <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 1")
topPeakNames2 <- paste0(topPeaks[[1]]$seqnames, ":", topPeaks[[1]]$start, "-", topPeaks[[1]]$end)[1:500]

markers <- getMarkerFeatures(ArchRProj = projAll, useMatrix = "PeakMatrix", groupBy = "Stage", testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), useGroups = "torpedo", bgdGroups = "heart")
topPeaks <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 1")
topPeakNames3 <- paste0(topPeaks[[1]]$seqnames, ":", topPeaks[[1]]$start, "-", topPeaks[[1]]$end)[1:500]

markers <- getMarkerFeatures(ArchRProj = projAll, useMatrix = "PeakMatrix", groupBy = "Stage", testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), useGroups = "heart", bgdGroups = "globular")
topPeaks <- getMarkers(markers, cutOff = "FDR <= 0.05 & Log2FC >= 1")
topPeakNames4 <- paste0(topPeaks[[1]]$seqnames, ":", topPeaks[[1]]$start, "-", topPeaks[[1]]$end)[1:500]

topPeakNames <- union(topPeakNames1,topPeakNames2,topPeakNames3,topPeakNames4)

peakMat <- getMatrixFromProject(projAll, useMatrix = "PeakMatrix")
peak_names <- paste0(seqnames(rowRanges(peakMat)), ":", 
                     start(rowRanges(peakMat)), "-", 
                     end(rowRanges(peakMat)))

mat <- assay(peakMat)
rownames(mat) <- peak_names
topPeakNames <- intersect(topPeakNames, rownames(mat))

# 子集化矩阵
peakMat_sub <- mat[topPeakNames, ]

reducedDims <- prcomp(t(peakMat_sub), scale. = TRUE)$x[, 1:10]

#projAll <- addReducedDims(projAll, reducedDims = reducedDims, name = "TopPeakPCA")
projAll@reducedDims[[2]][[1]] <- reducedDims
head(projAll@reducedDims[[2]][[1]])

trajectory_order <- c("globular", "heart", "torpedo", "bent", "cotyledon")
projAll <- addTrajectory(
  ArchRProj = projAll,
  name = "EndospermTrajectory",
  groupBy = "Stage",
  trajectory = trajectory_order,
  force = TRUE,
  reducedDims = "Harmony"  # 你新加的降维结果名
)
pdf("EndospermTrajectory_umap-peak.pdf", width = 6, height = 6)
p <- plotTrajectory(projAll,embedding = "UMAP-peak", trajectory = "EndospermTrajectory", colorBy = "cellColData", name = "EndospermTrajectory")
print(p)
dev.off()
dev.off()

#====================================================================================================================#
###                                        trajectory heatmap                                                      ###
#====================================================================================================================#

trajMM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, varCutOff = 0.8, pal = paletteContinuous(set = "solarExtra"),labelTop =200)
trajGSM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
TF <- c("AT3G26744","AT1G75080","AT3G19820","AT1G48930","AT3G48740","AT1G68560","AT1G21970","AT1G62340","AT2G47650","AT3G62830","AT3G08900","AT1G08560","AT1G30690","AT5G50750","AT1G71250","AT4G04460","AT2G25540","AT5G13170","AT4G33600","AT5G50260")

gene <- c("AT1G71250","AT1G21970","AT1G08560","AT1G30690","AT1G48930","AT3G62830","AT5G13170","AT3G08900","AT1G68560","AT2G47650","AT2G25540","AT4G33600","AT3G19820","AT4G04460","AT1G62340","AT5G50750","AT3G48740","AT5G50260","AT2G19500","AT1G23320","AT2G35670","AT5G60440","AT2G35230","AT1G65330","AT1G67775","AT1G50650","AT1G62340","AT3G08900","AT2G47650","AT1G49770","AT1G68560","AT3G11520","AT1G48930","AT1G75080","AT5G07280","AT1G71250","AT3G26744","AT3G48740","AT1G11190","AT5G18270")
# select motif to plot
pattern <- paste(TF, collapse = "|")
# 查找 trajMM 中行名包含这些 TF ID 的 motif
motifs <- grep(pattern, rownames(trajMM), value = TRUE)
motifs_z <- motifs[grep("^z:", motifs)]
# 子集
trajGSM_sub <- trajGSM[motifs, ]
colData(trajGSM_sub) <- colData(trajGSM)
library(pheatmap)

mat <- assay(trajGSM_sub)  # 已经是 smoothMat
p2_2_2 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)    
mat <- assay(trajGIM_sub)       
p2_2_3 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)       
p2_2 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)            
p2_2 <- plotTrajectoryHeatmap(trajGSM_sub,  pal = paletteContinuous(set = "horizonExtra"))

gene <- c("AT1G71250","AT1G21970","AT1G08560","AT1G30690","AT1G48930","AT3G62830","AT5G13170","AT3G08900","AT1G68560","AT2G47650","AT2G25540","AT4G33600","AT3G19820","AT4G04460","AT1G62340","AT5G50750","AT3G48740","AT5G50260","AT2G19500","AT1G23320","AT2G35670","AT5G60440","AT2G35230","AT1G65330","AT1G67775","AT1G50650","AT1G62340","AT3G08900","AT2G47650","AT1G49770","AT1G68560","AT3G11520","AT1G48930","AT1G75080","AT5G07280","AT1G71250","AT3G26744","AT3G48740","AT1G11190","AT5G18270")
# select motif to plot
pattern <- paste(gene, collapse = "|")
# 查找 trajMM 中行名包含这些 TF ID 的 motif
motifs <- grep(pattern, rownames(trajGSM), value = TRUE)
#motifs_z <- motifs[grep("^z:", motifs)]
# 子集
trajGSM_sub <- trajGSM[motifs, ]
trajGIM_sub <- trajGIM[motifs, ]

colData(trajGSM_sub) <- colData(trajGSM)
library(pheatmap)
mat <- assay(trajGSM_sub)  # 已经是 smoothMat
rownames(mat) <- rownames(trajGSM_sub)  # 保证有行名
p2_2 <- pheatmap(mat, color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)     
p2_3 <- plotTrajectoryHeatmap(trajGIM_sub,  pal = paletteContinuous(set = "horizonExtra"))
p2_2 <- plotTrajectoryHeatmap(trajGSM_sub,  pal = paletteContinuous(set = "horizonExtra"))
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, p4, name = "Plot-EndospermTrajectory-Traj-Heatmaps.pdf", ArchRProj = projAll, addDOC = FALSE, width = 6, height = 8)
pdf("Endosperm_trajectory_heatmap_out2.pdf", width = 6, height = 8)
print(p1)
print(p2)
print(p3)
print(p4)
print(p2_2_2)
print(p2_2_3)
dev.off()

# 提取该基因/TF的轨迹数据（这里假设是log2Norm后的矩阵smoothMat）
mat <- assay(trajMM_sub, "smoothMat")  # 1行，n列（时间点）
# 转成长格式方便绘图
df <- data.frame(
  Time = colnames(mat),
  Value = as.numeric(mat[1, ])
)
pdf("Endosperm_trajectory_heatmap_TF.pdf")
ggplot(df, aes(x = Time, y = Value, group = 1)) +
  geom_line() +
  geom_point() +
  ggtitle(rownames(trajMM_sub)) +
  theme_minimal() +
  xlab("Trajectory") +
  ylab("Motif activity (smoothed)")
dev.off()  
  
  
pdf("Endosperm_trajectory_heatmap_TF.pdf")
plotTrajectoryHeatmap(trajMM_sub, pal = paletteContinuous(set = "solarExtra"))
dev.off()
# select some motif to plot 
motifs_of_interest <- c("AT3G26744","AT1G75080","AT3G19820","AT1G48930","AT3G48740","AT1G68560","AT1G21970","AT1G62340","AT2G47650","AT3G62830","AT3G08900","AT1G08560","AT1G30690","AT5G50750","AT1G71250","AT4G04460","AT2G25540","AT5G13170","AT4G33600","AT5G50260") 
time_labels <- projAll$Stage
pseudotime_bins <- colnames(trajMM)
plotTrajectory(projAll, trajectory = "EndospermTrajectory", colorBy = "cellColData", name = "Stage")
# Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
idxToRemove <- grep(pattern = "deviations", x = corGSM_MM[["correlatedMappings"]]$name2)
if(length(idxToRemove > 1)){
  corGSM_MM[["correlatedMappings"]] <- corGSM_MM[["correlatedMappings"]][-idxToRemove,]
}
corGSM_MM[["correlatedMappings"]]
trajGSM2 <- trajGSM[corGSM_MM[["correlatedMappings"]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[["correlatedMappings"]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, force = TRUE)

ComplexHeatmap::draw(ht1 + ht2)

corGIM_MM <- correlateTrajectories(trajGIM, trajMM)

idxToRemove2 <- grep(pattern = "deviations", x = corGIM_MM[["correlatedMappings"]]$name2)
if(length(idxToRemove2 > 1)){
  corGIM_MM[["correlatedMappings"]] <- corGIM_MM[["correlatedMappings"]][-idxToRemove2,]
}
corGIM_MM[[1]]


trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)

ComplexHeatmap::draw(ht1 + ht2)

#====================================================================================================================#
###                               from 13000 cells to trajectory heatmap                                           ###
#====================================================================================================================#

 
projAll <- readRDS("ArchR_EndospermD_harmony_rmTorpedo2_motif.rds")
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", threads = 1)
projAll <- addGeneScoreMatrix(input = projAll, force = TRUE)
projAll <- addImputeWeights(projAll)

# anchor ##
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig4.endo_traj/Fig4a.endo_recluster/EndoCells_recluster.rds")
projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = TRUE, force = TRUE, groupRNA = "endotype0514", nameCell = "predictedCell_endo", nameGroup = "predictedGroup_endo", nameScore = "predictedScore_endo", threads = 1)
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor2.rds")
# projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor.rds")

# plot 看预测结果
p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedScore_endo", embedding = "UMAP")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup_endo", embedding = "UMAP")

pdf("ArchR_EndospermD_UMAP_reanchor_plot1.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "UMAP - Endosperm Development")
#grid.arrange(p1, p2, ncol = 2, top = "UMAP - Endosperm Development")
dev.off()

# 过滤
idxPass <- which(projAll$predictedScore_endo > 0.6)
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]

# 查看各 cluster 的时间分布，辅助选择轨迹顺序
table(projAll$Stage, projAll$predictedGroup)

# 假设经验上知道 trajectory 顺序是以下几个 cluster
trajectory_order <- c("globular", "heart", "torpedo", "bent", "cotyledon")

projAll <- addTrajectory(
  ArchRProj = projAll,
  name = "EndospermTrajectory",
  groupBy = "Stage",
  trajectory = trajectory_order,
  postFilterQuantile = 0.9,
  #reducedDims = "IterativeLSI",
  embedding = "UMAP", 
  force = TRUE
)

pdf("EndospermTrajectory_umap_reanchor_plot2.pdf", width = 6, height = 6)
p <- plotTrajectory(projAll, trajectory = "EndospermTrajectory", colorBy = "cellColData", name = "EndospermTrajectory")
print(p)
dev.off()

# 结合伪时间和真实细胞

cell_df <- as.data.frame(getCellColData(projAll))

# 假设你的数据框是 cell_df，包含 pseudotime 和 Stage 两列
# 1. 将 pseudotime 分成 50 个bin
cell_df <- cell_df %>%
  mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))

# 2. 统计每个 bin 内不同 Stage 细胞数量
bin_stage_counts <- cell_df %>%
  group_by(bin, Stage) %>%
  summarise(count = n(), .groups = 'drop')

# 3. 计算每个 bin 细胞总数
bin_totals <- bin_stage_counts %>%
  group_by(bin) %>%
  summarise(total = sum(count), .groups = 'drop')

# 4. 合并并计算比例
bin_stage_props <- bin_stage_counts %>%
  left_join(bin_totals, by = "bin") %>%
  mutate(prop = count / total)

# 5. 画堆积柱状图显示比例
pdf("Endosperm_trajectory_box1_umap-reanchor_plot3.pdf")
ggplot(bin_stage_props, aes(x = bin, y = prop, fill = Stage)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "Stage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

cell_df$predictedGroup_add <- paste(cell_df$predictedGroup_endo,cell_df$Stage,sep="_")
# 假设伪时间列是 pseudotime，先做bin
cell_df <- cell_df %>%
  mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))

# 统计每个bin中不同predictedGroup_add的数量
bin_group_counts <- cell_df %>%
  group_by(bin, predictedGroup_add) %>%
  summarise(count = n(), .groups = 'drop')

# 计算每个bin总数
bin_totals <- bin_group_counts %>%
  group_by(bin) %>%
  summarise(total = sum(count), .groups = 'drop')

# 计算比例
bin_group_props <- bin_group_counts %>%
  left_join(bin_totals, by = "bin") %>%
  mutate(prop = count / total)

# 画堆积百分比柱状图

pdf("Endosperm_trajectory_box2_predictedGroup_umap-reanchor_plot4.pdf")
ggplot(bin_group_props, aes(x = bin, y = prop, fill = predictedGroup_add)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "PredictedGroup_Stage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# 假设伪时间列是 pseudotime，先做bin
cell_df <- cell_df %>%
  mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))

# 统计每个bin中不同predictedGroup_add的数量
bin_group_counts <- cell_df %>%
  group_by(bin, predictedGroup_endo) %>%
  summarise(count = n(), .groups = 'drop')

# 计算每个bin总数
bin_totals <- bin_group_counts %>%
  group_by(bin) %>%
  summarise(total = sum(count), .groups = 'drop')

# 计算比例
bin_group_props <- bin_group_counts %>%
  left_join(bin_totals, by = "bin") %>%
  mutate(prop = count / total)

# 画堆积百分比柱状图
pdf("Endosperm_trajectory_box2_predictedGroup_endo_umap-reanchor_plot5.pdf")
ggplot(bin_group_props, aes(x = bin, y = prop, fill = predictedGroup_endo)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "PredictedGroup_Stage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#====================================================================================================================#
###                                plot trajectory heatmap - use old motif matrix                                  ###
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
motifs <- grep(pattern, rownames(traGSM), value = TRUE)
# motifs_z <- motifs[grep("^z:", motifs)]
trajGSM_sub <- trajGSM[motifs, ]
trajGIM_sub <- trajGIM[motifs, ]

mat <- assay(trajGSM_sub) 
p2_2 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)    
mat <- assay(trajGIM_sub)       
p3_2 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)       

pdf("EndospermTrajectory_umap_reanchor_plot6.pdf", width = 6, height = 8)
print(p1)
print(p2)
print(p2_2)
print(p3)
print(p3_2)
dev.off()

# combine motifMatrix and geneInteration matrix
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
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
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, force = TRUE)

# ComplexHeatmap::draw(ht1 + ht2)
pdf("EndospermTrajectory_umap_reanchor_plot7.pdf", width = 6, height = 8)
print(ht1)
print(ht2)
dev.off()

# combine motifMatrix and geneScore matrix
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
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
dev.off()

#====================================================================================================================#
###                  plot trajectory heatmap - use new motif matrix based on pseudotime trajectory                 ###
#====================================================================================================================#
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", threads = 1)
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "Stage", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addPeakMatrix(projAll)


# 添加motif
meme_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/Ath_TF_binding_motifs.meme" 
motifs <- read_meme(meme_file)
# 转换为PWMatrix格式
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
projAll <- addMotifAnnotations(ArchRProj = projAll, motifPWMs = pwmlist, name = "Motif")

# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)

# 保存
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor2_motif2.rds")
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
motifs <- grep(pattern, rownames(traGSM), value = TRUE)
# motifs_z <- motifs[grep("^z:", motifs)]
trajGSM_sub <- trajGSM[motifs, ]
trajGIM_sub <- trajGIM[motifs, ]

mat <- assay(trajGSM_sub) 
p2_2 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)    
mat <- assay(trajGIM_sub)       
p3_2 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)       

pdf("EndospermTrajectory_umap_reanchor_plot7.pdf", width = 6, height = 8)
print(p1)
print(p2)
print(p2_2)
print(p3)
print(p3_2)
dev.off()

# combine motifMatrix and geneInteration matrix
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
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
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
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

dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()


#====================================================================================================================#
###                                     RNA and ATAC data integration                                              ###
#====================================================================================================================#

############################## combine different stages snATAC data
arrow_files_all <- list.files(path = ".", pattern = "\\.arrow$", recursive = F, full.names = TRUE)
projAll <- ArchRProject(
   ArrowFiles = arrow_files_all,       
   geneAnnotation = gene_annotation,
   genomeAnnotation = genome_annotation,
   outputDirectory = "Stage_combined",   
   showLogo = F,           
   copyArrows = F,  
   threads = 8
 )
 
# Assume Sample is something like "A250415004_bent1"
projAll$Stage <- sub(".*_(globular|heart|torpedo|bent|cotyledon)[0-9]*", "\\1", projAll$Sample)
# table(projAll$Stage)
# bent	cotyledon	globular	heart	torpedo 
# 9342465	2452967	963371	896144	6958372

globular <- readRDS("rds/ArchR_globular_filterd.rds")
heart <- readRDS("rds/ArchR_heart_filterd.rds")
torpedo <- readRDS("rds/ArchR_torpedo_filterd.rds")
bent <- readRDS("rds/ArchR_bent_filterd.rds")
cotyledon <- readRDS("rds/ArchR_cotyledon_filterd.rds")

endo_cells_all <- c(globular$cellNames, heart$cellNames, torpedo$cellNames, bent$cellNames, cotyledon$cellNames)
# 2. 构建 projAll 后，筛选保留这些细胞
projAll <- subsetArchRProject(
  ArchRProj = projAll,
  cells = endo_cells_all,
  outputDirectory = "Stage_combined",
  dropCells = FALSE,
  force = TRUE)
projAll <- addGeneScoreMatrix(input = projAll, force = TRUE)

# Consider removing the A250410002_torpedo2 library
cells_to_keep <- projAll$cellNames[projAll$Sample != "A250410002_torpedo2"]
projAll <- projAll[cells_to_keep, ]
projAll <- addIterativeLSI(projAll, useMatrix = "GeneScoreMatrix", name = "IterativeLSI")
projAll <- addUMAP(projAll, reducedDims = "IterativeLSI", name = "UMAP")

saveRDS(projAll, "ArchR_Stage_combined_UMAP_rmTorpedo2.rds")

p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
# p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
pdf("ArchR_Stage_combined_UMAP_plot1.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, ncol = 3, top = "snATAC UMAP - Stage combined")
dev.off()
dev.off()
dev.off()


projAll <- readRDS("ArchR_Stage_combined_UMAP_rmTorpedo2.rds")
# use peak matrix
# 单独根据每一组细胞（例如聚类）鉴定peak，这样可以分析出不同组的特异性peak
# 拟混池重复
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", threads = 1)
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "Stage", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addPeakMatrix(projAll, force=TRUE)
# getAvailableMatrices(projAll)

# 保存
saveRDS(projAll, "ArchR_Stage_combined_UMAP_rmTorpedo2_peak.rds")

# 添加motif
meme_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/Ath_TF_binding_motifs.meme" 
motifs <- read_meme(meme_file)
corrected_pwms <- lapply(motifs, function(m) {
  matrix_data <- as.matrix(m@motif) 
  new("PWMatrix", ID = m@name, name = ifelse(length(m@altname) > 0, m@altname, m@name), profileMatrix = matrix_data, strand = "+")})
corrected_pwms <- lapply(corrected_pwms, function(m) {
  m@profileMatrix <- sweep(m@profileMatrix, 2, colSums(m@profileMatrix), "/")
  return(m)})
pwmlist <- do.call(PWMatrixList, corrected_pwms)
ids <- sapply(pwmlist, function(x) x@ID)
names(pwmlist) <- ids

projAll <- addMotifAnnotations(ArchRProj = projAll, motifPWMs = pwmlist, name = "Motif")

# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)

# 保存
projAll <- addIterativeLSI(projAll, useMatrix = "PeakMatrix", name = "IterativeLSI_peak")
projAll <- addUMAP(projAll, reducedDims = "IterativeLSI_peak", name = "UMAP_peak")


projAll <- addIterativeLSI(projAll, useMatrix = "MotifMatrix", name = "IterativeLSI_motif")
projAll <- addUMAP(projAll, reducedDims = "IterativeLSI_motif", name = "UMAP_motif")

saveRDS(projAll, "ArchR_Stage_combined_UMAP_rmTorpedo2_peak_motif.rds")

p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP_peak")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP_peak")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "nFrags", embedding = "UMAP_peak")
# p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP_peak")
pdf("ArchR_Stage_combined_UMAP_plot3.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, ncol = 3, top = "snATAC UMAP_peak - Stage combined")
dev.off()


p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP_motif")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP_motif")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "nFrags", embedding = "UMAP_motif")
# p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP_motif")
pdf("ArchR_Stage_combined_UMAP_plot4.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, ncol = 3, top = "snATAC _motif - Stage combined")
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()

########################## combine different stages snRNA data
rds_files <- c(
  "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/2.snRNA_results/20250430.final_cell_anno/04.72HAP.rds",
  "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/2.snRNA_results/20250430.final_cell_anno/05.heart.rds",
  "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/2.snRNA_results/20250430.final_cell_anno/06.torpedo.rds",
  "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/2.snRNA_results/20250430.final_cell_anno/07.bent.rds",
  "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/2.snRNA_results/20250430.final_cell_anno/08.Cot.rds"
)

seurat_list <- lapply(rds_files, readRDS)
stage_names <- c("globular", "heart", "torpedo", "bent", "cotyledon")
combined_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = stage_names,
  merge.data = TRUE
)

combined_seurat$Stage <- mapvalues(
  combined_seurat$Organ,
  from = c("cotyledon", "embryo_72h", "heart", "torpedo", "walking_stick"),
  to   = c("cotyledon", "globular",   "heart", "torpedo", "bent")
)

combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:30)

saveRDS(combined_seurat, "snRNA_stage_combined_UMAP_seurat.rds")

p1 <- DimPlot(combined_seurat, reduction = "pca", group.by = "Stage", label = TRUE) +  ggtitle("snRNA UMAP colored by stage")
p2 <- DimPlot(combined_seurat, reduction = "pca", group.by = "Cell_type", label = TRUE) +  ggtitle("snRNA UMAP colored by Cell_type")
p3 <- DimPlot(combined_seurat, reduction = "pca", group.by = "Group", label = TRUE) +  ggtitle("snRNA UMAP colored by Group")


pdf("ArchR_Stage_combined_UMAP_plot5.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, ncol = 3, top = "snRNA pca - Stage combined")
dev.off()

################################ integrate snRNA and snATAC data

combined_seurat <- readRDS("snRNA_stage_combined_UMAP_seurat.rds")
projAll <- readRDS("ArchR_Stage_combined_UMAP_rmTorpedo2.rds")
projAll <- addGeneScoreMatrix(input = projAll, force = TRUE)

geneScoreMat <- getMatrixFromProject(ArchRProj = projAll, useMatrix = "GeneScoreMatrix")
geneScoreMat_sparse <- assay(geneScoreMat)
rownames(geneScoreMat_sparse) <- rowData(geneScoreMat)$name
cells <- colnames(geneScoreMat_sparse)
pbmc.atac <- CreateSeuratObject(counts = geneScoreMat_sparse, assay = "ATAC")

pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)
pbmc.atac <- RunPCA(pbmc.atac, npcs = 30)
pbmc.atac <- RunSVD(pbmc.atac, assay = "ATAC", reduction.key = "LSI_", reduction.name = "lsi")
pbmc.atac <- FindNeighbors(pbmc.atac, reduction = "lsi", dims = 1:30)
pbmc.atac <- FindClusters(pbmc.atac, resolution = 8)
pbmc.atac <- RunUMAP(pbmc.atac,  reduction = "lsi",dims = 1:30)

# 识别anchor，使用全部的细胞和特征跑起来非常慢，随机抽取10000个细胞进行测试
set.seed(123) 
pbmc.rna <- subset(combined_seurat, cells = sample(colnames(combined_seurat), 20000))
pbmc.atac <- subset(pbmc.atac, cells = sample(colnames(pbmc.atac), 20000))
VariableFeatures(pbmc.rna) <- head(VariableFeatures(pbmc.rna), 1500)

transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), reference.assay = "RNA", query.assay = "ATAC", reduction = "cca")


# 使用TransferData函数进行数据转移映射
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$Cell_type, weight.reduction = pbmc.atac[["pca"]], dims = 2:30) 
# 添加预测信息 
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions) 
# 映射之后进行评估
# 方式一
pdf("ArchR_Stage_combined_UMAP_plot6.pdf",width = 8, height = 6)
hist(pbmc.atac$prediction.score.max) 
abline(v = 0.5, col = "red") 
dev.off()

table(pbmc.atac$prediction.score.max > 0.5) 

# 数据可视化: 画映射RNA注释之后的ATAC数据集
pbmc.atac.filtered <- subset(pbmc.atac, subset = prediction.score.max > 0.5) 
pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = unique(pbmc.rna$Cell_type))
pdf("ArchR_Stage_combined_UMAP_plot7.pdf",width = 8, height = 6)
DimPlot(pbmc.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + NoLegend() + scale_colour_hue(drop = FALSE) 
DimPlot(pbmc.rna, group.by = "Cell_type", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + NoLegend() 
dev.off()

# 共嵌入可视化
genes.use <- VariableFeatures(pbmc.rna) 
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ] 
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["pca"]], dims = 2:30) 
pbmc.atac[["RNA"]] <- imputation 
# 进行merge合并 
coembed <- merge(x = pbmc.rna, y = pbmc.atac) 
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE) 
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE) 
coembed <- RunUMAP(coembed, dims = 1:30) 
coembed$celltype <- ifelse(!is.na(coembed$Cell_type), coembed$Cell_type, coembed$predicted.id) 
coembed$tech <- ifelse(!is.na(coembed$Cell_type), "RNA", "ATAC") 

saveRDS(coembed, "snRNA_stage_combined_ATAC-RNA.rds")

pdf("ArchR_Stage_combined_UMAP_plot8.pdf",width = 10, height = 6)
DimPlot(coembed, group.by = "tech") 
DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE) 
DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() 
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()
dev.off()

pdf("RNA-ATAC-integration-test_cluster.pdf",width = 10, height = 6)
DimPlot(inter, group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
dev.off()

# 相关性可视化
predictions <- table(pbmc.atac$seurat_clusters, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

pdf("RNA-ATAC-integration-test6-lsi8.pdf",width = 10, height = 6)
ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# 整合ATAC注释的cluster的类型
cM <- as.matrix(confusionMatrix(pbmc.atac$seurat_clusters, pbmc.atac$predicted.id))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
match <- data.frame(celltype = preClust, cluster_number = rownames(cM), stringsAsFactors = FALSE)
# 把match表转成一个vector，方便后面查
cluster_to_celltype <- setNames(match$celltype, match$cluster_number)
# 给pbmc.atac加一列
pbmc.atac$cluster_celltype <- cluster_to_celltype[as.character(pbmc.atac$seurat_clusters)]
# 保存
saveRDS(pbmc.atac, "PCA_snATAC_root_combined.rds")

predictions <- table(pbmc.atac$atac_celltype, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)
# 细胞类型相关性可视化
pdf("RNA-ATAC-integration-test8-match.pdf",width = 10, height = 6)
ggplot(predictions, aes(Var2, Var1, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", low = "#ffffc8", high = "#7d0025") + xlab("Predicted cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()
row_percent <- cM / rowSums(cM)
# 找出最大值所在的列名和对应占比
result <- data.frame(RowName = rownames(row_percent), MaxColumn = colnames(row_percent)[apply(row_percent, 1, which.max)], MaxPercentage = apply(row_percent, 1, max))

# 画密度图
pbmc.atac$annotation_correct <- pbmc.atac$atac_celltype == pbmc.atac$predicted.id
correct <- length(which(pbmc.atac$atac_celltype == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$atac_celltype != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
pdf("RNA-ATAC-integration-test9.pdf",width = 10, height = 6)
ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) + geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct", labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct", labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
dev.off()

#画桑葚图
sankey_data <- data.frame(ATAC = pbmc.atac$atac_celltype, RNA = pbmc.atac$predicted.id)
sankey_data <- sankey_data[sankey_data$ATAC != "Unknown" & sankey_data$RNA != "Unknown", ]
sankey_data <- as.data.frame(table(sankey_data))
pdf("RNA-ATAC-integration-test10.pdf",width = 10, height = 6)
ggplot(sankey_data, aes(axis1 = ATAC, axis2 = RNA, y = Freq)) + geom_alluvium(aes(fill = ATAC), width = 1/12) + geom_stratum(width = 1/12, color = "grey") + geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + scale_x_discrete(limits = c("ATAC", "RNA"), expand = c(.05, .05)) + theme_minimal() + ggtitle("ATAC to RNA cell type mapping (Sankey Diagram)")
dev.off()



################################################################# foot print ############################################################
projAll <- readRDS("rds/ArchR_EndospermD_rmTorpedo_reAnchor2_motif2.rds")
projAll <- addImputeWeights(projAll)
#====================================================================================================================#
#                                               motif足迹分析                                                     ####
#====================================================================================================================#
motifPositions <- getPositions(projAll)
# 提取部分感兴趣的TF motifs用于展示
motifs <- c("AT1G75080")
# markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
# markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
# markerMotifs <- markerMotifs[c(1,2,3,14,15,16)]
# 为了准确找到TF足迹，我们需要大量的reads。因此，细胞需要进行分组生成拟混池ATAC-seq谱才能用于TF足迹分析。这些拟混池谱之前在peak鉴定时就已经保存为分组覆盖文件。 如果没有在ArchRProject添加分组覆盖信息，则运行如下命令
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage",force = TRUE)
# 通过positions参数来选择motif
seFoot <- getFootprints(ArchRProj = projAll, positions = motifPositions[motifs], groupBy = "Stage")
# 当我们获取了这些足迹，我们可以使用plotFootprints()函数进行展示。该函数能够同时以多种方式对足迹进行标准化。
# Tn5偏好的足迹标准化
# 减去Tn5偏好
pdf("Endosperm_footprint_plot4.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, plotName = "Footprints-Subtract-Bias-Stage", ArchRProj = projAll, normMethod = "Subtract", addDOC = FALSE, smoothWindow = 5)
dev.off()

# 这些图保存在ArchRProject的outputDirectory。如果需要绘制所有motif, 可以将其返回为ggplot2对象，需要注意这个ggplot对象会非常大
# 除以Tn5偏好
# pdf("Endosperm_footprint_plot5.pdf", width = 20, height = 20)
# plotFootprints(seFoot = seFoot, ArchRProj = projAll, normMethod = "Divide", plotName = "Footprints-Divide-Bias", addDOC = FALSE, smoothWindow = 5)
# dev.off()
# 无Tn5偏好标准化的足迹
# pdf("Endosperm_footprint_plot6.pdf", width = 20, height = 20)
# plotFootprints(seFoot = seFoot, ArchRProj = projAll,normMethod = "None", plotName = "Footprints-No-Normalization", addDOC = FALSE, smoothWindow = 5)
# dev.off()
# 与TSS的距离
# seTSS <- getFootprints(ArchRProj = projAll, positions = GRangesList(TSS = getTSS(projAll)), groupBy = "Stage", flank = 2000)
# plotFootprints(
#   seFoot = seTSS,
#   ArchRProj = projHeme5, 
#   normMethod = "None",
#   plotName = "TSS-No-Normalization",
#   addDOC = FALSE,
#   flank = 2000,
#   flankNorm = 100
# )

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
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor2_motif2_track.rds")

# 正向TF-调控因子鉴定
projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor2_motif2_track.rds")
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
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
p_GSM <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) + geom_point() + theme_ArchR() + geom_vline(xintercept = 0, lty = "dashed") + scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) + xlab("Correlation To Gene Score") + ylab("Max TF Motif Delta") + scale_y_continuous(expand = c(0,0), limits = c(0, max(corGSM_MM$maxDelta)*1.05)) + geom_text_repel(
    data = subset(corGSM_MM, TFRegulator == "YES"),
    aes(label = GeneScoreMatrix_name),
    color = "black",  # Label color
    size = 3,         # Text size
    max.overlaps = Inf
  )
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
p_GIM <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("NO" = "darkgrey", "YES" = "firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(corGIM_MM$maxDelta) * 1.05)) +
  # Add text labels only for points where TFRegulator == "YES"
  geom_text_repel(
    data = subset(corGIM_MM, TFRegulator == "YES"),
    aes(label = GeneIntegrationMatrix_name),
    color = "black",  # Label color
    size = 3,         # Text size
    max.overlaps = Inf
  )
  
pdf("Endosperm-Tracks-Marker-Genes-with-Peak2Gene-heatmap2.pdf",width = 15,height = 10)

print(p_GSM)
print(p_GIM)
dev.off()

projAll <- addImputeWeights(projAll)
p1 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "MotifMatrix", name = "z:AT1G75080", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))
p2 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "GeneScoreMatrix", name = "AT1G75080", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))
p3 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "GeneIntegrationMatrix", name = "AT1G75080", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))

pdf("AT1G75080_plot2.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, ncol = 3, top = "AT1G75080 combined")
dev.off()

