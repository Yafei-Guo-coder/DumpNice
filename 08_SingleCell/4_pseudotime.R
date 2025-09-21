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
filter_chr <- function(gr) {
  keepSeqlevels(gr, setdiff(seqlevels(gr), c("ChrC", "ChrM")), pruning.mode = "coarse")
}

# 对 genome_annotation 的每一部分进行清洗
genome_annotation_clean <- genome_annotation
genome_annotation_clean[[2]] <- filter_chr(genome_annotation[[2]])
blacklist <- genome_annotation[[3]]
blacklist_filtered <- keepSeqlevels(blacklist, paste0("Chr", 1:5), pruning.mode = "coarse")
genome_annotation_clean[[3]] <- blacklist_filtered

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

######################################################################################################### create archR project ############################################################################################

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

####################################################################################################################### motif #######################################################################################################

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

projAll <- readRDS("rds_Endosperm/ArchR_EndospermD_harmony_rmTorpedo2_motif.rds")
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", threads = 1)
projAll <- addGeneScoreMatrix(input = projAll, force = TRUE)
projAll <- addImputeWeights(projAll)

# anchor ##
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig4.endo_traj/Fig4a.endo_recluster/EndoCells_recluster.rds")
projAll2 <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = TRUE, force = TRUE, groupRNA = "endotype0514", nameCell = "predictedCell_endo", nameGroup = "predictedGroup_endo", nameScore = "predictedScore_endo", threads = 1)
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor2.rds")
# projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor.rds")

# anchor 3 ##
projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor2_motif.rds")
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
projAll2 <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = FALSE, force = TRUE, groupRNA = "endotype2025", nameCell = "predictedCell_endo2", nameGroup = "predictedGroup_endo2", nameScore = "predictedScore_endo2", threads = 1)
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor3_true_motif.rds")

# 对应上面的anchor 3, 采取另外一种方式，直接根据最匹配的细胞进行类型映射
# ArchR_EndospermD_rmTorpedo_reAnchor2_motif3.rds is created by 10_encoDevelop_motif.R in Rstudio
projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor2_motif3.rds")
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
map_vec <- seRNA@meta.data$endotype2025
names(map_vec) <- seRNA@meta.data$cellID
projAll@cellColData$endotype2025 <- map_vec[ projAll@cellColData$predictedCell_endo ]
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor3_motif3.rds")

# plot 看预测结果
p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedScore_endo", embedding = "UMAP")
#p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup_endo", embedding = "UMAP")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "endotype2025", embedding = "UMAP")

pdf("ArchR_EndospermD_UMAP_reanchor3_plot1.pdf", width = 13, height = 6)
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
pdf("Endosperm_trajectory_box1_umap-reanchor3_plot3.pdf")
ggplot(bin_stage_props, aes(x = bin, y = prop, fill = Stage)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "Stage") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

cell_df$predictedGroup_add <- paste(cell_df$endotype2025,cell_df$Stage,sep="_")
# cell_df$predictedGroup_add <- paste(cell_df$predictedGroup_endo,cell_df$Stage,sep="_")
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

pdf("Endosperm_trajectory_box2_predictedGroup_umap-reanchor3_plot4.pdf")
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

################################ plot trajectory heatmap - use old motif matrix #####################################################
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
p2_2 <- pheatmap(scaled_mat, color = paletteContinuous(set = "horizonExtra"),
                 cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)
mat <- assay(trajGIM_sub)
scaled_mat <- t(scale(t(mat)))
p3_2 <- pheatmap(scaled_mat, color = paletteContinuous(set = "horizonExtra"),
                 cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)

pdf("EndospermTrajectory_umap_reanchor_plot6.pdf", width = 6, height = 8)
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
#

################################ plot trajectory heatmap - use new motif matrix based on pseudotime trajectory #################################

projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", force=TRUE, threads = 1)
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "Stage", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addPeakMatrix(projAll,force = TRUE)


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
projAll <- addMotifAnnotations(ArchRProj = projAll, force = TRUE,motifPWMs = pwmlist, name = "Motif")

# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll,force = TRUE)
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
motifs <- grep(pattern, rownames(trajGSM), value = TRUE)
# motifs_z <- motifs[grep("^z:", motifs)]
trajGSM_sub <- trajGSM[motifs, ]
trajGIM_sub <- trajGIM[motifs, ]
mat <- assay(trajGSM_sub)
scaled_mat <- t(scale(t(mat)))
p2_2 <- pheatmap(scaled_mat, color = paletteContinuous(set = "horizonExtra"),
                 cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)
mat <- assay(trajGIM_sub)
scaled_mat <- t(scale(t(mat)))
p3_2 <- pheatmap(scaled_mat, color = paletteContinuous(set = "horizonExtra"),
                 cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)

pdf("EndospermTrajectory_umap_reanchor_plot7.pdf", width = 6, height = 8)
print(p1)
print(p2)
print(p2_2)
print(p3)
print(p3_2)
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


# 提取鱼雷胚的胚乳细胞
proj <- readRDS("rds_Endosperm/ArchR_EndospermD_rmTorpedo_reAnchor2_motif2.rds")
# 过滤
idxPass <- which(proj$Sample %in% c("A250410001_torpedo1","A250410003_torpedo3"))
cellsPass <- proj$cellNames[idxPass]
projAll <- proj[cellsPass, ]

