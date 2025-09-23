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
# 对应上面的 anchor2, 因为上一版anchor2 没有剔除snRNA的多余时期(2_28h|3_48h)的细胞类型, 采取另外一种方式，直接根据最匹配的细胞进行类型映射
# ArchR_EndospermD_rmTorpedo_reAnchor2_motif3.rds is created by 10_encoDevelop_motif.R in Rstudio

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
  outputDirectory = "Merged_Endosperm_0911",   
  showLogo = F,           
  copyArrows = T,  
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
  outputDirectory = "Merged_Endosperm_Filtered_0911",
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
saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/ED_0911/ArchR_EndospermD_UMAP_rmTorpedo2.rds")

projAll <- addIterativeLSI(projAll, useMatrix = "MotifMatrix", name = "IterativeLSI")

projAll <- addUMAP(projAll, reducedDims = "IterativeLSI", name = "UMAP")
projAll <- addTSNE(ArchRProj = projAll, reducedDims = "IterativeLSI", name = "TSNE", perplexity = 30)
p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Clusters_merge", embedding = "UMAP")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/ED_0911/pdf/ArchR_EndospermD_UMAP.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "UMAP - Endosperm Development")
dev.off()

p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "TSNE")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "TSNE")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Clusters_merge", embedding = "TSNE")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup", embedding = "TSNE")
pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/ED_0911/pdf/ArchR_EndospermD_TSNE.pdf", width = 13, height = 6)
grid.arrange(p1, p2, p3, p4, ncol = 4, top = "UMAP - Endosperm Development")
dev.off()

#====================================================================================================================#
#                                                        motif 富集                                                   #
#====================================================================================================================#
# # 计算motif富集矩阵，然后构建LSI，随后做伪时序分析
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2.rds")
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", threads = 1)
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
# projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "predictedGroup", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "Stage", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addPeakMatrix(projAll)

# 添加motif
# projHeme5 <- readRDS("ArchR_EndospermD_harmony_rmTorpedo2_peak.rds")
meme_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/motif/Ath_TF_binding_motifs_addV2.meme" 
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
projAll <- addMotifAnnotations(ArchRProj = projAll, motifPWMs = pwmlist, name = "Motif", force = TRUE)

# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)

# 保存
saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_motif_stage.rds")


#====================================================================================================================#
#                                         使用ArchR进行RNA和ATAC数据映射                                                #
#====================================================================================================================#
# 使用更新过的snRNA数据的注释进行ATAC数据细胞类型映射
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds_Endosperm/ArchR_EndospermD_harmony_rmTorpedo2_motif_V2.rds")
projAll <- addImputeWeights(projAll)
# new RNA annotation
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
cell_ids <- seRNA@meta.data$cellID
cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
seRNA_filtered <- subset(seRNA, cells = cells_to_keep)
# old RNA annotation
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig4.endo_traj/Fig4a.endo_recluster/EndoCells_recluster.rds")
cell_ids <- seRNA@meta.data$cellID

# old
cell_ids <- rownames(seRNA@meta.data)
cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
seRNA_filtered <- subset(seRNA, cells = cells_to_keep)


# 无约束整合
# new
projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = TRUE, force = TRUE, groupRNA = "endotype2025", nameCell = "predictedCell_endo2025", nameGroup = "predictedGroup_endo2025", nameScore = "predictedScore_endo2025", threads = 1)
# old
projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = FALSE, force = TRUE, groupRNA = "endotype0514", nameCell = "predictedCell_endo2025", nameGroup = "predictedGroup_endo2025", nameScore = "predictedScore_endo2025", threads = 1)

# 约束整合
# new
groupList <- SimpleList(
    globular = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "globular"],
        RNA = seRNA_filtered$cellID[seRNA_filtered$stage == "04.72HAP"]
    ),
    heart = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "heart"],
        RNA = seRNA_filtered$cellID[seRNA_filtered$stage == "05.heart"]
    ),
    torpedo = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "torpedo"],
        RNA = seRNA_filtered$cellID[seRNA_filtered$stage == "06.torpedo"]
    ),
    bent = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "bent"],
        RNA = seRNA_filtered$cellID[seRNA_filtered$stage == "07.bent"]
    ),
    cotyledon = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "cotyledon"],
        RNA = seRNA_filtered$cellID[seRNA_filtered$stage == "08.Cot"]
    )
)

# old
groupList <- SimpleList(
    globular = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "globular"],
        RNA = rownames(seRNA_filtered@meta.data)[seRNA_filtered$stage == "04.72HAP"]
    ),
    heart = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "heart"],
        RNA = rownames(seRNA_filtered@meta.data)[seRNA_filtered$stage == "05.heart"]
    ),
    torpedo = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "torpedo"],
        RNA = rownames(seRNA_filtered@meta.data)[seRNA_filtered$stage == "06.torpedo"]
    ),
    bent = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "bent"],
        RNA = rownames(seRNA_filtered@meta.data)[seRNA_filtered$stage == "07.bent"]
    ),
    cotyledon = SimpleList(
        ATAC = projAll$cellNames[projAll$Stage == "cotyledon"],
        RNA = rownames(seRNA_filtered@meta.data)[seRNA_filtered$stage == "08.Cot"]
    )
)

# new
projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA_filtered, addToArrow = FALSE, groupList = groupList, groupRNA = "endotype2025", nameCell = "predictedCell_Co", nameGroup = "predictedGroup_Co", k.anchor = 50, nameScore = "predictedScore_Co", threads = 5)

# old
projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA_filtered, addToArrow = FALSE, groupList = groupList, groupRNA = "endotype0514", nameCell = "predictedCell_Co", nameGroup = "predictedGroup_Co", nameScore = "predictedScore_Co", threads = 5)

saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2.rds")

# # > table(projAll$predictedGroup_endo2025)
# #     CZE     ESR ESR_PEN     MCE     PEN 
# #    1535    4287    1495       2    5870 

#====================================================================================================================#
#                                         使用Signac进行RNA和ATAC数据映射                                               #
#====================================================================================================================#
# # 换一种方式，使用Signac进行映射
# # 处理RNA数据
# seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
# cell_ids <- seRNA@meta.data$cellID
# cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
# pbmc.rna <- subset(seRNA, cells = cells_to_keep)
# pbmc.rna <- NormalizeData(pbmc.rna)
# pbmc.rna <- FindVariableFeatures(pbmc.rna)
# pbmc.rna <- ScaleData(pbmc.rna)
# pbmc.rna <- RunPCA(pbmc.rna)
# pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

# # 处理atac数据
# projHeme2 <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds_Endosperm/ArchR_EndospermD_harmony_rmTorpedo2_motif_V2.rds")
# geneScoreMat <- getMatrixFromProject(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix")
# geneScoreMat_sparse <- assay(geneScoreMat)
# rownames(geneScoreMat_sparse) <- rowData(geneScoreMat)$name
# cells <- colnames(geneScoreMat_sparse)
# pbmc.atac <- CreateSeuratObject(counts = geneScoreMat_sparse, assay = "ATAC")

# pbmc.atac <- NormalizeData(pbmc.atac)
# pbmc.atac <- FindVariableFeatures(pbmc.atac)
# pbmc.atac <- ScaleData(pbmc.atac)
# pbmc.atac <- RunPCA(pbmc.atac, npcs = 30)
# # 使用Lsi聚类
# pbmc.atac <- RunSVD(pbmc.atac, assay = "ATAC", reduction.key = "LSI_", reduction.name = "lsi")
# pbmc.atac <- FindNeighbors(pbmc.atac, reduction = "lsi", dims = 1:30)
# pbmc.atac <- FindClusters(pbmc.atac, resolution = 8)
# pbmc.atac <- RunUMAP(pbmc.atac,  reduction = "lsi",dims = 1:30)
# # 在最后画图

# # 取子集计算（需要确认一下细胞的数量）
# DefaultAssay(pbmc.rna) <- "RNA"
# DefaultAssay(pbmc.atac) <- "ATAC"

# pbmc.reference.small <- subset(pbmc.rna, cells = sample(colnames(pbmc.rna), 2000))
# pbmc.reference.small <- NormalizeData(pbmc.reference.small)
# pbmc.reference.small <- FindVariableFeatures(pbmc.reference.small)
# pbmc.reference.small <- ScaleData(pbmc.reference.small)
# pbmc.reference.small <- RunPCA(pbmc.reference.small)
# pbmc.query.small <- subset(pbmc.atac, cells = sample(colnames(pbmc.atac), 2000))
# pbmc.query.small <- RunTFIDF(pbmc.query.small)
# pbmc.query.small <- FindTopFeatures(pbmc.query.small, min.cutoff = 'q0')
# pbmc.query.small <- RunSVD(pbmc.query.small)

# pbmc.reference.small <- pbmc.rna
# pbmc.query.small <- pbmc.atac
p1 <- DimPlot(pbmc.rna, group.by = "endotype2025", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
pdf("anchor_ED/pdf/0914_test_plot1.pdf")
p1 + p2
dev.off()
# # anchors
# transfer.anchors <- FindTransferAnchors(reference = pbmc.reference.small, query = pbmc.query.small, k.anchor = 50, features = VariableFeatures(pbmc.reference.small), reference.assay = "RNA", query.assay = "ATAC", reduction = "cca")
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), k.anchor = 50,
    reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
# # Transfer labels
# predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = pbmc.reference.small$endotype2025, weight.reduction = pbmc.query.small[['lsi']], dims = 2:30)
# atac <- AddMetaData(pbmc.query.small, metadata = predicted.labels)
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$sendotype2025,
    weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]],
    dims = 2:30)
pbmc.atac[["RNA"]] <- imputation

coembed <- merge(x = pbmc.rna, y = pbmc.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$endotype2025), coembed$endotype2025, coembed$predicted.id) 
coembed$tech <- ifelse(!is.na(coembed$endotype2025), "RNA", "ATAC") 

#coembed <- addHarmony(ArchRProj = coembed, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "assay")
pdf("anchor_ED/pdf/0914_test5.pdf", width = 20, height = 6)
DimPlot(coembed, group.by = c("celltype","tech","endotype2025","predicted.id"))
dev.off()
dev.off()
dev.off()
# # 共嵌入可视化
# genes.use <- VariableFeatures(pbmc.reference.small)
# refdata <- GetAssayData(pbmc.reference.small, assay = "RNA", slot = "data")[genes.use, ] 
# imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]], dims = 2:30) 
# atac[["RNA"]] <- imputation 

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$seurat_annotations,
    weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)

# # 进行merge合并 
# coembed <- merge(x = pbmc.reference.small, y = atac) 
# coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE) 
# coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE) 
# coembed <- RunUMAP(coembed, dims = 1:30) 
# coembed$celltype <- ifelse(!is.na(coembed$endotype2025), coembed$endotype2025, coembed$predicted.id) 
# coembed$tech <- ifelse(!is.na(coembed$endotype2025), "RNA", "ATAC") 
# saveRDS(atac, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_signac.rds")
# saveRDS(coembed, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_coembed2.rds")

# # 提取锚点的配对信息
# anchor_pairs <- transfer.anchors@anchors
# # 保存为 txt 或 csv
# write.table(anchor_pairs, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/transfer_anchors_ED2.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# # 画图
# pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/pdf/RNA-ATAC-integration-ED_develop2.pdf",width = 10, height = 6)
# DimPlot(pbmc.rna, group.by = "endotype2025", label = TRUE) + NoLegend() + ggtitle("RNA")
# DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
# DimPlot(coembed, group.by = "tech") 
# DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE) 
# DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() 
# dev.off()
# dev.off()
# dev.off()
# dev.off()
# dev.off()

  proj_sub <- addIterativeLSI(proj_sub, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2)
  proj_sub <- addHarmony(proj_sub, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample")
  proj_sub <- addUMAP(proj_sub, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30)
  proj_sub <- addClusters(proj_sub, reducedDims = "Harmony", name = "Clusters", resolution = resolution)
  p1 <- plotEmbedding(proj_sub, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p2 <- plotEmbedding(proj_sub, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAPHarmony")
  p3 <- plotEmbedding(proj_sub, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
  
  saveRDS(proj_sub, paste0("ArchR_", tissue_name, "_harmony.rds"))
  pdf(paste0("UMAP_", tissue_name, "_Harmony.pdf"), width = 10, height = 6)
  grid.arrange(p1, p2, p3, ncol = 3, top = paste("UMAP -", tissue_name))
  dev.off()

#====================================================================================================================#
#                                                 计算细胞类型相关性                                                    #
#====================================================================================================================#
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2.rds")
cells.keep <- rownames(projAll)[projAll$predictedScore_endo2025 > 0.6]
projAll <- projAll[cells.keep, ] 


seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
cell_ids <- seRNA@meta.data$cellID
cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
seRNA_filtered <- subset(seRNA, cells = cells_to_keep)

gim <- getMatrixFromProject(projAll, useMatrix = "GeneIntegrationMatrix")
gim_genes <- rowData(gim)$name
gene_scores <- assay(gim, "GeneIntegrationMatrix")
rownames(gene_scores) <- rowData(gim)$name

# 方式一，RNA高变基因以及跟ATAC overlap的基因，采用这种方法
marker <- VariableFeatures(seRNA_filtered)

# # 方式二，找RNA差异表达的基因，效果不好，不采用这种方法
# Idents(seRNA_filtered) <- seRNA_filtered$endotype2025
# seRNA_filtered.markers <- FindAllMarkers(seRNA_filtered, only.pos = TRUE)
# marker <- seRNA_filtered.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
# # 找差异开放的基因
# projAll$endotype2025_merged <- projAll$endotype2025
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
genes_of_interest <- head(genes_of_interest, 1000)

# RNA 原始数据（可以用 scale.data 或 data）
rna_mat <- GetAssayData(seRNA_filtered, assay = "RNA", slot = "data")[genes_of_interest, ]

# 可以按细胞类型汇总平均表达
rna_avg <- aggregate(t(as.matrix(rna_mat)), by = list(celltype = seRNA_filtered$endotype2025), FUN = mean)
rownames(rna_avg) <- rna_avg$celltype
rna_avg <- t(rna_avg[,-1])

# 提取基因和细胞
atac_mat <- gene_scores[genes_of_interest, ]

# 按 ATAC 细胞类型平均
atac_avg <- aggregate(t(as.matrix(atac_mat)), by = list(celltype = projAll$predictedGroup_Co), FUN = mean)
rownames(atac_avg) <- atac_avg$celltype
atac_avg <- t(atac_avg[,-1])
cor_matrix1 <- cor(rna_avg, atac_avg, method = "pearson")  
# 行：RNA 细胞类型，列：ATAC 细胞类型

row_order <- c("CZE","ESR","PEN","ESR_PEN","CZE_PEN","MCE","MCE-like")
col_order <- c("CZE","ESR","PEN","ESR_PEN")

cor_matrix_ordered <- cor_matrix1[row_order, col_order]
# cor_matrix_ordered2 <- cor_matrix[row_order, col_order]
# my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pdf("anchor_stage/pdf/celltype_correlation_heatmap_co.pdf")
pheatmap(cor_matrix1, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, main = "RNA vs ATAC cell type correlation")
dev.off()

#====================================================================================================================#
#                                                         构建伪时序                                                   #
#====================================================================================================================#
# projAll <- addIterativeLSI(projAll, useMatrix = "PeakMatrix", name = "IterativeLSI",force = TRUE)
# projAll <- addUMAP(projAll, reducedDims = "IterativeLSI", name = "UMAP",force = TRUE)
# projAll <- addTSNE(ArchRProj = projAll, reducedDims = "IterativeLSI", name = "TSNE", perplexity = 30)
# saveRDS(projAll, "ArchR_EndospermD_UMAP_rmTorpedo2.rds")

# plot 看预测结果
p1 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Stage", embedding = "UMAP")
p2 <- plotEmbedding(projAll, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP")
p3 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAP")
p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "predictedScore_Co", embedding = "UMAP")
p5 <- plotEmbedding(projAll, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
# p4 <- plotEmbedding(projAll, colorBy = "cellColData", name = "endotype2025", embedding = "UMAP")

# cells_back <- read.table("cells_bent_PEN.txt", stringsAsFactors = FALSE, comment.char = "")[,1]
# projAll$highlight <- ifelse(projAll$cellNames %in% cells_back, "bent_PEN", "other")
p6 <- plotEmbedding(projAll, colorBy = "cellColData", name = "highlight", embedding = "UMAP")
pdf("anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_predictedScore_Co_subcell_k50_combined_cluster.pdf", width = 16, height = 6)
grid.arrange(p1, p2, p3, p4, p5, ncol = 5, top = "UMAP - Endosperm Development")
#grid.arrange(p1, p2, ncol = 2, top = "UMAP - Endosperm Development")
dev.off()
dev.off()

pdf("anchor_ED/pdf/bent_PEN_ED.pdf")
grid.arrange(p6,top = "UMAP - Endosperm Development")
#grid.arrange(p1, p2, ncol = 2, top = "UMAP - Endosperm Development")
dev.off()
dev.off()

# 后面的伪时序有点异常，分析原因是映射的噪音造成的，采取过滤掉低质量的映射细胞来处理
idxPass <- which(projAll$predictedScore_Co > 0.6)
cellsPass <- projAll$cellNames[idxPass]
projAll <- projAll[cellsPass, ]

# 查看各 cluster 的时间分布，辅助选择轨迹顺序
table(projAll$Stage, projAll$predictedScore_Co)
# 假设经验上知道 trajectory 顺序是以下几个 cluster
trajectory_order <- c("globular", "heart", "torpedo", "bent", "cotyledon")
projAll <- addTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", groupBy = "Stage", trajectory = trajectory_order,  preFilterQuantile = 0.5, postFilterQuantile = 0.8, embedding = "UMAP", force = TRUE)
projAll <- addTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", groupBy = "Stage", trajectory = trajectory_order,  preFilterQuantile = 0.9,  embedding = "UMAP", force = TRUE)
saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_trajectory.rds")

pdf("anchor_ED/pdf/ArchR_EndospermD_UMAP_reanchor2_predictedScore_Co_subcell_traj_combined.pdf", width = 6, height = 6)
p <- plotTrajectory(projAll, trajectory = "EndospermTrajectory", colorBy = "cellColData", name = "EndospermTrajectory")
print(p)
dev.off()

# 结合伪时间和真实细胞
cell_df <- as.data.frame(getCellColData(projAll))
# 假设你的数据框是 cell_df，包含 pseudotime 和 Stage 两列
# 1. 将 pseudotime 分成 50 个bin
cell_df <- cell_df %>% mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))
# 2. 统计每个 bin 内不同 Stage 细胞数量
bin_stage_counts <- cell_df %>% group_by(bin, Stage) %>% summarise(count = n(), .groups = 'drop')
# 3. 计算每个 bin 细胞总数
bin_totals <- bin_stage_counts %>% group_by(bin) %>% summarise(total = sum(count), .groups = 'drop')
# 4. 合并并计算比例
bin_stage_props <- bin_stage_counts %>% left_join(bin_totals, by = "bin") %>% mutate(prop = count / total)
# 5. 画堆积柱状图显示比例
pdf("anchor_ED/pdf/Endosperm_trajectory_box1_umap-reanchor2_subcell_predictedGroup_Co_combined.pdf")
ggplot(bin_stage_props, aes(x = bin, y = prop, fill = Stage)) + geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "Stage") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

cell_df$predictedGroup_add <- paste(cell_df$predictedScore_Co,cell_df$Stage,sep="_")
# cell_df$predictedGroup_add <- paste(cell_df$predictedGroup_endo,cell_df$Stage,sep="_")
# 假设伪时间列是 pseudotime，先做bin
cell_df <- cell_df %>% mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))
# 统计每个bin中不同predictedGroup_add的数量
bin_group_counts <- cell_df %>% group_by(bin, predictedGroup_add) %>% summarise(count = n(), .groups = 'drop')
# 计算每个bin总数
bin_totals <- bin_group_counts %>% group_by(bin) %>% summarise(total = sum(count), .groups = 'drop')
# 计算比例
bin_group_props <- bin_group_counts %>% left_join(bin_totals, by = "bin") %>% mutate(prop = count / total)
# 画堆积百分比柱状图
pdf("anchor_ED/pdf/Endosperm_trajectory_box2_predictedGroup_umap-reanchor2_subcell_predictedGroup_Co_combined.pdf")
ggplot(bin_group_props, aes(x = bin, y = prop, fill = predictedGroup_add)) + geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "PredictedGroup_Stage") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# 假设伪时间列是 pseudotime，先做bin
cell_df <- cell_df %>% mutate(bin = cut(EndospermTrajectory, breaks = 50, include.lowest = TRUE))
# 统计每个bin中不同predictedGroup_add的数量
bin_group_counts <- cell_df %>% group_by(bin, predictedScore_Co) %>% summarise(count = n(), .groups = 'drop')
# 计算每个bin总数
bin_totals <- bin_group_counts %>% group_by(bin) %>% summarise(total = sum(count), .groups = 'drop')
# 计算比例
bin_group_props <- bin_group_counts %>% left_join(bin_totals, by = "bin") %>% mutate(prop = count / total)
# 画堆积百分比柱状图
pdf("pdf/Endosperm_trajectory_box3_predictedGroup_endo_umap-reanchor2_subcell_predictedGroup_Co_combined.pdf")
ggplot(bin_group_props, aes(x = bin, y = prop, fill = predictedScore_Co)) + geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "Pseudotime bin", y = "Proportion of cells", fill = "PredictedGroup_Stage") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
p2_2 <- pheatmap(scaled_mat, color = paletteContinuous(set = "horizonExtra"), cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)
mat <- assay(trajGIM_sub)
scaled_mat <- t(scale(t(mat)))
p3_2 <- pheatmap(scaled_mat, color = paletteContinuous(set = "horizonExtra"), cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)

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
projAll <- addPeakMatrix(projAll, force = TRUE)

# 添加 motif
meme_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/motif/Ath_TF_binding_motifs_add.meme" 
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
projAll <- addMotifAnnotations(ArchRProj = projAll, force = TRUE,motifPWMs = pwmlist, name = "Motif")
# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll,force = TRUE)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)

# 保存
# saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_trajectory_motif.rds")
gene <- c("AT1G71250","AT1G21970","AT1G08560","AT1G30690","AT1G48930","AT3G62830","AT5G13170","AT3G08900","AT1G68560","AT2G47650","AT2G25540","AT4G33600","AT3G19820","AT4G04460","AT1G62340","AT5G50750","AT3G48740","AT5G50260","AT2G19500","AT1G23320","AT2G35670","AT5G60440","AT2G35230","AT1G65330","AT1G67775","AT1G50650","AT1G62340","AT3G08900","AT2G47650","AT1G49770","AT1G68560","AT3G11520","AT1G48930","AT1G75080","AT5G07280","AT1G71250","AT3G26744","AT3G48740","AT1G11190","AT5G18270")

trajectory_order <- c("globular", "heart", "torpedo", "bent", "cotyledon")
projAll <- addTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", groupBy = "Stage", trajectory = trajectory_order, postFilterQuantile = 0.9, embedding = "UMAP", force = TRUE)

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
#
################################################################# foot print ############################################################
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_trajectory_motif.rds")

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
pdf("pdf/Endosperm_footprint_Subtract.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, plotName = "Footprints-Subtract-Bias", ArchRProj = projAll, normMethod = "Subtract", addDOC = FALSE, smoothWindow = 5)
dev.off()
# 这些图保存在ArchRProject的outputDirectory。如果需要绘制所有motif, 可以将其返回为ggplot2对象，需要注意这个ggplot对象会非常大
# 除以Tn5偏好
pdf("pdf/Endosperm_footprint_Divide.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, ArchRProj = projAll, normMethod = "Divide", plotName = "Footprints-Divide-Bias", addDOC = FALSE, smoothWindow = 5)
dev.off()
# 无Tn5偏好标准化的足迹
pdf("pdf/Endosperm_footprint_None.pdf", width = 20, height = 20)
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

#====================================================================================================================#
#                                          正向TF-调控因子鉴定                                                       ####
#====================================================================================================================#
projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_trajectory_motif.rds")
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
library(ggrepel)
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
  
pdf("anchor_ED/pdf/Endosperm-Tracks-Marker-Genes-with-Peak2Gene-heatmap2_predictedGroup_Co_k50_09.pdf",width = 15,height = 10)
print(p_GSM)
print(p_GIM)
dev.off()
