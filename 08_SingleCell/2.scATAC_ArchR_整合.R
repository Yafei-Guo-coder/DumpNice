library(matrixStats, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib")
packageVersion("matrixStats")
library(ArchR)
library(BSgenome.Athaliana.TAIR10.02022009)
library(parallel)
library(ggalt)
library(clustree)
library(dplyr)
library(Seurat)
library(patchwork)
library(gridExtra)
library(grid)
.libPaths("/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib")
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
library(viridis)
library(pheatmap)
library(Seurat)
library(Signac)
library(cowplot)
library(parallel)
library(ComplexHeatmap)

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 参数检查
if (length(args) < 4) {
  stop("Usage: Rscript integrate_by_group.R <input_directory> <TSS_enrich> <nFrags> <resolution>")
}

# 解析参数
input_dir <- args[1] 
TSS_enrich <- as.numeric(args[2])
nFrags <- as.numeric(args[3])
resolution <- as.numeric(args[4])

#====================================================================================================================#
#                                        创建基因组注释及ArchR对象                                                ####
#====================================================================================================================#
# 从BSgenome对象中创建基因组注释
genome_annotation <- createGenomeAnnotation(genome = BSgenome.Athaliana.TAIR10.02022009)
blacklist_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/BlackList_TAIR10.bed"
bl <- read.table(blacklist_file)
bl <- bl[,1:3]
colnames(bl) <- c("chr", "start", "end")
bl$strand <- "+" 
blacklist <- makeGRangesFromDataFrame(bl, keep.extra.columns = T)
genome_annotation$blacklist <- blacklist
# 定义工具脚本路径
gene_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.gene.txt"
trans_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.transcript.txt"
exon_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.exon.txt"
archr_utils <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/archr_utils.R"
# 加载工具脚本
source(archr_utils)
gene_annotation <- loadGeneinfo(gene = gene_file, exon = exon_file, trans = trans_file)
gene_annotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]] <- unlist(lapply(gene_annotation@listData[["exons"]]@elementMetadata@listData[["gene_id"]], function(x){unlist(strsplit(x, ":exon:"))[1]}))
gene_annotation@listData[["exons"]]@elementMetadata@listData[["symbol"]] <- unlist(lapply(gene_annotation@listData[["exons"]]@elementMetadata@listData[["symbol"]], function(x){unlist(strsplit(x, ":exon:"))[1]}))
rscripts <- list.files(path = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/ArchR_R", pattern = ".R$", recursive = F, full.names = F)
for (i in rscripts) {source(paste0("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/ArchR_R/",i))}
# 获取所有子目录名称
subdirs <- list.dirs(input_dir, recursive = F)
all_samples <- character()
all_fragments <- character()
for (dir in subdirs) {
  fragments <- list.files(path = dir, pattern = ".tsv.gz$", recursive = T, full.names = F)
  samples <- unlist(lapply(fragments, function(x){
    unlist(strsplit(x, ".", fixed = T))[1]
  }))
  all_samples <- c(all_samples, samples)
  all_fragments <- c(all_fragments, paste0(dir, "/", fragments))
}
cat("All samples:\n")
for (sample in all_samples) {cat(sample, "\n")}
cat("\nAll fragments:\n")
for (fragment in all_fragments) {cat(fragment, "\n")}
addArchRThreads(threads = 5)
# 创建 Arrow 文件
arrowFiles <- createArrowFiles(inputFiles = all_fragments,       
                               sampleNames = all_samples,        
                               outputNames = all_samples,       
                               minTSS = 0,                  
                               minFrags = 0,                 
                               excludeChr = c("ChrC", "ChrM"), 
                               addTileMat = T,             
                               addGeneScoreMat = T,       
                               subThreading = F,        
                               geneAnnotation = gene_annotation,
                               genomeAnnotation = genome_annotation) 

arrowFiles <- list.files(pattern = ".arrow", recursive = F, full.names = F)
# 创建一个ArchRProject对象
proj_all <- ArchRProject(
  ArrowFiles = arrowFiles,         
  geneAnnotation = gene_annotation, 
  genomeAnnotation = genome_annotation,
  outputDirectory = "seed",       
  showLogo = F,                   
  copyArrows = F, 
  threads = 8
)

projHeme1 <- readRDS("snATAC_raw.rds")

# r <- sample(projHeme1$cellNames, 160000) 
# projHeme1 <- projHeme1[which(projHeme1$cellNames %in% r), ]
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
df_list <- split(df, df$Sample)

plot_list <- lapply(names(df_list), function(sname) {
  sample_df <- df_list[[sname]]
  p <- ggPoint(
    x = sample_df[, "log10(nFrags)"],
    y = sample_df[, "TSSEnrichment"],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment"
  ) + ggtitle(paste("Sample:", sname))
  return(p)
})

# 额外QC图
p2 <- plotGroups(ArchRProj = projHeme1, groupBy = "Sample", colorBy = "cellColData",
                 name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p4 <- plotGroups(ArchRProj = projHeme1, groupBy = "Sample", colorBy = "cellColData",
                 name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

df2 <- as.data.frame(df)
colnames(df2)[1] <- "log10nFrags"

p5 <- ggplot(df2, aes(x = log10nFrags, y = Sample, fill = Sample)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "A") +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("") +
  ggtitle("Boxplot of log10(nFrags) by Sample")

p6 <- ggplot(df2, aes(x = TSSEnrichment, y = Sample, fill = Sample)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "A") +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("") +
  ggtitle("Boxplot of TSS Enrichment by Sample")

# TSS富集质量图
p7 <- plotFragmentSizes(ArchRProj = projHeme1)
p8 <- plotTSSEnrichment(ArchRProj = projHeme1)


pdf("QC_All_Plots_Annotated_raw.pdf", width = 6, height = 6)

for (i in seq_along(plot_list)) {
  grid::grid.newpage()
  grid::grid.text(paste0("Scatter: ", names(df_list)[i]), x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
  print(plot_list[[i]])
}

for (plt in list(p2, p4, p5, p6)) {
  grid::grid.newpage()
  grid::grid.text(plt$labels$title, x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
  print(plt)
}

# Fragment size and TSS enrichment
grid::grid.newpage()
grid::grid.text("Fragment Size Distribution", x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
print(p7)

grid::grid.newpage()
grid::grid.text("TSS Enrichment Profile", x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
print(p8)
dev.off()

saveRDS(proj_all, "snATAC_raw.rds")

# 过滤以及可视化
projHeme1 <- readRDS("snATAC_raw.rds")
idxPass <- which(projHeme1$TSSEnrichment >= TSS_enrich & projHeme1$nFrags>= nFrags)
cellsPass <- projHeme1$cellNames[idxPass]
projHeme1 <- projHeme1[cellsPass, ]
projHeme1 <- addDoubletScores(input = projHeme1, k = 10, knnMethod = "UMAP", LSIMethod = 1)
projHeme1 <- filterDoublets(projHeme1)

saveRDS(projHeme1, "snATAC_filterd.rds")

df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
df_list <- split(df, df$Sample)

plot_list <- lapply(names(df_list), function(sname) {
  sample_df <- df_list[[sname]]
  p <- ggPoint(
    x = sample_df[, "log10(nFrags)"],
    y = sample_df[, "TSSEnrichment"],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment"
  ) + ggtitle(paste("Sample:", sname))
  return(p)
})

# 额外QC图
p2 <- plotGroups(ArchRProj = projHeme1, groupBy = "Sample", colorBy = "cellColData",
                 name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p4 <- plotGroups(ArchRProj = projHeme1, groupBy = "Sample", colorBy = "cellColData",
                 name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

df2 <- as.data.frame(df)
colnames(df2)[1] <- "log10nFrags"

p5 <- ggplot(df2, aes(x = log10nFrags, y = Sample, fill = Sample)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "A") +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("") +
  ggtitle("Boxplot of log10(nFrags) by Sample")

p6 <- ggplot(df2, aes(x = TSSEnrichment, y = Sample, fill = Sample)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "A") +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("") +
  ggtitle("Boxplot of TSS Enrichment by Sample")

# TSS富集质量图
p7 <- plotFragmentSizes(ArchRProj = projHeme1)
p8 <- plotTSSEnrichment(ArchRProj = projHeme1)

# 所有图画入PDF，并添加标题注释
pdf("QC_All_Plots_Annotated_filter.pdf", width = 6, height = 6)

# 每个样本的散点图
for (i in seq_along(plot_list)) {
  grid::grid.newpage()
  grid::grid.text(paste0("Scatter: ", names(df_list)[i]), x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
  print(plot_list[[i]])
}

# Violin 和 boxplot 图
for (plt in list(p2, p4, p5, p6)) {
  grid::grid.newpage()
  grid::grid.text(plt$labels$title, x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
  print(plt)
}

# Fragment size and TSS enrichment
grid::grid.newpage()
grid::grid.text("Fragment Size Distribution", x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
print(p7)

grid::grid.newpage()
grid::grid.text("TSS Enrichment Profile", x = 0.5, y = 0.97, gp = gpar(fontsize = 14, fontface = "bold"))
print(p8)
dev.off()

#====================================================================================================================#
#                                           细胞聚类及批次效应消除                                                ####
#====================================================================================================================#

# projHeme1 <- readRDS("snATAC_filterd.rds")
all_samples <- unique(projHeme1$Sample)

# 从样本名中提取组织类型，例如 A250408001_heart1 -> heart
get_tissue_from_sample <- function(sample_name) {
  strsplit(sample_name, "_")[[1]][2] %>% gsub("[0-9]+$", "", .)
}

# 获取组织到样本的映射
sample_tissue_map <- setNames(sapply(all_samples, get_tissue_from_sample), all_samples)
tissue_groups <- split(names(sample_tissue_map), sample_tissue_map)
for (tissue_name in names(tissue_groups)) {
  cat("Processing tissue:", tissue_name, "\n")
  samples <- tissue_groups[[tissue_name]]
  
  idx <- which(projHeme1$Sample %in% samples)
  cells <- projHeme1$cellNames[idx]
  proj_sub <- projHeme1[cells, ]
  
  # LSI -> Harmony -> UMAP
  proj_sub <- addIterativeLSI(proj_sub, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2)
  proj_sub <- addHarmony(proj_sub, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample")
  proj_sub <- addUMAP(proj_sub, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30)
  proj_sub <- addClusters(proj_sub, reducedDims = "Harmony", name = "Clusters", resolution = resolution)
  # 绘图
  p1 <- plotEmbedding(proj_sub, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p2 <- plotEmbedding(proj_sub, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAPHarmony")
  p3 <- plotEmbedding(proj_sub, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
  
  # 保存结果和图
  saveRDS(proj_sub, paste0("ArchR_", tissue_name, "_harmony.rds"))
  pdf(paste0("UMAP_", tissue_name, "_Harmony.pdf"), width = 10, height = 6)
  grid.arrange(p1, p2, p3, ncol = 3, top = paste("UMAP -", tissue_name))
  dev.off()
}

#====================================================================================================================#
#                                            鉴定并测试marker基因                                                 ####
#====================================================================================================================#

# 自动识别 .rds 文件
files <- list.files(pattern = "^ArchR_.*_harmony\\.rds$")

for (f in files) {
  # 提取组织名
  tissue_name <- sub("ArchR_(.*)_harmony\\.rds", "\\1", f)
  cat("Processing Marker Genes for:", tissue_name, "\n")
  
  # 读取对象
  projHeme <- readRDS(f)
  
  # 计算 marker
  markersGS <- getMarkerFeatures(
    ArchRProj = projHeme,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  # 筛选 marker 基因，输出为列表
  markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
  
  # 保存每个 cluster 的 marker TSV
  for (i in seq_along(markerList)) {
    df <- as.data.frame(markerList[[i]])
    tsv_file <- paste0("markerList_", tissue_name, "_cluster", i, ".tsv")
    write.table(df, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  # 画热图（使用更严格的 cutoff）
  heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    transpose = TRUE
  )
  
  # 保存热图为 PDF
  pdf_file <- paste0("GeneScores_MarkerHeatmap_", tissue_name, ".pdf")
  plotPDF(heatmapGS,
          name = gsub(".pdf", "", pdf_file),
          width = 8,
          height = 6,
          ArchRProj = projHeme,
          addDOC = FALSE)
  dev.off()
}

