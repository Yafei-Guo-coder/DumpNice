#============================================================================================#
#                                           加载包                                         ####
#====================================================================================================================#
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
#devtools::install_version("matrixStats", version="1.5.0", lib = "./lib")
#install.packages("universalmotif", lib = "./lib")
library(universalmotif)
library(tidydr)
library(ggplot2)
library(chromVAR)
library(TFBSTools)
library(SummarizedExperiment)
library(motifmatchr)
library(Matrix)
library(rhdf5)
library(readxl)
#options(mc.cores = 1) 
#====================================================================================================================#
#                                    读root, silique数据                                   ####
#====================================================================================================================#

silique <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/rds/PCA_snATAC_motifs_silique.rds")
root <- readRDS("PCA_snATAC-Integration-root-0316-1.8.rds")
arrows_root <- getArrowFiles(root)
arrows_silique <- getArrowFiles(silique)

#添加基因组
genome_annotation <- createGenomeAnnotation(genome = BSgenome.Athaliana.TAIR10.02022009)
blacklist_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/BlackList_TAIR10.bed"
bl <- read.table(blacklist_file)[,1:3]
colnames(bl) <- c("chr", "start", "end")
bl$strand <- "*"  # 更推荐用 "*" 而不是 "+"
blacklist <- makeGRangesFromDataFrame(bl, keep.extra.columns = TRUE)
genome_annotation$blacklist <- blacklist

# 定义工具脚本路径
gene_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.gene.txt"
trans_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.transcript.txt"
exon_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.exon.txt"
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
for (i in rscripts) { source(paste0("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/ArchR_R/",i))}

# 获取所有子目录名称
subdirs <- list.dirs("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202407.At_PCA/0.data/snATAC/fragments", recursive = F)

# 初始化空向量用于存储所有样本名称和片段文件路径
all_samples <- character()
all_fragments <- character()

# 遍历每个子目录
for (dir in subdirs) {
  # 列出当前子目录下所有以 ".tsv.gz" 结尾的片段文件（不递归）
  fragments <- list.files(path = dir, pattern = ".tsv.gz$", recursive = F, full.names = F)
  # 提取每个片段文件的样本名称（假设样本名称是文件名的第一部分）
  samples <- unlist(lapply(fragments, function(x){
    unlist(strsplit(x, ".", fixed = T))[1]}))
  # 将样本名称和片段文件路径添加到总向量中
  all_samples <- c(all_samples, samples)
  all_fragments <- c(all_fragments, paste0(dir, "/", fragments))
}

cat("All samples:\n")
for (sample in all_samples) {cat(sample, "\n")}

cat("\nAll fragments:\n")
for (fragment in all_fragments) {cat(fragment, "\n")}

proj_combined <- ArchRProject(ArrowFiles = c(arrows_root, arrows_silique), geneAnnotation = gene_annotation, genomeAnnotation = genome_annotation, outputDirectory = "root-silique", showLogo = F, copyArrows = F, threads = 8)
proj_combined$organ <- ifelse(grepl("Silique", proj_combined$cellNames), "Silique", "Root")
proj_combined <- addIterativeLSI(ArchRProj = proj_combined, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2)
proj_combined <- addHarmony(ArchRProj = proj_combined, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "organ")
proj_combined <- addUMAP(ArchRProj = proj_combined, reducedDims = "Harmony", name = "UMAP")

# 从 root 和 silique 提取 cell_type 和 cellNames
root_meta <- data.frame(cellNames = root$cellNames, cell_type = root$predictedGroup_1.8)
silique_meta <- data.frame(cellNames = silique$cellNames, cell_type = silique$predictedGroup_Un)
# 合并两份 metadata
cell_type_meta <- rbind(root_meta, silique_meta)
# 确保 rownames 是 cellNames
rownames(cell_type_meta) <- cell_type_meta$cellNames
# 对 proj_combined 的细胞做匹配
common_cells <- intersect(proj_combined$cellNames, rownames(cell_type_meta))
# 构建 cell_type 列（未匹配的为 NA）
cell_type_full <- rep(NA, length(proj_combined$cellNames))
names(cell_type_full) <- proj_combined$cellNames
cell_type_full[common_cells] <- cell_type_meta[common_cells, "cell_type"]
# 添加到 proj_combined 中
proj_combined <- addCellColData(ArchRProj = proj_combined, data = cell_type_full, name = "celltype", cells = names(cell_type_full), force = TRUE)
keep_cells <- which(!is.na(proj_combined$celltype))
# 取出这些细胞名
keep_cell_names <- proj_combined$cellNames[keep_cells]
# 子集化 ArchR 对象
proj_combined_filtered <- subsetCells(ArchRProj = proj_combined, cellNames = keep_cell_names)
proj_combined_filtered <- addUMAP(ArchRProj = proj_combined_filtered, reducedDims = "Harmony", force = TRUE)
pdf("root-silique-celltype.pdf",width = 10,height = 6)
plotEmbedding(ArchRProj = proj_combined_filtered, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
dev.off()
pdf("root-silique-organ.pdf",width = 10,height = 6)
plotEmbedding(ArchRProj = proj_combined_filtered, colorBy = "cellColData", name = "organ", embedding = "UMAP")
dev.off()
saveRDS(proj_combined_filtered, "PCA_snATAC_root_silique_celltype.rds")
#=========================================================================================####
