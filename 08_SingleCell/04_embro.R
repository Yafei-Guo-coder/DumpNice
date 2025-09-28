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
file <- list.files(path = ".", pattern = "\\.arrow$", full.names = TRUE)[-8]
proj_merged <- ArchRProject(
  geneAnnotation = gene_annotation,
  genomeAnnotation = genome_annotation,
  ArrowFiles = file,
  outputDirectory = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/Merged_proj",
  copyArrows = TRUE
)

saveArchRProject(
  ArchRProj = proj_merged,
  outputDirectory = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/Merged_proj",
  load = FALSE
)

#====================================================================================================================#
###                                        create archR project                                                   ####
#====================================================================================================================#
projAll <- loadArchRProject(path = "./Merged_proj")

globular <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_globular_anchor_refine.rds")
# 23146 cells
heart <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_heart_anchor_refine.rds")
# 15877 cells
torpedo <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_torpedo_anchor_refine.rds")
# 11819 cells
bent <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_bent_anchor_refine.rds")
# 21953 cells
cotyledon <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/rds/ArchR_cotyledon_anchor_refine.rds")
# 38038 cells
projAll <- projAll[c(globular$cellNames, heart$cellNames, torpedo$cellNames, bent$cellNames, cotyledon$cellNames),]

# 合并所有 predictedGroup 信息
cell_meta <- data.frame(
  cellNames = c(globular$cellNames,
                heart$cellNames,
                torpedo$cellNames,
                bent$cellNames,
                cotyledon$cellNames),
  predictedGroup = c(globular$predictedGroup,
                     heart$predictedGroup,
                     torpedo$predictedGroup,
                     bent$predictedGroup,
                     cotyledon$predictedGroup),
  stage = c(rep("globular", length(globular$cellNames)),
            rep("heart", length(heart$cellNames)),
            rep("torpedo", length(torpedo$cellNames)),
            rep("bent", length(bent$cellNames)),
            rep("cotyledon", length(cotyledon$cellNames))),
  stringsAsFactors = FALSE
)
cell_meta$stage_predictedGroup <- paste(cell_meta$stage, cell_meta$predictedGroup, sep = "_")
# 添加到 projAll
# 添加到 projAll
projAll <- addCellColData(
  ArchRProj = projAll,
  data = cell_meta$predictedGroup,
  cells = cell_meta$cellNames,
  name = "predictedGroup",
  force = TRUE
)

projAll <- addCellColData(
  ArchRProj = projAll,
  data = cell_meta$stage,
  cells = cell_meta$cellNames,
  name = "stage",
  force = TRUE
)

projAll <- addCellColData(
  ArchRProj = projAll,
  data = cell_meta$stage_predictedGroup,
  cells = cell_meta$cellNames,
  name = "stage_predictedGroup",
  force = TRUE
)
saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_anchor_all.rds")

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

# 添加motif
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
projAll <- addMotifAnnotations(ArchRProj = projAll, motifPWMs = pwmlist, name = "Motif")

# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)

saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_anchor_all_predictedGroupMotif.rds")


# 上面添加了peakmatrix和motifmatrix，接下来鉴定差异的peak和motif，这里是鉴定每个分类（前景）和其他所有分类（背景）的差异，也可以指定两个分类变量做差异
# 热图
markerPeaks <- getMarkerFeatures(
  ArchRProj = projAll, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
heatmapPeaks <- markerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
pdf("Peak-Marker-Heatmap.pdf", width = 12, height = 5)
print(heatmapPeaks)
dev.off()
# plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projAll, addDOC = FALSE)
# 火山图
# pma <- plotMarkers(seMarker = markerPeaks, name = "General embryo", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pv <- plotMarkers(seMarker = markerPeaks, name = "General embryo", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pdf("General-embryo-Markers-Volcano.pdf")
print(pv)
dev.off()

# plotPDF(pv, name = "General-embryo-Markers-Volcano", width = 5, height = 5, ArchRProj = projAll, addDOC = FALSE)
# 在基因组上的track
TF <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/txt/ruiying_TF.txt", header=F, stringsAsFactors=F)[,1]
p <- plotBrowserTrack(
  ArchRProj = projAll, 
  groupBy = "predictedGroup", 
  geneSymbol = c("AT3G26744"),
  features =  getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["General embryo"],
  upstream = 5000,
  downstream = 3000
)
pdf("Plot-Tracks-With-Features.pdf",width = 12, height = 8)
print(p[[1]])
dev.off()
# plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projAll, addDOC = FALSE)

# 上调的motif，这里可以看marker peaks，就是每种分类特异性的peaks，也可以看差异peaks，比如说两种分类的差异peak，指定一组peak，就可以看相对的motif富集的上调和下调
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = projAll,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5" # "FDR <= 0.1 & Log2FC <= -0.5"下调
)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = T)
# plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projAll, addDOC = FALSE)
pdf("Motifs-Enriched-Marker-Heatmap.pdf", width = 12, height = 5)
print(heatmapEM)
dev.off()

# 画motif的logo
pwm <- getPeakAnnotation(projAll, "Motif")$motifs[["AT5G47670"]]
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}
ppm <- PWMatrixToProbMatrix(pwm)
colSums(ppm) %>% range
library(ggseqlogo)
ggseqlogo(ppm, method = "bits")
ggseqlogo(ppm, method = "prob")

# 保存
saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/ED_0911/ArchR_EndospermD_UMAP_rmTorpedo2_motif.rds")

#====================================================================================================================#
####                                              trajectory                                                       ###
#====================================================================================================================#
# anchor ##
projAll <- addImputeWeights(projAll)
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
cell_ids <- seRNA@meta.data$cellID
cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
seRNA_filtered <- subset(seRNA, cells = cells_to_keep)
projAll <- addGeneIntegrationMatrix(ArchRProj = projAll, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA_filtered, addToArrow = TRUE, force = TRUE, groupRNA = "endotype2025", nameCell = "predictedCell_endo2025", nameGroup = "predictedGroup_endo2025", nameScore = "predictedScore_endo2025", threads = 1)
saveRDS(projAll, "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/ED_0911/ArchR_EndospermD_UMAP_rmTorpedo2_motif_reanchor2.rds")

#====================================================================================================================#
####                                         cell type corelation heatmap                                          ###
#====================================================================================================================#
seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/EndoCells_recluster_pseudotime.rds")
cell_ids <- seRNA@meta.data$cellID
cells_to_keep <- cell_ids[!grepl("^(2_28h|3_48h)", cell_ids)]
seRNA_filtered <- subset(seRNA, cells = cells_to_keep)

projAll <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/rds/ArchR_EndospermD_rmTorpedo2_reAnchor2_trajectory.rds")