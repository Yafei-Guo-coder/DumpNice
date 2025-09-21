# 新测序文件
# /jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L01/V350319993_L01_read_1.fq.gz
# /jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L01/V350319993_L01_read_2.fq.gz
# /jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L02/V350319993_L02_read_1.fq.gz
# /jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L02/V350319993_L02_read_2.fq.gz
# /jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L03/V350319993_L03_read_1.fq.gz
# /jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L03/V350319993_L03_read_2.fq.gz
# /jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L04/V350319993_L04_read_1.fq.gz
# /jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L04/V350319993_L04_read_2.fq.gz

# run.sh
/usr/bin/java -jar /jdfsbjcas1/ST_BJ/PUB/Tool/bin/cromwell-35.jar run /jdfsbjcas1/ST_BJ/PUB/Pipeline/DNBelabC4scATAC/V3/scATAC_nonmodel_local_V3.0_BJ.wdl --inputs iDrop_scATAC_default.json

(base) [guoyafei1@bjcas-compute-2-4 V3]$ cat iDrop_scATAC_default.json
{
"iDrop_scATAC_default_2.Data" : [ [ "/jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L01/V350319993_L01_read_1.fq.gz", "/jdfsbjcas1/BJ_DATA/MGISEQ2000/R10040100200001/V350319993/L01/V350319993_L01_read_2.fq.gz" ] ],
"iDrop_scATAC_default_2.Outdir" : "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/02_fq/L01/",
"iDrop_scATAC_default_2.SampleID" : "L01",
"iDrop_scATAC_default_2.reference" : "Araport11",
"iDrop_scATAC_default_2.ProjectID" : "P21Z28400N0234",
"iDrop_scATAC_default_2.readStructure" : "/jdfsbjcas1/ST_BJ/PUB/Pipeline/DNBelabC4scATAC/script/C4scATAClib_seqT1_R1_70_R2_50.json"
}


# 在集群上跑出报告之后，使用ArchR进行过滤和聚类
# mv /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/02_fq/L01/d2cfile/L01.fragments.tsv.gz* /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/02_ATACseed_test
# mv /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/02_fq/L02/d2cfile/L02.fragments.tsv.gz* /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/02_ATACseed_test
# mv /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/02_fq/L03/d2cfile/L03.fragments.tsv.gz* /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/02_ATACseed_test
# mv /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/02_fq/L04/d2cfile/L04.fragments.tsv.gz* /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/02_ATACseed_test

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
#====================================================================================================================#
#                                        创建基因组注释及ArchR对象                                                ####
#====================================================================================================================#
# 从BSgenome对象中创建基因组注释
genome_annotation <- createGenomeAnnotation(genome = BSgenome.Athaliana.TAIR10.02022009)
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

# 获取所有子目录名称
subdirs <- list.dirs("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/02_ATACseed_test/fragments", recursive = F)
# 初始化空向量用于存储所有样本名称和片段文件路径
all_samples <- character()
all_fragments <- character()
# 遍历每个子目录
for (dir in subdirs) {
  # 列出当前子目录下所有以 ".tsv.gz" 结尾的片段文件（不递归）
  fragments <- list.files(path = dir, pattern = ".tsv.gz$", recursive = T, full.names = F)
  # 提取每个片段文件的样本名称（假设样本名称是文件名的第一部分）
  samples <- unlist(lapply(fragments, function(x){
    unlist(strsplit(x, ".", fixed = T))[1]
  }))
  # 将样本名称和片段文件路径添加到总向量中
  all_samples <- c(all_samples, samples)
  all_fragments <- c(all_fragments, paste0(dir, "/", fragments))
}
cat("All samples:\n")
for (sample in all_samples) {cat(sample, "\n")}
cat("\nAll fragments:\n")
for (fragment in all_fragments) {cat(fragment, "\n")}
# 设置 ArchR 使用的线程数为 8
addArchRThreads(threads = 5)
# 创建 Arrow 文件
arrowFiles <- createArrowFiles(inputFiles = all_fragments,       # 输入的片段文件
                               sampleNames = all_samples,        # 样本名称
                               outputNames = all_samples,        # 输出文件名称
                               minTSS = 0,                   # 最小 TSS 值
                               minFrags = 0,                 # 最小片段数
                               excludeChr = c("ChrC", "ChrM"), # 排除的染色体
                               addTileMat = T,               # 是否添加 Tile 矩阵
                               addGeneScoreMat = T,          # 是否添加基因评分矩阵
                               subThreading = F,             # 是否使用子线程 # 当设置线程数大于样本数时会启用这个参数，以确认是否使用多线程并行处理同一个样本，而多线程并行处理同一个样本对系统要求较高
                               geneAnnotation = gene_annotation, # 基因注释
                               genomeAnnotation = genome_annotation) # 基因组注释
# 只需要创建一次
# 列出当前目录下所有以 ".arrow" 结尾的 Arrow 文件（不递归）
arrowFiles <- list.files(pattern = ".arrow", recursive = F, full.names = F)
# 创建一个ArchRProject对象
proj_all <- ArchRProject(
  ArrowFiles = arrowFiles,          # 输入的Arrow文件列表
  geneAnnotation = gene_annotation, # 基因注释对象
  genomeAnnotation = genome_annotation, # 基因组注释对象
  outputDirectory = "seed",       # 输出目录，用于存储项目的输出文件
  showLogo = F,                     # 是否显示ArchR的logo，设置为FALSE表示不显示
  copyArrows = F,                    # 是否复制Arrow文件，设置为FALSE表示不复制 # 推荐用T，如果你修改了Arrow文件，可以保留一个原始副本以备后用
  threads = 8
)

saveRDS(proj_all, "snATAC_seed_test.rds")

#######################################################################

projHeme1 <- readRDS("snATAC_seed_test.rds")
# r <- sample(projHeme1$cellNames, 200000) 
# projHeme1 <- projHeme1[which(projHeme1$cellNames %in% r), ]
# 原始数据可视化

df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment","Sample"))
heart1 <- df[which(df$Sample == names(table(projHeme1$Sample)[1])),]
torpedo1 <- df[which(df$Sample == names(table(projHeme1$Sample)[2])),]
torpedo3 <- df[which(df$Sample == names(table(projHeme1$Sample)[3])),]
cotyledon2 <- df[which(df$Sample == names(table(projHeme1$Sample)[4])),]

heart1_plot <- ggPoint(x = heart1[,1], y = heart1[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") 
torpedo1_plot <- ggPoint(x = torpedo1[,1], y = torpedo1[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") 
torpedo3_plot <- ggPoint(x = torpedo3[,1], y = torpedo3[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") 
cotyledon2_plot <- ggPoint(x = cotyledon2[,1], y = cotyledon2[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") 

p2 <- plotGroups(ArchRProj = projHeme1, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p4 <- plotGroups(ArchRProj = projHeme1, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
df2 <- as.data.frame(df)
colnames(df2)[1] <-"log10nFrags"

# TSS富集
p7 <- plotFragmentSizes(ArchRProj = projHeme1)
p8 <- plotTSSEnrichment(ArchRProj = projHeme1)

pdf("raw_test_plot.pdf",width = 6,height = 8)
print(heart1_plot)
print(torpedo1_plot)
print(torpedo3_plot)
print(cotyledon2_plot)

print(p2)
print(p4)
print(p7)
print(p8)
dev.off()

# 过滤以及可视化
projHeme1 <- readRDS("snATAC_seed_test.rds")
idxPass <- which(projHeme1$TSSEnrichment >= 1 & projHeme1$nFrags>= 1500)
cellsPass <- projHeme1$cellNames[idxPass]
projHeme2 <- projHeme1[cellsPass, ]
projHeme2 <- addDoubletScores(input = projHeme2, k = 10, knnMethod = "UMAP", LSIMethod = 1)
projHeme2 <- filterDoublets(projHeme2)
saveRDS(projHeme2, "PCA_snATAC_filter_TSS_doublets_test.rds")
projHeme1 <- projHeme2
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment","Sample"))
heart1 <- df[which(df$Sample == names(table(projHeme1$Sample)[1])),]
torpedo1 <- df[which(df$Sample == names(table(projHeme1$Sample)[2])),]
torpedo3 <- df[which(df$Sample == names(table(projHeme1$Sample)[3])),]
cotyledon2 <- df[which(df$Sample == names(table(projHeme1$Sample)[4])),]

heart1_plot <- ggPoint(x = heart1[,1], y = heart1[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") 
torpedo1_plot <- ggPoint(x = torpedo1[,1], y = torpedo1[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") 
torpedo3_plot <- ggPoint(x = torpedo3[,1], y = torpedo3[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") 
cotyledon2_plot <- ggPoint(x = cotyledon2[,1], y = cotyledon2[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") 

p2 <- plotGroups(ArchRProj = projHeme1, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p4 <- plotGroups(ArchRProj = projHeme1, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
df2 <- as.data.frame(df)
colnames(df2)[1] <-"log10nFrags"

# TSS富集
p7 <- plotFragmentSizes(ArchRProj = projHeme1)
p8 <- plotTSSEnrichment(ArchRProj = projHeme1)

pdf("filter_test_plot.pdf",width = 6,height = 8)
print(heart1_plot)
print(torpedo1_plot)
print(torpedo3_plot)
print(cotyledon2_plot)

print(p2)
print(p4)
print(p7)
print(p8)
dev.off()

> table(projHeme1$Sample)
  L01   L02   L03   L04 
25607  1394 14947  8509 

> table(projHeme2$Sample)
  L01   L02   L03   L04 
10455   519  8784  6119 

##################################################### 聚类 ###########################################

library(pheatmap)
for (j in c(1.5,1.0,0.5)) {
      projHeme2 <- readRDS("PCA_snATAC_filter_TSS_doublets_test.rds")
      projHeme2 <- addIterativeLSI(ArchRProj = projHeme2, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 4, clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10), varFeatures = 15000, dimsToUse = 1:30)
      projHeme2 <- addClusters(input = projHeme2, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = j)
      projHeme2 <- addUMAP(ArchRProj = projHeme2, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine")
      projHeme2 <- addHarmony(ArchRProj = projHeme2, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample")
      projHeme2 <- addUMAP(ArchRProj = projHeme2, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, metric = "cosine")
      p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
      p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
      p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
      p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
      outfile2 <- paste("cluster_resolut_",j,".pdf",sep="")
      pdf(outfile2,width = 12,height = 12)
        grid.arrange(p1,p2,p3,p4,nrow=2)
      dev.off()
}



