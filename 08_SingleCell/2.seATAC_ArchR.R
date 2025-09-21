#====================================================================================================================#
#                                                  加载包                                                         ####
#====================================================================================================================#
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
library(Seurat)
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
library(SeuratObject)
.libPaths("/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
# devtools::install_version("matrixStats", version="1.5.0", lib = "./lib")
# install.packages("universalmotif", lib = "./lib")
library(universalmotif)
library(tidydr)
# options(mc.cores = 1) 
# 参考网站 https://www.jianshu.com/p/1c3852f123e7 (简书：使用ArchR分析单细胞ATAC-seq数据(第二章))
# 一些基本操作,ArchR并不会直接提取数据，而是构建一个新的ArchRProject对象
# 根据位置
# projHeme1[1:100, ]
# 根据细胞名
# projHeme1[projHeme1$cellNames[1:100], ]
# sample name
# idxSample <- BiocGenerics::which(projHeme1$Sample %in% "scATAC_BMMC_R1")
# cellsSample <- projHeme1$cellNames[idxSample]
# projHeme1[cellsSample, ]
# TSS enrichment score
# idxPass <- which(projHeme1$TSSEnrichment >= 8)
# cellsPass <- projHeme1$cellNames[idxPass]
# projHeme1[cellsPass, ]
# 另一种添加方式
# bioNames <- bioNames[1:10]
# cellNames <- projHeme1$cellNames[1:10]
# projHeme1 <- addCellColData(ArchRProj = projHeme1, data = paste0(bioNames),cells = cellNames, name = "bioNames2")
# getCellColData(projHeme1, select = c("bioNames", "bioNames2"))
# df <- getCellColData(projHeme1, select = "nFrags")
# df <- getCellColData(projHeme1, select = c("log10(nFrags)", "nFrags - 1"))
# 可以通过下面这行代码了解下ArchRProject在R中会使用多少内存
# paste0("Memory Size = ", round(object.size(projHeme1) / 10^6, 3), " MB")
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
#gene_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.gene.txt"
#trans_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.transcript.txt"
#exon_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/00.ref/Arabidopsis/Araport11/snATAC_ref.ArchR/Araport11.Mar92021.exon.txt"

gene_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/00_ref/Arabidopsis/Araport11/Araport11.Mar92021.gene.noCM.txt"
trans_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/00_ref/Arabidopsis/Araport11/Araport11.Mar92021.transcript.noCM.txt"
exon_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/00_ref/Arabidopsis/Araport11/Araport11.Mar92021.exon.noCM.txt"


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
subdirs <- list.dirs("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment", recursive = F)
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
all_fragments <- all_fragments[2:3]
cat("All samples:\n")
for (sample in all_samples) {cat(sample, "\n")}
cat("\nAll fragments:\n")
for (fragment in all_fragments) {cat(fragment, "\n")}
# 设置 ArchR 使用的线程数为 8
addArchRThreads(threads = 8)
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
projAll <- ArchRProject(
  ArrowFiles = arrowFiles,          # 输入的Arrow文件列表
  geneAnnotation = gene_annotation, # 基因注释对象
  genomeAnnotation = genome_annotation, # 基因组注释对象
  outputDirectory = "Development",       # 输出目录，用于存储项目的输出文件
  showLogo = F,                     # 是否显示ArchR的logo，设置为FALSE表示不显示
  copyArrows = F,                    # 是否复制Arrow文件，设置为FALSE表示不复制 # 推荐用T，如果你修改了Arrow文件，可以保留一个原始副本以备后用
  threads = 8
)
saveRDS(projAll, "snATAC_torpedo1_bent1.rds")
#====================================================================================================================#
#                                        质控以及基本数据统计可视化                                               ####
#====================================================================================================================#
projHeme1 <- readRDS("PCA_snATAC_all.rds")
idxPass <- which(projHeme1$TSSEnrichment >= 1.5 & projHeme1$nFrags>= 1500)
cellsPass <- projHeme1$cellNames[idxPass]
projHeme2 <- projHeme1[cellsPass, ]
doubScores <- addDoubletScores(input = projHeme2, k = 10, knnMethod = "UMAP", LSIMethod = 1)
saveRDS(doubScores, "PCA_snATAC_filter_TSS.rds")
projHeme3 <- filterDoublets(doubScores)
#saveRDS(projHeme3, "PCA_snATAC_filter_TSS_doublets.rds")
projHeme2 <- readRDS("PCA_snATAC_filter_TSS_doublets.rds")
projHeme2$bioNames <- gsub("_AT","",projHeme2$Sample)

df <- getCellColData(projHeme2, select = c("log10(nFrags)", "TSSEnrichment","PromoterRatio","bioNames"))
silique <- df[which(df$bioNames %in% c("A220221001_Silique_0221","A220322001_Silique_0322","A220322002_Silique_0322"),),]
silique1 <- df[which(df$bioNames == "A220221001_Silique_0221"),]
silique2 <- df[which(df$bioNames == "A220322001_Silique_0322"),]
silique3 <- df[which(df$bioNames == "A220322002_Silique_0322"),]
root <- df[which(df$bioNames %in% c("A220316001_7d-Root_0316","A220316002_7d-Root_0316","A220322003_13d-Root_0322","A220322004_13d-Root_0322","A220711005_7d-Root_0711","A220711006_7d-Root_0711"),),]
root1 <- df[which(df$bioNames == "A220316001_7d-Root_0316"),]
root2 <- df[which(df$bioNames == "A220316002_7d-Root_0316"),]
root3 <- df[which(df$bioNames == "A220322003_13d-Root_0322"),]
root4 <- df[which(df$bioNames == "A220322004_13d-Root_0322"),]
root5 <- df[which(df$bioNames == "A220711005_7d-Root_0711"),]
root6 <- df[which(df$bioNames == "A220711006_7d-Root_0711"),]
shoot <- df[which(df$bioNames %in% c("A221008001_Shoot_1008","A221008002_Shoot_1008"),),]
shoot1 <- df[which(df$bioNames == "A221008001_Shoot_1008"),]
shoot2 <- df[which(df$bioNames == "A221008002_Shoot_1008"),]
# 原始细胞数
# silique:35219; shoot:18274; root:55146; silique1:13884; silique2:9888; silique3:11447; root1:6564; root2:6570; root3:8314; root4:10077; root5:11859; root6:11762; shoot1:8289; shoot2:9985
# 画密度图
silique1_plot <- ggPoint(x = silique1[,1], y = silique1[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
silique2_plot <- ggPoint(x = silique2[,1], y = silique2[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
silique3_plot <- ggPoint(x = silique3[,1], y = silique3[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
silique_plot <- ggPoint(x = silique[,1], y = silique[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
root1_plot <- ggPoint(x = root1[,1], y = root1[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
root2_plot <- ggPoint(x = root2[,1], y = root2[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
root3_plot <- ggPoint(x = root3[,1], y = root3[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
root4_plot <- ggPoint(x = root4[,1], y = root4[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
root5_plot <- ggPoint(x = root5[,1], y = root5[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
root6_plot <- ggPoint(x = root6[,1], y = root6[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
root_plot <- ggPoint(x = root[,1], y = root[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
shoot1_plot <- ggPoint(x = shoot1[,1], y = shoot1[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
shoot2_plot <- ggPoint(x = shoot2[,1], y = shoot2[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
shoot_plot <- ggPoint(x = shoot[,1], y = shoot[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500)+0.5, quantile(df[,1], probs = 0.999)), ylim = c(0, quantile(df[,2], probs = 0.999))) 
dev.off()
pdf("silique_raw_plot1.pdf",width = 21,height = 7)
grid.arrange(silique1_plot, silique2_plot, silique3_plot, nrow=1)
dev.off()
pdf("root_raw_plot1.pdf",width = 21,height = 15)
grid.arrange(root1_plot, root2_plot, root3_plot, root4_plot,root5_plot,root6_plot,nrow=2)
dev.off()
pdf("shoot_raw_plot1.pdf",width = 15,height = 7)
grid.arrange(shoot1_plot, shoot2_plot,nrow=1)
dev.off()
pdf("all_raw_plot1.pdf",width = 15,height = 7)
grid.arrange(silique_plot,root_plot,shoot_plot,nrow=1)
dev.off()
# 画山脊图
p1 <- plotGroups(ArchRProj = projHeme2, groupBy = "Sample", colorBy = "cellColData", orderGroups = c("A220316001_7d-Root_AT_0316","A220316002_7d-Root_AT_0316","A220711005_7d-Root_AT_0711","A220711006_7d-Root_AT_0711","A220322003_13d-Root_AT_0322", "A220322004_13d-Root_AT_0322", "A221008001_Shoot_AT_1008","A221008002_Shoot_AT_1008","A220221001_Silique_AT_0221","A220322001_Silique_AT_0322","A220322002_Silique_AT_0322"), name = "TSSEnrichment", plotAs = "ridges")
p2 <- plotGroups(ArchRProj = projHeme2, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", orderGroups = c("A220316001_7d-Root_AT_0316","A220316002_7d-Root_AT_0316","A220711005_7d-Root_AT_0711","A220711006_7d-Root_AT_0711","A220322003_13d-Root_AT_0322", "A220322004_13d-Root_AT_0322", "A221008001_Shoot_AT_1008","A221008002_Shoot_AT_1008","A220221001_Silique_AT_0221","A220322001_Silique_AT_0322","A220322002_Silique_AT_0322"), plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p3 <- plotGroups(ArchRProj = projHeme2, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "ridges")
p4 <- plotGroups(ArchRProj = projHeme2, groupBy = "Sample", colorBy = "cellColData", orderGroups = c("A220316001_7d-Root_AT_0316","A220316002_7d-Root_AT_0316","A220711005_7d-Root_AT_0711","A220711006_7d-Root_AT_0711","A220322003_13d-Root_AT_0322", "A220322004_13d-Root_AT_0322", "A221008001_Shoot_AT_1008","A221008002_Shoot_AT_1008","A220221001_Silique_AT_0221","A220322001_Silique_AT_0322","A220322002_Silique_AT_0322"), name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
dev.off()
pdf("filter_plot1.pdf",width = 8,height = 5)
grid.arrange(p1,p3,nrow=1)
dev.off()
dev.off()
pdf("filter_plot2.pdf",width = 9,height = 4.5)
grid.arrange(p2,p4,nrow=1)
dev.off()
# 画箱线图 
df2 <- as.data.frame(df)
colnames(df2)[1] <-"log10nFrags"
df2$bioNames <- factor(df2$bioNames,levels = c("A221008002_Shoot_1008","A221008001_Shoot_1008","A220322002_Silique_0322","A220322001_Silique_0322","A220221001_Silique_0221","A220711006_7d-Root_0711","A220711005_7d-Root_0711","A220322004_13d-Root_0322","A220322003_13d-Root_0322","A220316002_7d-Root_0316","A220316001_7d-Root_0316"))
p1 <- ggplot(df2, aes(x=log10nFrags, y=bioNames, fill=bioNames)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") + theme_classic() + theme(legend.position="none", plot.title = element_text(size=11)) + xlab("")
p2 <- ggplot(df2, aes(x=TSSEnrichment, y=bioNames, fill=bioNames)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") + theme_classic() + theme(legend.position="none", plot.title = element_text(size=11)) + xlab("")
p3 <- ggplot(df2, aes( x=PromoterRatio, y=bioNames, fill=bioNames)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") + theme_classic() + theme(legend.position="none", plot.title = element_text(size=11)) + xlab("")
pdf("all_filter_boxplot4.pdf", width = 20, height = 7)
grid.arrange(p1,p2,p3,nrow=1)
dev.off()
# 画片段长度和TSS富集图 
p5 <- plotFragmentSizes(ArchRProj = projHeme2)
p6 <- plotTSSEnrichment(ArchRProj = projHeme2)
dev.off()
pdf("filter_plot3.pdf",width = 6,height = 8)
grid.arrange(p5, p6, nrow=2)
dev.off()
# 数据保存
saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme-filter", load = FALSE)

# do QC plot
cotyledon_atac <- readRDS("rds/ArchR_bent_filterd.rds")
bent_atac <- readRDS("rds/ArchR_cotyledon_filterd.rds")
globular_atac <- readRDS("rds/ArchR_globular_filterd.rds")
heart_atac <- readRDS("rds/ArchR_heart_filterd.rds")
torpedo_atac <- readRDS("rds/ArchR_torpedo_filterd.rds")

globular <- getCellColData(globular_atac, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
globular$type <- "globular"
heart <- getCellColData(heart_atac, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
heart$type <- "heart"
torpedo <- getCellColData(torpedo_atac, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
torpedo$type <- "torpedo"
bent <- getCellColData(bent_atac, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
bent$type <- "bent"
cotyledon <- getCellColData(cotyledon_atac, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
cotyledon$type <- "cotyledon"


all <- as.data.frame(rbind(globular,heart,torpedo,bent,cotyledon))
names(all) <- c("log10nFrags","TSSEnrichment","Sample","type")
all2 <- all[which(all$Sample != "A250410002_torpedo2"),]

write.table(all2,"QC_plot.txt",sep="\t",row.names= F, quote=F)
data$Sample <- factor(data$Sample, levels = c("A250401001_globular1", "A250401002_globular2", "A250401003_globular3", "A250408001_heart1", "A250408002_heart2", "A250408003_heart3", "A250410001_torpedo1","A250410003_torpedo3",  "A250415004_bent1", "A250415005_bent2", "A250415006_bent3", "A250415007_bent4", "A250415001_cotyledon1", "A250415002_cotyledon2", "A250415003_cotyledon3"))



ggplot(data, aes(x = Sample, y = TSSEnrichment, fill = type)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_classic() +
  scale_fill_manual(values = rev(viridis::viridis(9, alpha = 0.6)[2:6])) +
  labs(x = "Batch", y = "", fill = "Organ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=1.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(limits = c(1, 5))

"#44015499" "#472D7B99" "#3B528B99" "#2C728E99" "#21908C99" "#27AD8199" "#5DC86399" "#AADC3299" "#FDE72599"
#====================================================================================================================#
#                                           细胞聚类及批次效应消除                                                ####
#====================================================================================================================#
#去掉两个文库，"A220711005_7d-Root_AT_0711","A220711006_7d-Root_AT_0711"，这俩数据跟其他的文库合不上
root <- c("A220316001_7d-Root_AT_0316", "A220316002_7d-Root_AT_0316","A220322003_13d-Root_AT_0322","A220322004_13d-Root_AT_0322")
shoot <- c("A221008001_Shoot_AT_1008","A221008002_Shoot_AT_1008") 
silique <- c("A220221001_Silique_AT_0221","A220322001_Silique_AT_0322","A220322002_Silique_AT_0322")
#root_7d <- c("A220316001_7d-Root_AT_0316", "A220316002_7d-Root_AT_0316","A220711005_7d-Root_AT_0711","A220711006_7d-Root_AT_0711")
#root_13d <- c("A220322003_13d-Root_AT_0322","A220322004_13d-Root_AT_0322")
#tissue <- list(root_a,root_b,root_c,root_d,root_e,root_f,root_7d,root_13d)
#names(tissue) <- c("root_a","root_b","root_c","root_d","root_e","root_f","root_7d","root_13d")
tissue <- list(root, shoot, silique)
names(tissue) <- c("root_4lib","shoot","silique")
library(pheatmap)
#for (j in c(0.8, 0.6, 0.4)) {
for (j in c(1.8,1.2,1.0,0.8)) {  
  i=1
  projHeme2 <- loadArchRProject(path = "Save-ProjHeme-filter")
  idxSample <- BiocGenerics::which(projHeme2$Sample %in% tissue[[i]])
  #idxSample <- BiocGenerics::which(projHeme2$Sample %in% root_7d)
  cellsSample <- projHeme2$cellNames[idxSample]
  projHeme <- projHeme2[cellsSample, ]
  projHeme2 <- addIterativeLSI(ArchRProj = projHeme, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 4, clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10), varFeatures = 15000, dimsToUse = 1:30)
  projHeme2 <- addClusters(input = projHeme2, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = j)
  cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
  cM <- cM / Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), border_color = "black")
  projHeme2 <- addUMAP(ArchRProj = projHeme2, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine")
  projHeme2 <- addHarmony(ArchRProj = projHeme2, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample")
  projHeme2 <- addUMAP(ArchRProj = projHeme2, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, metric = "cosine")
  outfile <- paste("PCA_snATAC_cluster_",names(tissue)[i],"_resolut_",j,".rds",sep="")
  saveRDS(projHeme2, outfile)
  dev.off()
  dev.off()
  dev.off()
  outfile1 <- paste("filter_seurat_heatmap_",names(tissue)[i],"_resolut_",j,".pdf",sep="")
  pdf(outfile1,width = 8,height = 8)
  print(p)
  dev.off()
  #测试
  #projHeme1 <- readRDS("PCA_snATAC_cluster_root_resolut_0.8.rds")
  cellsSample <- projHeme1$cellNames[BiocGenerics::which(projHeme1$Sample %in% root_a)]
  projHeme_7d_a <- projHeme1[cellsSample, ]
  cellsSample <- projHeme1$cellNames[BiocGenerics::which(projHeme1$Sample %in% root_b)]
  projHeme_7d_b <- projHeme1[cellsSample, ]
  cellsSample <- projHeme1$cellNames[BiocGenerics::which(projHeme1$Sample %in% root_c)]
  projHeme_7d_c <- projHeme1[cellsSample, ]
  cellsSample <- projHeme1$cellNames[BiocGenerics::which(projHeme1$Sample %in% root_d)]
  projHeme_7d_d <- projHeme1[cellsSample, ]
  cellsSample <- projHeme1$cellNames[BiocGenerics::which(projHeme1$Sample %in% root_e)]
  projHeme_7d_e <- projHeme1[cellsSample, ]
  p1 <- plotEmbedding(ArchRProj = projHeme_7d_a, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p2 <- plotEmbedding(ArchRProj = projHeme_7d_a, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
  p3 <- plotEmbedding(ArchRProj = projHeme_7d_b, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p4 <- plotEmbedding(ArchRProj = projHeme_7d_b, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
  p5 <- plotEmbedding(ArchRProj = projHeme_7d_c, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p6 <- plotEmbedding(ArchRProj = projHeme_7d_c, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
  p7 <- plotEmbedding(ArchRProj = projHeme_7d_d, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p8 <- plotEmbedding(ArchRProj = projHeme_7d_d, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
  p9 <- plotEmbedding(ArchRProj = projHeme_7d_e, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p10 <- plotEmbedding(ArchRProj = projHeme_7d_e, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
  pdf("sample_root_combination.pdf",width = 20,height = 14)
  grid.arrange(p1,p3,p5,p7,p9,nrow=2)
  dev.off()
  pdf("cluster_root_combination.pdf",width = 20,height = 14)
  grid.arrange(p2,p4,p6,p8,p10, nrow=2)
  dev.off()
  p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
  p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
  p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
  dev.off()
  outfile2 <- paste("filter_umap_cluster_",names(tissue)[i],"_resolut_",j,".pdf",sep="")
  pdf(outfile2,width = 10,height = 6)
  grid.arrange(p1,p2,nrow=1)
  dev.off()
  outfile3 <- paste("filter_UMAPHarmony_cluster_",names(tissue)[i],"_resolut_",j,".pdf",sep="")
  pdf(outfile3,width = 10,height = 6)
  grid.arrange(p3,p4,nrow=1)
  dev.off()
}
root <- readRDS("PCA_snATAC_cluster_root_resolut_0.8.rds")
cluster_labels <- getCellColData(projHeme2, "Clusters")
cluster_labels <- as.numeric(as.factor(cluster_labels))
embedding_matrix <- getReducedDims(projHeme2, reducedDims = "IterativeLSI")  # 或 "PCA", "UMAP"，"UMAPHarmony"
distance_matrix <- dist(embedding_matrix)
sil <- silhouette(cluster_labels, distance_matrix)
# 画出轮廓系数
plot(sil, col = 1:max(cluster_labels), border = NA)
#====================================================================================================================#
#                                            鉴定并测试marker基因                                                 ####
#====================================================================================================================#
# silique
# 计算基因得分，并且鉴定标记基因
projHeme2 <- readRDS("PCA_snATAC_cluster_silique_resolut_0.8.rds")
markersGS <- getMarkerFeatures(ArchRProj = projHeme1, useMatrix = "GeneScoreMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
# 参考陈迪俊的文章的阈值parameter “cutOff = FDR < = 0.01 & Log2FC > = 1”
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# 写出marker基因的列表
for (i in seq_along(markerList)) {
  df <- as.data.frame(markerList[[i]])
  write.table(df, paste0("markerList_silique_resolut0.8_", i, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
# 画热图
markerGenes  <- c("AT1G05220","AT1G03103","AT1G23910","AT1G04357")
heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", labelMarkers = markerGenes,transpose = TRUE)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)
dev.off()
# 画散点图
projHeme2 <- addImputeWeights(projHeme2)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAPHarmony", imputeWeights = getImputeWeights(projHeme2))
pdf("Plot-UMAPHarmony-Marker-Genes-W-Imputation.pdf",width = 35,height = 35)
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[9]],p1[[10]],p1[[11]],p1[[12]],p1[[13]],p1[[14]],p1[[15]],p1[[16]],p1[[17]],p1[[18]],p1[[19]],p1[[20]],p1[[21]],p1[[22]],p1[[23]],p1[[24]],p1[[25]],nrow=5)
dev.off()
# 画基因组通道
p3 <- plotBrowserTrack(ArchRProj = projHeme2, groupBy = "Clusters", geneSymbol = markerGenes, upstream = 30000, downstream = 30000)
pdf("Plot-Tracks-Marker-Genes.pdf",width = 55,height = 55)
grid.arrange(p3[[1]],p3[[2]],p3[[3]],p3[[4]],p3[[5]],p3[[6]],p3[[7]],p3[[8]],p3[[9]],p3[[10]],p3[[11]],p3[[12]],p3[[13]],p3[[14]],p3[[15]],p3[[16]],p3[[17]],p3[[18]],p3[[19]],p3[[20]],p3[[21]],p3[[22]],p3[[23]],p3[[24]],p3[[25]],nrow=5)
dev.off()

# root
projHeme2 <- readRDS("PCA_snATAC_cluster_root_4lib_resolut_0.8.rds")
markersGS <- getMarkerFeatures(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
for (i in seq_along(markerList)) {
  df <- as.data.frame(markerList[[i]])
  write.table(df, paste0("markerList_root_resolut0.8_", i, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
markerGenes  <- c("AT1G26870","AT1G61590","AT2G48130","AT4G25250")
p1 <- plotGroups(ArchRProj = projHeme2, groupBy = "Clusters", colorBy = "GeneScoreMatrix", name = markerGenes, title = "Marker GeneScoreMatrix by Clusters")
dev.off()
pdf("Plot-UMAPHarmony-Marker-Genes-W-Imputation.pdf",width = 10,height = 10)
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],nrow=2)
dev.off()
heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", labelMarkers = markerGenes, transpose = TRUE)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-C2", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)
dev.off()
projHeme2 <- addImputeWeights(projHeme2)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAPHarmony", imputeWeights = getImputeWeights(projHeme2))
pdf("Plot-UMAPHarmony-Marker-Genes-W-Imputation-C2.pdf",width = 35,height = 35)
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[9]],p1[[10]],p1[[11]],p1[[12]],p1[[13]],p1[[14]],p1[[15]],p1[[16]],p1[[17]],p1[[18]],p1[[19]],p1[[20]],p1[[21]],p1[[22]],p1[[23]],p1[[24]],p1[[25]],nrow=5)
dev.off()
p3 <- plotBrowserTrack(ArchRProj = projHeme2, groupBy = "Clusters", geneSymbol = markerGenes, upstream = 30000, downstream = 30000)
pdf("Plot-Tracks-Marker-Genes.pdf",width = 55,height = 55)
grid.arrange(p3[[1]],p3[[2]],p3[[3]],p3[[4]],p3[[5]],p3[[6]],p3[[7]],p3[[8]],p3[[9]],p3[[10]],p3[[11]],p3[[12]],p3[[13]],p3[[14]],p3[[15]],p3[[16]],p3[[17]],p3[[18]],p3[[19]],p3[[20]],p3[[21]],p3[[22]],p3[[23]],p3[[24]],p3[[25]],nrow=5)
dev.off()

#====================================================================================================================#
#                                         单细胞RNA数据细胞类型映射                                               ####
#====================================================================================================================#
seRNA <- readRDS("Silique_final202502_scRNA.rds")
seRNA <- readRDS("../root_final202502.rds")
# seRNA <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202407.At_PCA/2.Cluster2024/snRNA_202502.summary/rds.rmChrMC_cluster.20250307/silique_final202502.rds")
projHeme2 <- loadArchRProject("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/Save-ProjHeme-filter/")
# 无约束整合
projHeme2 <- addGeneIntegrationMatrix(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "Harmony", seRNA = seRNA_6, addToArrow = FALSE, force = TRUE, groupRNA = "cell_type", nameCell = "predictedCell_1.2", nameGroup = "predictedGroup_1.2", nameScore = "predictedScore_1.2", threads = 1)

# 约束整合
# cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_1.2))
# preClust <- colnames(cM)[apply(cM, 1 , which.max)]
# cbind(preClust, rownames(cM)) 
# unique(unique(projHeme2$predictedGroup_1.2))

# 画细胞类型散点图
pal <- paletteDiscrete(values = names(table(seRNA$cell_type)))
cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_1.2))
row_percent <- cM / rowSums(cM)
result <- data.frame(RowName = rownames(row_percent), MaxColumn = colnames(row_percent)[apply(row_percent, 1, which.max)], MaxPercentage = apply(row_percent, 1, max))
p1 <- plotEmbedding(projHeme2, colorBy = "cellColData", embedding = "UMAPHarmony", name = "Clusters")
dev.off()
p2 <- plotEmbedding(projHeme2, colorBy = "cellColData", name = "predictedGroup_1.2", embedding = "UMAPHarmony", pal = pal)
dev.off()
p3 <- DimPlot(seRNA, group.by = "cell_type", reduction = "umap", label = TRUE, pt.size = 1) + scale_color_manual(values = pal) + NoLegend()
pdf("Plot-UMAP-RNA-cellclass-Integration-root-0316-1.8.pdf",width = 20,height = 12)
grid.arrange(p1, p2, p3, nrow=1)
dev.off()
pdf("Plot-UMAP-RNA-cellclass-Integration-silique.pdf",width = 20,height = 12)
grid.arrange(p1, p2,p3, nrow=1)
dev.off()
pdf("Plot-UMAP-RNA-cellclass-Integration-silique2.pdf",width = 6,height = 6)
grid.arrange(p3,nrow=1)
dev.off()
# 统计每个sample中的细胞类型的占比
data <- as.data.frame(cbind(projHeme2$Sample,projHeme2$predictedGroup_1.2))
data <- as.data.frame(cbind(seRNA$orig.ident,seRNA$cell_type))
colnames(data) <- c("sample","cell_type")
df_plot <- data %>% group_by(sample, cell_type) %>% summarise(count = n(), .groups = "drop") %>% group_by(sample) %>% mutate(percentage = count / sum(count) * 100)
# 画堆叠柱状图
pal <- paletteDiscrete(values = names(table(seRNA$cell_type)))
p <- ggplot(df_plot, aes(x = sample, y = percentage, fill = cell_type)) + geom_bar(stat = "identity") + scale_y_continuous(labels = scales::percent_format(scale = 1)) + scale_fill_manual(values = pal) + labs(x = "Sample", y = "Percentage", fill = "Cell Type") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf("Plot-UMAP-RNA-cellproportion-Integration-silique.pdf",width = 6,height = 6)
grid.arrange(p,nrow=1)
dev.off()
pdf("Plot-UMAP-RNA-cellproportion-Integration-root-seRNA.pdf",width = 6,height = 6)
grid.arrange(p,nrow=1)
dev.off()
pdf("Plot-UMAP-RNA-cellproportion-Integration-root-predictedScore_1.2.pdf",width = 6,height = 6)
grid.arrange(p1,nrow=1)
dev.off()
# 保存
saveRDS(projHeme2, "PCA_snATAC-Integration-root-0316-1.8.rds")

# 可视化marker基因的散点图
options(mc.cores = 1) 
seRNA <- readRDS("Silique_final202502_scRNA.rds")
# projHeme2 <- readRDS("PCA_snATAC_cluster_silique_resolut_0.8.rds")
projHeme2 <- loadArchRProject("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/Save-ProjHeme-filter/")
projHeme3 <- addImputeWeights(projHeme3)
markerGenes  <- c("AT1G05220", "AT1G03103", "AT1G04357", "AT1G01140")
p1 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "GeneIntegrationMatrix", name = markerGenes, continuousSet = "horizonExtra", embedding = "UMAPHarmony", imputeWeights = getImputeWeights(projHeme3))
p2 <- plotEmbedding(ArchRProj = projHeme3, colorBy = "GeneScoreMatrix", continuousSet = "horizonExtra", name = markerGenes, embedding = "UMAPHarmony", imputeWeights = getImputeWeights(projHeme3))
# 用cowplot将这些标记基因绘制在一起
p1c <- lapply(p1, function(x){x + guides(color = FALSE, fill = FALSE) + theme_ArchR(baseSize = 6.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())})
p2c <- lapply(p2, function(x){x + guides(color = FALSE, fill = FALSE) + theme_ArchR(baseSize = 6.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())})
dev.off()
pdf("Plot-UMAP-RNA-cellclass-Integration-GeneIntegrationMatrix-silique1.pdf",width = 7.5,height = 6)
grid.arrange(p1c[[1]],p1c[[2]],p1c[[3]],p1c[[4]],p1c[[5]],nrow=2)
dev.off()
pdf("Plot-UMAP-RNA-cellclass-Integration-GeneIntegrationMatrix-silique2.pdf",width = 7.5,height = 6)
grid.arrange(p2c[[1]],p2c[[2]],p2c[[3]],p2c[[4]],p2c[[5]],nrow=2)
dev.off()

# 保存
projHeme3$Clusters2 <- projHeme3$predictedGroup
saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3-Ingegration-Silique", load = FALSE)

#====================================================================================================================#
#                                       合并不同组织(root,silique)数据                                             ####
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

#====================================================================================================================#
#                               peak鉴定及差异分析 & 差异peak中的motif富集分析                                    ####
#====================================================================================================================#
# 单独根据每一组细胞（例如聚类）鉴定peak，这样可以分析出不同组的特异性peak
# root-silique
proj <- readRDS("PCA_snATAC_root_silique_celltype.rds")
# 拟混池重复
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "celltype", threads = 1)
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "celltype", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
proj <- addPeakMatrix(proj)
# getAvailableMatrices(proj)
# 另一种方法鉴定Peaks w/TileMatrix
# projHemeTmp <- addReproduciblePeakSet(ArchRProj = projHeme4, groupBy = "Clusters2", peakMethod = "Tiles", method = "p")
# 保存
saveRDS(proj, "PCA_snATAC_root_silique_peak.rds")
projHeme5 <- proj
# 计算每组分类的特异性peak并画热图
markersPeaks <- getMarkerFeatures(ArchRProj = projHeme5, useMatrix = "PeakMatrix", groupBy = "celltype", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
pdf("Plot-peaksMatrix-root_silique1.pdf", width = 8, height = 6)
markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", transpose = TRUE)
dev.off()
# returnGR=TRUE, 可以用getMarkers()再返回一个GRangesList
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR=TRUE)
# 画某一个分类（Cotyledons）的特异性peak
pdf("Plot-peaksMatrix--root_silique2.pdf",width = 10,height = 6)
markerPlot(seMarker = markersPeaks, name = "Cotyledons", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
markerPlot(seMarker = markersPeaks, name = "Cotyledons", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
dev.off()
markerGenes  <- c("AT1G05220", "AT1G03103", "AT1G04357", "AT1G01140")
# 画基因的peak track
p <- plotBrowserTrack(ArchRProj = projHeme5, groupBy = "celltype", geneSymbol = markerGenes, features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Cotyledons"], upstream = 20000, downstream = 20000)
pdf("Plot-peaksMatrix--root_silique3.pdf",width = 13,height = 30)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],nrow=5)
dev.off()
# 画两个分类之间（Cotyledons和Trichoblast）的差异peak
markerTest <- getMarkerFeatures(ArchRProj = projHeme5, useMatrix = "PeakMatrix", groupBy = "celltype", testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), useGroups = "Cotyledons", bgdGroups = "Trichoblast")
pdf("Plot-peaksMatrix--root_silique4.pdf",width = 12,height = 5)
markerPlot(seMarker = markerTest, name = "Cotyledons", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
markerPlot(seMarker = markerTest, name = "Cotyledons", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
dev.off()

# 添加motif
projHeme5 <- readRDS("PCA_snATAC_root_silique_peak.rds")
meme_file <- "../motif_AT/Ath_TF_binding_motifs.meme" 
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
# peaks <- getPeakSet(projHeme5)
# motif_hits <- matchMotifs(pwmlist, peaks, genome = BSgenome.Athaliana.TAIR10.02022009)
# motif_matrix <- assay(motif_hits, "motifMatches")
# overlapMotifs <- motifmatchr::matchMotifs(pwmlist, peaks, genome = "BSgenome.Athaliana.TAIR10.02022009")
# colnames(overlapMotifs) <- paste0("Motif_", seq_len(ncol(overlapMotifs)))
projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifPWMs = pwmlist, name = "Motif")

# Motif偏离，不基于细胞分类情况
projHeme5 <- addBgdPeaks(projHeme5)
# 根据所有的motif注释计算每个细胞的偏离值
projHeme5 <- addDeviationsMatrix(ArchRProj = projHeme5, peakAnnotation = "Motif", force = TRUE)
# 保存
saveRDS(projHeme5, "PCA_snATAC_root_silique_motif.rds")
# getAvailableMatrices(projHeme5)
projHeme5 <- readRDS("PCA_snATAC_root_silique_motif.rds")

plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
dev.off()
pdf("motif1-root-0316.pdf",width = 5,height = 5)
grid.arrange(plotVarDev,nrow=1)
dev.off()

# 画motif富集热图
motif_matrix <- TransferDataFromProject(ArchRProj = projHeme5, useMatrix = "MotifMatrix")
deviations <- assay(motif_matrix)
celltype <- projHeme5$celltype
deviation_by_celltype <- sapply(unique(celltype), function(ct) {
  cells <- which(celltype == ct)
  rowMeans(deviations[, cells])
})
mat_z <- t(scale(t(deviation_by_celltype)))
pdf("motif_heatmap1.pdf",width = 30,height = 20)
Heatmap(mat_z)
dev.off()

# 把motif聚类成家族进行富集分析
match_TF <- read.table("../motif_AT/TF_family.txt", header=T, stringsAsFactors = F)
match_TF <- match_TF[match_TF$TF %in% rownames(deviations), ]
deviation_mat_filtered <- deviations[match_TF$TF, ]
rownames(match_TF) <- match_TF$TF  
TF_family <- match_TF[rownames(deviation_mat_filtered), "Family"]
celltypes <- projHeme5$celltype
celltype_levels <- unique(celltypes)
dev_by_celltype <- sapply(celltype_levels, function(ct) {
  cells <- projHeme5$cellNames[projHeme5$celltype == ct]
  rowMeans(deviation_mat_filtered[, cells, drop=FALSE])
})
rownames(dev_by_celltype) <- rownames(deviation_mat_filtered)
dev_family <- as.data.frame(dev_by_celltype)
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

# 提取每种细胞类型的前2000个ACRs进行motif富集分析
# 提取每个分类的特异peak
markersPeaks <- getMarkerFeatures(ArchRProj = projHeme5, useMatrix = "PeakMatrix", groupBy = "celltype", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
top2000_list <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
# 提取 GRanges 对象
top2000_gr_list <- lapply(top2000_list, function(df) {
  head(df[order(df$Log2FC, decreasing = TRUE), ], 2000)
})

top2000_gr_list2 <- lapply(top2000_gr_list, function(df) {
  GRanges(seqnames = df$seqnames,
          ranges = IRanges(start = df$start, end = df$end))
})

enrichMotifs <- lapply(top2000_gr_list2, function(gr) {
  peakAnnoEnrichment(
    seMarker = gr,             
    ArchRProj = projHeme5,     
    peakAnnotation = "Motif",   
    cutOff = "FDR <= 0.1 & Log2FC >= 1"  
  )
})

enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = projHeme5, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 1")
motif_mat <- getMatrixFromProject(projHeme5, useMatrix = "MotifMatrix")
deviation_mat <- assay(motif_mat)
motif_peaks <- rowRanges(motif_mat)

# 计算上述每个分类的特异peak的motif富集并画图
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = projHeme5, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 1")
pdf("motif_enrichment.pdf",width = 10,height = 5)
plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
dev.off()

# 提取两个分类之间（前景Cotyledons和背景Trichoblast）的差异peak
markerTest <- getMarkerFeatures(ArchRProj = projHeme5, useMatrix = "PeakMatrix", groupBy = "celltype", testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), useGroups = "Cotyledons", bgdGroups = "Trichoblast")
# 计算上述peak中的前景分类（Cotyledons）中显著上调的motif并画图
motifsUp <- peakAnnoEnrichment(seMarker = markerTest, ArchRProj = projHeme5, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + geom_point(size = 1) + ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), size = 1.5, nudge_x = 2, color = "black") + theme_ArchR() + ylab("-log10(P-adj) Motif Enrichment") + xlab("Rank Sorted TFs Enriched") + scale_color_gradientn(colors = paletteContinuous(set = "comet"))
pdf("ggUp.pdf",width = 8,height = 8)
grid.arrange(ggUp,nrow=1)
dev.off()
# 计算上述peak中的前景分类（Cotyledons）中显著下调的motif并画图
motifsDo <- peakAnnoEnrichment(seMarker = markerTest, ArchRProj = projHeme5, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC <= -0.5")
df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + geom_point(size = 1) + ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), size = 1.5, nudge_x = 2, color = "black") + theme_ArchR() + ylab("-log10(FDR) Motif Enrichment") + xlab("Rank Sorted TFs Enriched") + scale_color_gradientn(colors = paletteContinuous(set = "comet"))
pdf("ggUpDo.pdf",width = 10,height = 5)
grid.arrange(ggUp,ggDo,nrow=1)
dev.off()

# 提取特异motif进行分析
motifs <- c("Motif_322", "Motif_406", "Motif_102", "Motif_128")
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
# 筛选motif
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "deviations:Motif_322"]
markerMotifs <- markerMotifs[c(1,2,3,14,15,16)]
# 画山脊图
p <- plotGroups(ArchRProj = projHeme5, groupBy = "Clusters2", colorBy = "MotifMatrix", name = markerMotifs, imputeWeights = getImputeWeights(projHeme5))
dev.off()
pdf("V2-motif2-root-0316.pdf",width = 30,height = 20)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],nrow=1)
dev.off()
# 画散点图
p3 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "MotifMatrix", name = sort(markerMotifs), embedding = "UMAPHarmony", imputeWeights = getImputeWeights(projHeme5))
dev.off()
pdf("V2-motif4-root-0316.pdf",width = 20,height = 14)
grid.arrange(p3[[1]],p3[[2]],p3[[3]],p3[[4]],p3[[5]],p3[[6]],nrow=2)
dev.off()
# 后面不知道是画什么图，也可以使用GeneIntegrationMatrix
markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
p <- plotEmbedding(ArchRProj = projHeme5, colorBy = "GeneScoreMatrix", name = sort(markerRNA), embedding = "UMAPHarmony", continuousSet = "blueYellow", imputeWeights = getImputeWeights(projHeme5))

#====================================================================================================================#
#                                               motif足迹分析                                                     ####
#====================================================================================================================#
projHeme5 <- readRDS("PCA_snATAC_motifs-root-0316.rds")
motifPositions <- getPositions(projHeme5)
# 提取部分感兴趣的TF motifs用于展示
motifs <- c("MA1236.1.ERF003","MA1275.1.DOF1.6","MA1210.2.HAT22")
# markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
# markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
# markerMotifs <- markerMotifs[c(1,2,3,14,15,16)]
# 为了准确找到TF足迹，我们需要大量的reads。因此，细胞需要进行分组生成拟混池ATAC-seq谱才能用于TF足迹分析。这些拟混池谱之前在peak鉴定时就已经保存为分组覆盖文件。 如果没有在ArchRProject添加分组覆盖信息，则运行如下命令
# projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2")
# 通过positions参数来选择motif
seFoot <- getFootprints(ArchRProj = projHeme5, positions = motifPositions[motifs], groupBy = "Clusters2")
# 当我们获取了这些足迹，我们可以使用plotFootprints()函数进行展示。该函数能够同时以多种方式对足迹进行标准化。
# Tn5偏好的足迹标准化
# 减去Tn5偏好
pdf("root-0316_foot1.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, plotName = "Footprints-Subtract-Bias", ArchRProj = projHeme5, normMethod = "Subtract", addDOC = FALSE, smoothWindow = 5)
dev.off()
# 这些图保存在ArchRProject的outputDirectory。如果需要绘制所有motif, 可以将其返回为ggplot2对象，需要注意这个ggplot对象会非常大
# 除以Tn5偏好
pdf("root-0316_foot2.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, ArchRProj = projHeme5, normMethod = "Divide", plotName = "Footprints-Divide-Bias", addDOC = FALSE, smoothWindow = 5)
dev.off()
# 无Tn5偏好标准化的足迹
pdf("root-0316_foot3.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot, ArchRProj = projHeme5,normMethod = "None", plotName = "Footprints-No-Normalization", addDOC = FALSE, smoothWindow = 5)
dev.off()
# 与TSS的距离？
seTSS <- getFootprints(ArchRProj = projHeme5, positions = GRangesList(TSS = getTSS(projHeme5)), groupBy = "Clusters2", flank = 2000)

#====================================================================================================================#
#                                        整合分析:peak关联gene调控                                                ####
#====================================================================================================================#
# 可以只用ATAC-seq数据进行分析，如识别peak之间的共开放性来预测调控相互作用，也可以整合scRNA-seq数据，如通过peak-基因的连锁分析预测增性子活性
# 查看peak之间的共开放相关性
projHeme5 <- addCoAccessibility(ArchRProj = projHeme5, reducedDims = "IterativeLSI")
cA <- getCoAccessibility(ArchRProj = projHeme5, corCutOff = 0.5, resolution = 1, returnLoops = FALSE)
# metadata(cA)[[1]]
# returnLoops的参数含义是什么？resolution的含义是什么？
cA <- getCoAccessibility(ArchRProj = projHeme5, corCutOff = 0.5, resolution = 1, returnLoops = TRUE)
cA <- getCoAccessibility(ArchRProj = projHeme5, corCutOff = 0.5, resolution = 1000, returnLoops = TRUE)
cA <- getCoAccessibility(ArchRProj = projHeme5, corCutOff = 0.5, resolution = 10000, returnLoops = TRUE)
# 在browser track中绘制基因区域在不同细胞类型中的共开放
markerGenes  <- c("AT1G05220", "AT1G03103","AT1G04357","AT1G01140")
p <- plotBrowserTrack(ArchRProj = projHeme5, groupBy = "Clusters2", geneSymbol = markerGenes, upstream = 50000, downstream = 50000, loops = getCoAccessibility(projHeme5))
pdf("root-Tracks-Marker-Genes-with-CoAccessibility.pdf",width = 25,height = 10)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],nrow=2)
dev.off()

# Peak2GeneLinkages peak关联基因分析
projHeme5 <- addPeak2GeneLinks(ArchRProj = projHeme5, reducedDims = "IterativeLSI")
p2g <- getPeak2GeneLinks(ArchRProj = projHeme5, corCutOff = 0.45, resolution = 1, returnLoops = FALSE)
# metadata(p2g)[[1]]
p2g <- getPeak2GeneLinks(ArchRProj = projHeme5, corCutOff = 0.45, resolution = 10000, returnLoops = TRUE)
# 在browser track中绘制每种细胞类型的peak-to-gene连接
markerGenes  <- c("AT1G05220", "AT1G03103", "AT1G04357", "AT1G01140")
p <- plotBrowserTrack(ArchRProj = projHeme5, groupBy = "Clusters2", geneSymbol = markerGenes, upstream = 50000, downstream = 50000, loops = getPeak2GeneLinks(projHeme5))
pdf("root-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf",width = 25,height = 10)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],nrow=2)
dev.off()

# 绘制Peak-to-gene连接热图
p <- plotPeak2GeneHeatmap(ArchRProj = projHeme5, groupBy = "Clusters2")
pdf("root-Tracks-Marker-Genes-with-Peak2Gene-heatmap.pdf",width = 25,height = 10)
# grid.arrange(p[[1]],p[[2]],nrow=1)
draw(p)
dev.off()
# 保存
saveRDS(projHeme5, "PCA_snATAC_foot-root-0316.rds")
saveArchRProject(ArchRProj = a, outputDirectory = "Save-ProjHeme-foot-root-0316", load = FALSE)

# 正向TF-调控因子鉴定
projHeme5 <- readRDS("PCA_snATAC_foot-root-0316.rds")
# 第一步: 鉴定偏离TF motif
seGroupMotif <- getGroupSE(ArchRProj = projHeme5, useMatrix = "MotifMatrix", groupBy = "Clusters2")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){rowMaxs(assay(seZ) - assay(seZ)[,x])}) %>% Reduce("cbind", .) %>% rowMaxs
# 第二步: 鉴定相关的TF motif和TF基因得分/表达值
corGSM_MM <- correlateMatrices(ArchRProj = projHeme5, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "MotifMatrix", reducedDims = "IterativeLSI")
corGIM_MM <- correlateMatrices(ArchRProj = projHeme5, useMatrix1 = "GeneIntegrationMatrix", useMatrix2 = "MotifMatrix", reducedDims = "IterativeLSI")
# 第三步: 在相关性DataFrame中添加极大偏差值
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
# 第四步: 鉴定正向TF调控因子并画图
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
p_GSM <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) + geom_point() + theme_ArchR() + geom_vline(xintercept = 0, lty = "dashed") + scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) + xlab("Correlation To Gene Score") + ylab("Max TF Motif Delta") + scale_y_continuous(expand = c(0,0), limits = c(0, max(corGSM_MM$maxDelta)*1.05))

corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
p_GIM <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) + geom_point() + theme_ArchR() + geom_vline(xintercept = 0, lty = "dashed") + scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) + xlab("Correlation To Gene Expression") + ylab("Max TF Motif Delta") + scale_y_continuous(expand = c(0,0), limits = c(0, max(corGIM_MM$maxDelta)*1.05))

# 画细胞分化轨迹
markersGS <- getMarkerFeatures(ArchRProj = projHeme5, useMatrix = "GeneScoreMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
projHeme5 <- addTrajectory(ArchRProj = projHeme5, name = "MyeloidU", groupBy = "Clusters2", trajectory = trajectory, embedding = "UMAPHarmony", force = TRUE)
head(projHeme5$MyeloidU[!is.na(projHeme5$MyeloidU)])
# [1] 46.07479 53.75027 44.82834 43.18828 47.49617 43.21015
p <- plotTrajectory(projHeme5, trajectory = "MyeloidU", colorBy = "cellColData", name = "MyeloidU")
pdf("root-Trajectory.pdf",width = 20,height = 10)
grid.arrange(p[[1]],p[[2]],nrow=1)
dev.off()

#====================================================================================================================#
#                                       使用marker基因进行细胞注释                                                ####
#====================================================================================================================#
projHeme2 <- loadArchRProject("../Save-ProjHeme3-Ingegration-root")
projHeme2 <- readRDS("PCA_snATAC-Integration-root-0316-1.2.rds")
projHeme2 <- saveRDS("PCA_snATAC-Integration-root-0316-1.2-mergeCell.rds")
seRNA <- readRDS("../root_final202502.rds")
seRNA <- readRDS("../root_final202502.rds")
seRNA <- readRDS("../../Cell/D6_root_rename.slim.RDS")
# 使用不同分辨率的聚类方式，和snRNA的数据进行整合
projHeme2 <- readRDS("../rds/PCA_snATAC_root")
# 此为无约束整合，使用RNA中注释出来的cell_type类型
projHeme2 <- addGeneIntegrationMatrix(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "Harmony", seRNA = seRNA, addToArrow = FALSE, force = TRUE, groupRNA = "cell_type", nameCell = "predictedCell", nameGroup = "predictedGroup", nameScore = "predictedScore", threads = 1)
# 通过每个分类的基因得分鉴定标记基因
markersGS <- getMarkerFeatures(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# 跟已知的marker基因进行overlap，已知marker基因来源于项目组收集以及cell文章使用的
root_marker <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/annotation/root_marker_merged.txt", header=T, sep="\t", stringsAsFactors = F)
root_marker$Validatedmarkergenes <- root_marker$GeneID
data <- data.frame(name = character(), seqnames = character(), start = integer(), end = integer(), strand = integer(), idx = integer(), Log2FC = numeric(), FDR = numeric(), MeanDiff = numeric(), cluster = numeric(), GeneName = character(), MarkerType = character(), References = character(), GeneSymbol = character(), CellType = character(), Organ = character(), Notes = character(), Validatedmarkergenes = character(), stringsAsFactors = FALSE)

for (i in seq_along(markerList)) {
  if(dim(markerList[[i]])[1] != 0){
    df <- as.data.frame(markerList[[i]])
    df$cluster <- i
    merged_data <- merge(df, root_marker, by.x = "name", by.y = "GeneID", all.x = TRUE)
    data <- rbind(data,merged_data)}}
# 写出文件
write.table(data, "root-0316_markerGene_solut1.8.txt", sep="\t",quote=F,row.names = F)
# 保存
saveRDS(projHeme2, "PCA_snATAC-Integration-root-0316-1.8.rds")

# 画图验证是否这些基因是在某个类型里面特异表达
projHeme2 <- addImputeWeights(projHeme2)
marker <- read.table("root_Validatedmarkergenes.txt", header=F,stringsAsFactors=F)
marker <- read.table("root_Validatedmarkergenes_resolu1.8.txt", header=F,stringsAsFactors=F)
markerGenes  <- marker[,1]

# 基因表达和基因得分在不同细胞类型中的散点图分布
# 基因表达量值来自于GeneIntegrationMatrix
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "GeneIntegrationMatrix", name = markerGenes, continuousSet = "horizonExtra", embedding = "UMAPHarmony", imputeWeights = getImputeWeights(projHeme2))
# 使用GeneScoreMatrix里的基因得分值
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "GeneScoreMatrix", continuousSet = "horizonExtra", name = markerGenes, embedding = "UMAPHarmony", imputeWeights = getImputeWeights(projHeme2))
pdf("Plot-UMAP-RNA-cellclass-Integration-GeneIntegrationMatrix-silique1.pdf",width = 7.5,height = 6)
grid.arrange(p1c[[1]],p1c[[2]],p1c[[3]],p1c[[4]],p1c[[5]],nrow=2)
dev.off()
pdf("Plot-UMAP-RNA-cellclass-Integration-GeneIntegrationMatrix-silique2.pdf",width = 7.5,height = 6)
grid.arrange(p2c[[1]],p2c[[2]],p2c[[3]],p2c[[4]],p2c[[5]],nrow=2)
dev.off()

# mark <- gene_scores[markerGenes,]
# data1 <- cbind(as.data.frame(t(mark)),projHeme2$Clusters, projHeme2$predictedGroup)
# 基因表达和基因得分在不同细胞类型中的箱线图分布
se3 <- getMatrixFromProject(ArchRProj = projHeme2, useMatrix = "GeneIntegrationMatrix")
gene_scores_2 <- assay(se3, "GeneIntegrationMatrix")
rownames(gene_scores_2) <- rowData(se3)$name
mark_2 <- gene_scores_2[markerGenes,]
data2 <- cbind(as.data.frame(t(mark_2)),projHeme2$Clusters, projHeme2$predictedGroup)
p3 <- plotGroups(ArchRProj = projHeme2, groupBy = "Clusters", colorBy = "GeneIntegrationMatrix", name = markerGenes, plotAs = "violin", title = "Marker GeneIntegrationMatrix by Clusters")
p4 <- plotGroups(ArchRProj = projHeme2, groupBy = "Clusters", colorBy = "GeneScoreMatrix", name = markerGenes, plotAs = "violin", title = "Marker GeneScoreMatrix by Clusters")
# 把p3和p4合并起来画图
pdf("combined_marker_1.8_plots.pdf", width = 10, height = 8)
# 计算需要多少行
n_plots <- 98
plots_per_row <- 2
n_rows <- ceiling(n_plots / plots_per_row)
# 循环绘制每一行的图形
for (i in 1:n_rows) {
  start_idx <- (i - 1) * plots_per_row + 1
  end_idx <- min(i * plots_per_row, n_plots)
  plots_row <- c(p3[start_idx:end_idx], p4[start_idx:end_idx])
  do.call(grid.arrange, c(plots_row, ncol = plots_per_row))}
dev.off()

# 其他的检验方式（细胞类型相关性）
projHeme2 <- readRDS("PCA_snATAC-Integration-root-0316-1.2.rds")
seRNA <- readRDS("../root_final202502.rds")
# 提取对应矩阵里面的数据（GeneScoreMatrix或者GeneIntegrationMatrix）画细胞类型的相关性热图，整合RNA数据之后再画一次相关性
se2 <- getMatrixFromProject(ArchRProj = projHeme2, useMatrix = "GeneIntegrationMatrix")
gene_mat <- assay(se2, "GeneIntegrationMatrix")
rownames(gene_mat) <- rowData(se2)$name
cells <- colnames(se2)
colnames(gene_mat) <- cells
gene_means <- rowMeans(gene_mat)
names(gene_means) <- rownames(gene_mat)
top_genes <- names(sort(gene_means, decreasing = TRUE))[1:3000]
top_expr_mat <- gene_mat[top_genes, ]
# 画样本的相关性
sample_ids <- projHeme2$Sample
# 画细胞类型的相关性
sample_ids <- projHeme2$Clusters
sample_expr <- sapply(unique(sample_ids), function(s) {
  col_idx <- which(sample_ids == s)
  Matrix::rowMeans(top_expr_mat[, col_idx, drop = FALSE])})
cor_mat <- cor(sample_expr, method = "spearman")
pdf("Plot-scATAC-root-0316-1.2-sample-heatmap.pdf",width = 8, height = 8)
pheatmap(cor_mat, main = "Sample Correlation (Top 5000 Genes)", cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 10, fontsize_col = 10)
dev.off()
pdf("Plot-scATAC-root-0316-1.2-celltype-heatmap-GeneIntegrationMatrix.pdf",width = 8, height = 8)
pheatmap(cor_mat, main = "Cell type Correlation (Top 5000 Genes)", cluster_rows = TRUE, cluster_cols = TRUE, fontsize_row = 10, fontsize_col = 10)
dev.off()

#====================================================================================================================#
#                                    比较我们的RNA注释和Cell的RNA注释                                             ####
#====================================================================================================================#
seRNA <- readRDS("../root_final202502.rds")
seRNA$tech <- "RNA"
seRNA$type <- seRNA$cell_type
projHeme2 <- readRDS("PCA_snATAC-Integration-root-0316-1.2-mergeCell.rds")
# 提取GeneScoreMatrix或者GeneIntegrationMatrix的数据
geneScoreMat <- getMatrixFromProject(ArchRProj = projHeme2, useMatrix = "GeneIntegrationMatrix")
geneScoreMat_sparse <- assay(geneScoreMat)
rownames(geneScoreMat_sparse) <- rowData(geneScoreMat)$name
cells <- colnames(geneScoreMat_sparse)
geneScore_seurat <- CreateAssayObject(counts = geneScoreMat_sparse)
seATAC <- CreateSeuratObject(counts = geneScoreMat_sparse)
seATAC[["ATAC"]] <- CreateAssayObject(counts = geneScoreMat_sparse)
seATAC$tech <- "ATAC"
seATAC$type <- projHeme2$predictedGroup
combined <- merge(seRNA, y = seATAC, add.cell.ids = c("RNA", "ATAC"))
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- FindMultiModalNeighbors(combined, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 1:30))
combined <- RunUMAP(combined, dims = 1:30)

pdf("Plot-UMAP-root-0316-RNA-ATAC-geneScore.pdf",width = 12,height = 6)
DimPlot(combined, group.by = "tech", reduction = "umap", pt.size = 0.5)
dev.off()
pdf("Plot-UMAP-root-0316-RNA-ATAC-geneIntegration.pdf",width = 12,height = 6)
DimPlot(combined, group.by = "tech", reduction = "umap", pt.size = 0.5)
dev.off()
pdf("Plot-UMAP-root-0316-RNA-ATAC-geneIntegration-celltype.pdf",width = 12,height = 6)
DimPlot(combined, group.by = "type", reduction = "umap", pt.size = 0.5)
dev.off()

geneScoreMat <- getMatrixFromProject(ArchRProj = projHeme2, useMatrix = "GeneIntegrationMatrix")
geneScoreMat_sparse <- assay(geneScoreMat)
rownames(geneScoreMat_sparse) <- rowData(geneScoreMat)$name
cells <- colnames(geneScoreMat_sparse)
geneScore_seurat <- CreateAssayObject(counts = geneScoreMat_sparse)
gene.activities <- GeneActivity(seATAC)
seATAC[["ACTIVITY"]] <- CreateAssayObject(counts = geneScoreMat_sparse)
DefaultAssay(seATAC) <- "ACTIVITY"

# RNA 数据预处理
seRNA <- NormalizeData(seRNA) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 30)
# ATAC 的 GeneActivity 做同样处理
seATAC <- NormalizeData(seATAC) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 30)
anchors <- FindIntegrationAnchors(object.list = list(seRNA, seATAC), anchor.features = intersect(VariableFeatures(seRNA), VariableFeatures(seATAC)), dims = 1:30)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(combined) <- "integrated"
combined$cell_type <- ifelse(is.na(combined$cell_type), combined$predicted_celltype, combined$cell_type)
combined <- ScaleData(combined) %>% FindVariableFeatures() %>% RunPCA(npcs = 30) %>% RunUMAP(dims = 1:30)
combined$type_seq <- ifelse(combined$orig.ident %in% c("A220316001", "A220316002", "A220322003", "A220322004"), "seATAC", "seRNA")

# 绘制UMAP图，按cell_type分组
pdf("seRNA_seATAC_integration.pdf",width = 21,height = 7)
DimPlot(combined, group.by = "cell_type", pt.size = 1)
# 按数据来源显示不同颜色
DimPlot(combined, group.by = "type_seq", pt.size = 1)
dev.off()

seATAC <- TransferData(anchorset = anchors, refdata = seRNA$celltype, dims = 1:30)
combined$celltype <- ifelse(combined$orig.ident == "seRNA", combined$celltype, combined$predicted.id)
DimPlot(combined, group.by = "celltype", label = TRUE)

# 提取 GeneScoreMatrix
geneScoreMat <- getMatrixFromProject(projHeme2, useMatrix = "GeneScoreMatrix")
# 提取 counts 矩阵
geneScoreMat_sparse <- assay(geneScoreMat)
# 转为稀疏矩阵
geneScoreMat_sparse <- as(geneScoreMat_sparse, "dgCMatrix")
# 提取元数据
gene.names  <- rowData(geneScoreMat)$name
cell.names <- colnames(geneScoreMat_sparse)
cell.meta <- as.data.frame(colData(geneScoreMat))
rownames(cell.meta) <- cell.names
rownames(geneScoreMat_sparse) <- gene.names
seATAC <- CreateSeuratObject(counts = geneScoreMat_sparse, assay = "ATAC", meta.data = cell.meta)
seATAC$cell_type <- projHeme2@cellColData$predictedGroup_1.2
combined <- merge(seRNA, y = seATAC, add.cell.ids = c("RNA", "ATAC"))
anchors <- FindTransferAnchors(reference = seRNA, query = seATAC, dims = 1:30)
# 进行细胞类型预测
predictions <- TransferData(anchorset = anchors, refdata = seRNA$cell_type, dims = 1:30)
# 将预测的细胞类型添加到 seATAC
seATAC <- AddMetaData(seATAC, metadata = predictions)
# RNA
seRNA <- NormalizeData(seRNA) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# 用 GeneScoreMatrix 作为类似 expression 的活动矩阵
DefaultAssay(seATAC) <- "ATAC"
seATAC <- NormalizeData(seATAC) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# seATAC@reductions
transfer_anchors <- FindTransferAnchors(reference = seRNA, query = seATAC, dims = 1:30)
predicted.labels <- TransferData(anchorset = transfer_anchors, refdata = seRNA$cell_type, dims = 1:30)
transfer_anchors <- FindTransferAnchors(reference = seRNA, query = seATAC, dims = 1:30, k.anchor = 10, k.filter = 200 )
seATAC$predicted_celltype <- predicted.labels$predicted.id
anchors <- FindIntegrationAnchors(object.list = list(seRNA, seATAC), dims = 1:30)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined) %>% RunPCA(npcs = 30) %>% RunUMAP(dims = 1:30)

# 获取共同的细胞名
common_cells <- intersect(colnames(seRNA), colnames(seATAC))
# 从 RNA 数据中提取细胞类型标签
rna_celltypes <- seRNA$celltype[common_cells]
# 将 RNA 的细胞类型标签赋值给 ATAC 数据
seATAC$celltype <- NA
seATAC$celltype[common_cells] <- rna_celltypes
# 将细胞类型标签添加到 combined 的 meta.data
combined$celltype <- NA
combined$celltype[common_cells] <- rna_celltypes
# 按细胞类型显示 UMAP 图
DimPlot(combined, group.by = "celltype", reduction = "umap", label = TRUE)
# 如果需要分不同来源的模态显示，可以使用 orig.ident
DimPlot(combined, group.by = "orig.ident", reduction = "umap", label = TRUE)
# 如果需要进一步的细胞类型分析和聚类，使用 FindClusters
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined)
# 可视化聚类结果
DimPlot(combined, group.by = "seurat_clusters", reduction = "umap", label = TRUE)
# 下面的代码需要整理
seRNA <- readRDS("../root_final202502.rds")
seRNA_6 <- readRDS("../../Cell/D6_root_rename.slim.RDS")
seRNA_11 <- readRDS("../../Cell/D11_root_rename.slim.RDS")

#整合cell文章的root数据
seurat.list <- list(seRNA_6, seRNA_11)
seurat.list <- Map(function(x, name) {
  x$sample <- name
  RenameCells(x, add.cell.id = name)
}, seurat.list, c("6day", "11day"))

combined <- Reduce(function(x, y) merge(x, y), seurat.list)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

p <- DimPlot(combined, reduction = "umap", group.by = "integrated_annotation", split.by = "sample", label=TRUE)
p <- DimPlot(combined, reduction = "umap", group.by = "sample", label=TRUE)

ggsave("Cell_RNA_sample.pdf", plot = p, width = 12, height = 8)

seurat.list <- list(seRNA_my,seRNA,seRNA2)
seurat.list <- Map(function(x, name) {
  x$sample <- name
  RenameCells(x, add.cell.id = name)
}, seurat.list, c("My", "Cell", "Cell"))

combined <- Reduce(function(x, y) merge(x, y), seurat.list)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

p <- DimPlot(combined, reduction = "umap", group.by = "final_annotation", split.by = "sample", label=TRUE)
ggsave("RNA.pdf", plot = p, width = 12, height = 8)

sample1_obj <- subset(combined, subset = sample == "My")
sample2_obj <- subset(combined, subset = sample == "Cell")

combined$cell_label <- NA
combined$cell_label[combined$sample == "My"] <- tolower(combined$cell_type[combined$sample == "My"])
combined$cell_label[combined$sample == "Cell"] <- tolower(combined$integrated_annotation[combined$sample == "Cell"])

p <- DimPlot(combined, group.by = "cell_label", split.by = "sample")

# 创建颜色表
library(scales)
color_pool <- hue_pal()(length(all_labels))
names(color_pool) <- all_labels

# Sample 1 里实际出现的细胞类型
used_labels1 <- unique(tolower(sample1_obj$cell_type))
colors1 <- color_pool[used_labels1]

# Sample 2 的
used_labels2 <- unique(tolower(sample2_obj$integrated_annotation))
colors2 <- color_pool[used_labels2]

sample1_obj$cell_label <- factor(tolower(sample1_obj$cell_type), levels = all_labels)
sample2_obj$cell_label <- factor(tolower(sample2_obj$integrated_annotation), levels = all_labels)

p1 <- DimPlot(sample1_obj, group.by = "cell_label") + scale_color_manual(values = colors1) + ggtitle("My")
p2<- DimPlot(sample2_obj, group.by = "cell_label") + scale_color_manual(values = colors2) + ggtitle("Cell")
p3 <- DimPlot(combined, group.by = "sample") 
#减少异质性
anchors <- FindIntegrationAnchors(object.list = list(sample1_obj, sample2_obj), dims = 1:30)
combined <- IntegrateData(anchorset = anchors, dims = 1:30)
pdf("RNA.pdf",width = 25,height = 8)
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

#====================================================================================================================#
#                                     使用Seurat进行RNA和ATAC数据映射                                             ####
#====================================================================================================================#
projHeme2 <- readRDS("rds/ArchR_bent_anchor.rds")
seRNA <- readRDS("../root_final202502.rds")

pbmc.rna <- seRNA
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

geneScoreMat <- getMatrixFromProject(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix")
geneScoreMat_sparse <- assay(geneScoreMat)
rownames(geneScoreMat_sparse) <- rowData(geneScoreMat)$name
cells <- colnames(geneScoreMat_sparse)
pbmc.atac <- CreateSeuratObject(counts = geneScoreMat_sparse, assay = "ATAC")

# 归一化
pbmc.atac <- NormalizeData(pbmc.atac)
# 选高变基因
pbmc.atac <- FindVariableFeatures(pbmc.atac)
# 标准化
pbmc.atac <- ScaleData(pbmc.atac)
# PCA降维
pbmc.atac <- RunPCA(pbmc.atac, npcs = 30)

#pbmc.atac <- FindNeighbors(pbmc.atac, reduction = "pca", dims = 1:30)
#pbmc.atac <- FindClusters(pbmc.atac, resolution = c(0.2, 0.5, 0.8, 1,1.2,1.6, 2, 4))
#pdf("RNA-ATAC-integration-test8-pca.pdf",width = 16, height = 20)
#clustree(pbmc.atac, prefix = "ATAC_snn_res.")
#dev.off()

# 使用Lsi聚类
pbmc.atac <- RunSVD(pbmc.atac, assay = "ATAC", reduction.key = "LSI_", reduction.name = "lsi")
pbmc.atac <- FindNeighbors(pbmc.atac, reduction = "lsi", dims = 1:30)
pbmc.atac <- FindClusters(pbmc.atac, resolution = 8)

# 测试精确度
# pbmc.atac <- FindClusters(pbmc.atac, resolution = c(0.2, 0.5, 0.8, 1,1.2,1.6, 2, 4))
# pdf("RNA-ATAC-integration-test7-lsi.pdf",width = 16, height = 20)
# clustree(pbmc.atac, prefix = "ATAC_snn_res.")
# dev.off()

# 然后UMAP
pbmc.atac <- RunUMAP(pbmc.atac,  reduction = "lsi",dims = 1:30)
# DimPlot(pbmc.atac, reduction = "umap")

# pbmc.atac <- RunTFIDF(pbmc.atac)
# pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
# pbmc.atac <- RunSVD(pbmc.atac)
# pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# 画图
# p1 <- DimPlot(pbmc.rna, group.by = "cell_type", label = TRUE) + NoLegend() + ggtitle("RNA")
# p2 <- DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
# pdf("RNA-ATAC-integration-test1.pdf",width = 8, height = 6)
# grid.arrange(p1,p2,nrow=1)
# dev.off()

# 通过基因上下游的peak估计基因表达的活性，如果是GeneScoreMatrix转换过来，是GeneScoreMatrix就代表了基因表达活性
# quantify gene activity
# gene.activities <- GeneActivity(pbmc.atac, features = VariableFeatures(pbmc.rna))
# add gene activities as a new assay
# pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# normalize gene activities
# DefaultAssay(pbmc.atac) <- "ACTIVITY"
# pbmc.atac <- NormalizeData(pbmc.atac)
# pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

# 另一种方法
# pbmc.atac$tech <- "atac" 
# DefaultAssay(pbmc.atac) <- "ACTIVITY" 
# pbmc.atac <- FindVariableFeatures(pbmc.atac) 
# pbmc.atac <- NormalizeData(pbmc.atac) 
# pbmc.atac <- ScaleData(pbmc.atac) 

# pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = geneScoreMat_sparse)
# DefaultAssay(pbmc.atac) <- "ACTIVITY"

# 对原始的peak矩阵进行降维处理，如果是从archR的GeneScoreMatrix转换过来的Seurat对象，则不需要做LSI降维，直接使用PCA降维
# DefaultAssay(pbmc.atac) <- "ATAC" #Error: Cannot find assay ATAC
# VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 100)) 
# 使用LSI方法进行降维
# pbmc.atac <- RunLSI(pbmc.atac, n = 50, scale.max = NULL) # could not find function "RunLSI"
# pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 1:50) 
# pbmc.atac <- RunPCA(pbmc.atac, npcs = 30)
# pbmc.atac <- RunUMAP(pbmc.atac, dims = 1:30)

# 画图
# p1 <- DimPlot(pbmc.atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq") 
# p2 <- DimPlot(pbmc.rna, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq") 
# pdf("RNA-ATAC-integration-test2.pdf",width = 8, height = 6)
# grid.arrange(p1,p2,nrow=1)
# dev.off()

# 识别anchor，使用全部的细胞和特征跑起来非常慢，随机抽取10000个细胞进行测试
set.seed(123) 
pbmc.rna_sub <- subset(pbmc.rna, cells = sample(colnames(pbmc.rna), 10000))
pbmc.atac_sub <- subset(pbmc.atac, cells = sample(colnames(pbmc.atac), 10000))
VariableFeatures(pbmc.rna_sub) <- head(VariableFeatures(pbmc.rna), 1500)

transfer.anchors <- FindTransferAnchors(reference = pbmc.rna_sub, query = pbmc.atac_sub, features = VariableFeatures(object = pbmc.rna_sub), reference.assay = "RNA", query.assay = "ATAC", reduction = "cca")

# 使用TransferData函数进行数据转移映射
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$cell_type, weight.reduction = pbmc.atac[["pca"]], dims = 2:30) 
# 添加预测信息 
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions) 
# 映射之后进行评估
# 方式一
pdf("RNA-ATAC-integration-test3.pdf",width = 8, height = 6)
hist(pbmc.atac$prediction.score.max) 
abline(v = 0.5, col = "red") 
dev.off()
table(pbmc.atac$prediction.score.max > 0.5) 

# 数据可视化: 画映射RNA注释之后的ATAC数据集
pbmc.atac.filtered <- subset(pbmc.atac, subset = prediction.score.max > 0.5) 
pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = unique(pbmc.rna$cell_type))
pdf("RNA-ATAC-integration-test4.pdf",width = 8, height = 6)
DimPlot(pbmc.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + NoLegend() + scale_colour_hue(drop = FALSE) 
DimPlot(pbmc.rna, group.by = "cell_type", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + NoLegend() 
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
coembed$celltype <- ifelse(!is.na(coembed$cell_type), coembed$cell_type, coembed$predicted.id) 
coembed$tech <- ifelse(!is.na(coembed$cell_type), "RNA", "ATAC") 

pdf("RNA-ATAC-integration-test5.pdf",width = 10, height = 6)
DimPlot(coembed, group.by = "tech") 
DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE) 
DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() 
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

#====================================================================================================================#
#                               Seurat对象转换ArchR对象画marker基因表达分布                                       ####
#====================================================================================================================#
projHeme2 <- readRDS("PCA_snATAC-Integration-root-0316-1.2.rds")
inter <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/annotation/PCA_snATAC_root_combined.rds")
# 提取Seurat的meta信息
meta_info <- inter@meta.data
# 确认名字对得上（非常重要）
all(rownames(meta_info) == projHeme2$cellNames)
# 应该是TRUE，如果是FALSE，需要按projHeme2$cellNames顺序重新排一下
# 如果需要重新排：
meta_info <- meta_info[match(projHeme2$cellNames, rownames(meta_info)), ]
# 把想要的列加到ArchR对象中
projHeme2 <- addCellColData(ArchRProj = projHeme2, data = meta_info$predicted.id, cells = projHeme2$cellNames, name = "predictedCelltype")
projHeme2 <- addCellColData(ArchRProj = projHeme2, data = as.character(meta_info$seurat_clusters), cells = projHeme2$cellNames, name = "seuratClusters")
projHeme2 <- addCellColData(ArchRProj = projHeme2, data = meta_info$atac_celltype, cells = projHeme2$cellNames, name = "atac_celltype")
# 找marker基因
markersGS <- getMarkerFeatures(ArchRProj = projHeme2, useMatrix = "GeneIntegrationMatrix", groupBy = "seuratClusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
# 每一个cluster里面的marker基因写到文件里
for (i in seq_along(markerList)) {
  df <- as.data.frame(markerList[[i]])
  write.table(df, paste0("markerList_root_resolut8_", i, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)}

root_marker <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/annotation/root_marker_merged.txt", header=T, sep="\t", stringsAsFactors = F)
root_marker$Validatedmarkergenes <- root_marker$GeneID
data <- data.frame(name = character(), seqnames = character(), start = integer(), end = integer(), strand = integer(), idx = integer(), Log2FC = numeric(), FDR = numeric(), MeanDiff = numeric(), cluster = numeric(), GeneName = character(), MarkerType = character(), References = character(), GeneSymbol = character(), CellType = character(), Organ = character(), Notes = character(), Validatedmarkergenes = character(), stringsAsFactors = FALSE)

for (i in seq_along(markerList)) {
  if(dim(markerList[[i]])[1] != 0){
    df <- as.data.frame(markerList[[i]])
    df$cluster <- i-1
    merged_data <- merge(df, root_marker, by.x = "name", by.y = "GeneID", all.x = TRUE)
    data <- rbind(data,merged_data)}}

write.table(data, "root-0316_markerGene_solut8.txt", sep="\t",quote=F,row.names = F)

# 下面做一些画图验证
marker <- read.table("root-0316_markerGene_solut8.gene", sep="\t",header=T,stringsAsFactors=F)
markerGenes  <- marker[,1]

p3 <- plotGroups(ArchRProj = projHeme2, groupBy = "predictedCelltype", colorBy = "GeneIntegrationMatrix", name = markerGenes, plotAs = "violin", title = "Marker GeneIntegrationMatrix by Clusters", title = tile)
p4 <- plotGroups(ArchRProj = projHeme2, groupBy = "predictedCelltype", colorBy = "GeneScoreMatrix", name = markerGenes, plotAs = "violin", title = tile)

# 把p3和p4合并起来画图
pdf("combined_marker_8_plots_all.pdf", width = 10, height = 8)
# 计算需要多少行
n_plots <- 1005
plots_per_row <- 2
n_rows <- ceiling(n_plots / plots_per_row)
# 循环绘制每一行的图形
for (i in 1:n_rows) {
  start_idx <- (i - 1) * plots_per_row + 1
  end_idx <- min(i * plots_per_row, n_plots)
  plots_row <- c(p3[start_idx:end_idx], p4[start_idx:end_idx])
  do.call(grid.arrange, c(plots_row, ncol = plots_per_row))}
dev.off()

a <- read.table("root-0316_markerGene_solut8.txt", sep="\t",header=T,stringsAsFactors=F)
cM <- as.matrix(confusionMatrix(projHeme2$seuratClusters, projHeme2$atac_celltype))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
b <- as.data.frame(cbind(preClust, rownames(cM)))
colnames(b) <- c("celltype","cluster")
b$cluster <- as.integer(b$cluster)
# unique(unique(projHeme2$predictedGroup))
b_new <- b %>% left_join(a, by = "cluster") 

write.table(b_new, "root-0316_markerGene_solut8.gene.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# markerList$C21 many features
markerGenes  <- c("AT1G31050","AT2G02130","AT1G05470","AT5G37660")
p1 <- plotGroups(ArchRProj = projHeme2, groupBy = "Clusters", colorBy = "GeneScoreMatrix", name = markerGenes, title = "Marker GeneScoreMatrix by Clusters")
dev.off()
pdf("Plot-UMAPHarmony-Marker-Genes-W-Imputation-C2.pdf",width = 10,height = 20)
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],nrow=4)
dev.off()
pdf("Plot-UMAPHarmony-Marker-Genes-W-Imputation-C21.pdf",width = 10,height = 10)
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],nrow=2)
dev.off()

heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", labelMarkers = markerGenes, transpose = TRUE)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-C2", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)
dev.off()

# 在图中显示
projHeme2 <- addImputeWeights(projHeme2)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAPHarmony", imputeWeights = getImputeWeights(projHeme2))
pdf("Plot-UMAPHarmony-Marker-Genes-W-Imputation-C2.pdf",width = 35,height = 35)
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[9]],p1[[10]],p1[[11]],p1[[12]],p1[[13]],p1[[14]],p1[[15]],p1[[16]],p1[[17]],p1[[18]],p1[[19]],p1[[20]],p1[[21]],p1[[22]],p1[[23]],p1[[24]],p1[[25]],nrow=5)
dev.off()

p3 <- plotBrowserTrack(ArchRProj = projHeme2, groupBy = "Clusters", geneSymbol = markerGenes, upstream = 30000, downstream = 30000)
pdf("Plot-Tracks-Marker-Genes.pdf",width = 55,height = 55)
grid.arrange(p3[[1]],p3[[2]],p3[[3]],p3[[4]],p3[[5]],p3[[6]],p3[[7]],p3[[8]],p3[[9]],p3[[10]],p3[[11]],p3[[12]],p3[[13]],p3[[14]],p3[[15]],p3[[16]],p3[[17]],p3[[18]],p3[[19]],p3[[20]],p3[[21]],p3[[22]],p3[[23]],p3[[24]],p3[[25]],nrow=5)
dev.off()

#====================================================================================================================#
#                                       使用ATAC数据进行cluster的注释                                             ####
#====================================================================================================================#
projHeme2 <- readRDS("PCA_snATAC-Integration-root-0316-1.2.rds")
inter <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/annotation/PCA_snATAC_root_combined.rds")

DefaultAssay(inter) <- "ATAC"
all.markers  <- FindAllMarkers(inter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
a <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/annotation/root_marker_merged.txt", sep="\t", header=T)
marker <- a[,1]
tit <- paste(a[,1],a[,6],sep="-")
pdf("markergene-root1.pdf",width = 55,height = 55)
VlnPlot(inter, features = marker)
dev.off()
# 每次画20个基因
markers_split <- split(marker, ceiling(seq_along(marker)/20))
# 循环绘图
plots_V <- lapply(markers_split, function(gene_set) {
  VlnPlot(inter, features = gene_set, ncol = 4) +
    ggtitle(tit) + # ncol可以调横向排版数量
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  #FeaturePlot(inter, features = gene_set, ncol = 4)
})

inter_sub <- subset(inter, idents = c(0,1, 3,4, 5,8,9,10,11:19,21,23,25,26:30,32,34,38,40,41,43:46,48,50,51,53,54,56,57,58,60,62,64:68))


plots_V <- mapply(function(gene, title) {
  VlnPlot(inter_sub, features = gene) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
}, marker, tit, SIMPLIFY = FALSE)

for (i in seq_along(plots_V)) {
  ggsave(filename = paste0("VlnPlot_", i, "_", marker[i], ".pdf"),
         plot = plots_V[[i]],
         width = 16, height = 5)
}

# 依次保存
for (i in seq_along(plots)) {
  ggsave(filename = paste0("VlinPlot_batch_", i, ".pdf"), plot = plots_V[[i]], limitsize = FALSE,width = 100, height = 60)
}

pdf("markergene-root2.pdf",width = 55,height = 55)
FeaturePlot(inter, features = marker)
dev.off()
DotPlot(scedata,features = markers)+coord_flip()
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
pdf("markergene-root3.pdf",width = 55,height = 55)
DoHeatmap(inter, features = top10$gene) + NoLegend()
dev.off()
#细胞群类重命名
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#====================================================================================================================#
#                                        根据ATAC开放手动注释cluster                                              ####
#====================================================================================================================#
inter <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/annotation/PCA_snATAC_root_combined.rds")
cM <- as.matrix(confusionMatrix(inter$seurat_clusters, inter$predicted.id))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
b <- as.data.frame(cbind(preClust, rownames(cM)))
new.cluster.ids <- c("0"="Unknown", "1"="Unknown", "2"="Endodermis", "3"="Vasculature", "4"="Vasculature", "5"="Unknown", "6"="Cortex", "7"="Cortex", "8"="Unknown", "9"="Unknown", "10"="Unknown", "11"="Cortex", "12"="Pericycle", "13"="Unknown", "14"="Vasculature", "15"="Unknown", "16"="Unknown", "17"="Unknown", "18"="Unknown", "19"="Unknown", "20"="Endodermis", "21"="Unknown",  "22"="Trichoblast", "23"="Unknown", "24"="Unknown", "25"="Unknown", "26"="Pericycle", "27"="Endodermis", "28"="Unknown", "29"="Pericycle", "30"="Unknown", "31"="Phloem", "32"="Unknown", "33"="Endodermis", "34"="Unknown", "35"="Atrichoblast", "36"="Atrichoblast", "37"="Endodermis", "38"="Vasculature", "39"="Endodermis", "40"="Unknown", "41"="Unknown",  "42"="Endodermis", "43"="Pericycle", "44"="Trichoblast", "45"="Root Cap,Lateral root cap", "46"="Unknown", "47"="Cortex", "48"="Endodermis", "49"="Atrichoblast", "50"="Vasculature", "51"="Unknown", "52"="Cortex", "53"="Pericycle", "54"="Cortex", "55"="Atrichoblast", "56"="Pericycle", "57"="Xylem  Pole Pericycle", "58"="Unknown", "59"="Cortex", "60"="Unknown", "61"="Phloem",  "62"="Unknown", "63"="Endodermis", "64"="Unknown", "65"="Unknown", "66"="Xylem", "67"="Root Cap,Lateral root cap", "68"="Endodermis")

inter <- RenameIdents(inter, new.cluster.ids)                        
inter$atac_manual_celltype <- inter@active.ident

#画桑葚图
sankey_data <- data.frame(ATAC = inter$atac_celltype, RNA = inter$predicted.id)
sankey_data <- sankey_data[sankey_data$ATAC != "Unknown" & sankey_data$RNA != "Unknown", ]
sankey_data <- as.data.frame(table(sankey_data))
pdf("RNA-ATAC-integration-test12.pdf",width = 10, height = 6)
ggplot(sankey_data, aes(axis1 = ATAC, axis2 = RNA, y = Freq)) + geom_alluvium(aes(fill = ATAC), width = 1/12) + geom_stratum(width = 1/12, color = "grey") + geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + scale_x_discrete(limits = c("ATAC", "RNA"), expand = c(.05, .05)) + theme_minimal() + ggtitle("ATAC to RNA cell type mapping (Sankey Diagram)")
dev.off()

# 相关性可视化
predictions <- table(pbmc.atac$seurat_clusters, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)
predictions <- as.data.frame(predictions)

pdf("RNA-ATAC-integration-test6-lsi8.pdf",width = 10, height = 6)
ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# 整合ATAC注释的cluster的类型
predictions <- table(pbmc.atac$atac_manual_celltype, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

common_levels <- unique(predictions$Var1)  # 或者 unique(predictions$Var2)，只要是一致顺序即可
predictions$Var1 <- factor(predictions$Var1, levels = common_levels)
predictions$Var2 <- factor(predictions$Var2, levels = common_levels)
pdf("RNA-ATAC-integration-test13-match.pdf",width = 10, height = 6)
ggplot(predictions, aes(Var2, Var1, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells", low = "#ffffc8", high = "#7d0025") + xlab("Predicted cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# 提取WOX5基因的开放
inter[["ATAC"]]@data["AT3G11260", ]
# 提取基因的表达
inter[["RNA"]]@data["AT3G11260", ]
# 提取QC细胞
qc_cells <- colnames(inter)[inter$predicted.id == "Quiescent center (QC)"]
# 提取在这洗细胞中的A基因的表达量
expr_values <- inter[["RNA"]]@data["AT3G11260", qc_cells]

# 提取基因表达量在前0.999%的细胞以及查看meta信息
gene_expr <- inter[["ATAC"]]@data["AT3G11260", ]
threshold <- quantile(gene_expr, 0.999) 
cell1 <- names(gene_expr[gene_expr >= threshold])
inter@meta.data[cell1, c("seurat_clusters","predicted.id","atac_celltype") ]

# 提取WOX5基因的开放
inter[["ATAC"]]@data["AT2G18480", ]
# 提取基因的表达
inter[["RNA"]]@data["AT2G18480", ]
# 提取QC细胞
qc_cells <- colnames(inter)[inter$predicted.id == "Quiescent center (QC)"]
# 提取在这洗细胞中的A基因的表达量
expr_values <- inter[["RNA"]]@data["AT2G18480", qc_cells]

# 提取基因表达量在前0.999%的细胞以及查看meta信息
gene_expr <- inter[["ATAC"]]@data["AT2G18480", ]
threshold <- quantile(gene_expr, 0.99) 
cell2 <- names(gene_expr[gene_expr >= threshold])
inter@meta.data[cell2,  c("seurat_clusters","predicted.id","atac_celltype")]
#================================================================================================================####

setwd("/Users/yafeiguo/Documents/1_arabidopsisDevelopment/6_seed/5.development")
data <- read.table("predictedScore.txt",sep="\t", header=T,stringsAsFactors = F)

