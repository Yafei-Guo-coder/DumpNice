#按照Shoot做一个测试
#原始文件路径：
#我的软链接目录：/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/arabidopsis/data/20250225.At_PCA/0.data/snATAC/fragments/Shoot/
#yangquanying的目录：/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangquanyin/shoot/A221008001_Shoot_AT_1008_process.rds
#A221008001_Shoot_AT_1008.fragments.tsv.gz  A221008001_Shoot_AT_1008.fragments.tsv.gz.tbi  A221008002_Shoot_AT_1008.fragments.tsv.gz  A221008002_Shoot_AT_1008.fragments.tsv.gz.tbi
#conda activate /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/conda/env/yafei
#输出文件目录：/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/

Silique/A220221001_Silique_AT_0221.fragments.tsv.gz
Silique/A220322001_Silique_AT_0322.fragments.tsv.gz
Silique/A220322002_Silique_AT_0322.fragments.tsv.gz

############################################################################################# 初步call peaks, 单样本处理 ##############################################################################
library(Signac)
library(Seurat)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(GenomeInfoDb)
library(future)
library(patchwork)
library(IRanges)

options(future.globals.maxSize = 200000 * 1024^2) 
# 读取 GTF 文件并转换为 GRanges 对象 ----

annotations <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangquanyin/Silique/gtf.rds")
seqlevelsStyle(annotations) <- 'UCSC'

seqlevels(annotations) <- sub("^Chr", "chr", seqlevels(annotations))
#处理注释

annotations <- annotations[!is.na(mcols(annotations)$gene_biotype)]
annotations <- annotations[annotations$gene_biotype == "protein_coding"]

annotations$type <- annotations$gene_biotype
Annotation(combined) <- annotations
pbmc <- TSSEnrichment(combined)

#处理片段shoot文件一
frags <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/arabidopsis/data/20250225.At_PCA/0.data/snATAC/fragments/Shoot/A221008001_Shoot_AT_1008.fragments.tsv.gz"
fragment.counts <- CountFragments(frags)
cells.use <- fragment.counts[fragment.counts$frequency_count > 1000, "CB"]
fragments_root1 <- CreateFragmentObject(
  path = frags,
  cells = cells.use,
  validate.fragments = FALSE
)
peaks_root1 <- CallPeaks(fragments_root1, macs2.path = "/jdfsbjcas1/ST_BJ/PUB/Tool/bin/macs2")
counts_root1 <- FeatureMatrix(
  fragments = fragments_root1,
  features = peaks_root1,
  cells = cells.use
)
#创建Seurat对象
atac <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = counts_root1,
    fragments = fragments_root1,
    annotation = annotations
  ),
  assay = "ATAC"
)

saveRDS(object = atac, file = "A221008001_Shoot_AT_1008.fragments_process.rds")
peaks <- granges(atac)
peaks <- as.data.frame(peaks)
write.table(x = peaks, file = "A221008001_Shoot_AT_1008.fragments_peaks.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#处理片段shoot文件二
frags <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/arabidopsis/data/20250225.At_PCA/0.data/snATAC/fragments/Shoot/A221008002_Shoot_AT_1008.fragments.tsv.gz"
fragment.counts <- CountFragments(frags)
cells.use <- fragment.counts[fragment.counts$frequency_count > 1000, "CB"]

fragments_root2 <- CreateFragmentObject(
  path = frags,
  cells = cells.use,
  validate.fragments = FALSE
)

peaks_root2 <- CallPeaks(fragments_root2, macs2.path = "/jdfsbjcas1/ST_BJ/PUB/Tool/bin/macs2")
counts_root2 <- FeatureMatrix(
  fragments = fragments_root2,
  features = peaks_root2,
  cells = cells.use
)

#创建Seurat对象
atac <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = counts_root2,
    fragments = fragments_root2,
    annotation = annotations
  ),
  assay = "ATAC"
)
saveRDS(object = atac, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/A221008002_Shoot_AT_1008.fragments_process.rds")
peaks <- granges(atac)
peaks <- as.data.frame(peaks)
write.table(x = peaks, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/A221008002_Shoot_AT_1008.fragments_peaks.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

######################################################################################### 整合样本peak ############################################################################################### 

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

# read in peak sets
#如果只有fragment文件，那么需要单独跑一边frag调用的macs2的峰输出后再读入峰文件
peaks_A220711006_7d_Root_AT_0711 <- read.table(
  file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/A221008001_Shoot_AT_1008.fragments_peaks.bed",
  col.names = c("chr", "start", "end", "width",	"strand"), header=T
)
peaks_A220711005_7d_Root_AT_0711 <- read.table(
  file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/A221008002_Shoot_AT_1008.fragments_peaks.bed",
  col.names = c("chr", "start", "end",  "width",	"strand"),header=T
)

# convert to genomic r
p.220711006 <- makeGRangesFromDataFrame(peaks_A220711006_7d_Root_AT_0711)
p.220711005 <- makeGRangesFromDataFrame(peaks_A220711005_7d_Root_AT_0711)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(p.220711006, p.220711005))
# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 8000 & peakwidths > 20]
combined.peaks
#读取碎片数据
frags_220711006 <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/arabidopsis/data/20250225.At_PCA/0.data/snATAC/fragments/Shoot/A221008001_Shoot_AT_1008.fragments.tsv.gz"
frags_220711005 <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/arabidopsis/data/20250225.At_PCA/0.data/snATAC/fragments/Shoot/A221008002_Shoot_AT_1008.fragments.tsv.gz"

annot <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangquanyin/Silique/gtf.rds")
seqlevelsStyle(annot) <- 'UCSC'

fragment220711006.counts <- CountFragments(frags_220711006)
fragment220711005.counts <- CountFragments(frags_220711005)

cells.use220711006 <- fragment220711006.counts[fragment220711006.counts$frequency_count > 1000, "CB"]
cells.use220711005 <- fragment220711005.counts[fragment220711005.counts$frequency_count > 1000, "CB"]

fragments220711006 <- CreateFragmentObject(
  path = frags_220711006,
  cells = cells.use220711006,
  validate.fragments = FALSE
)

fragments220711005 <- CreateFragmentObject(
  path = frags_220711005,
  cells = cells.use220711005,
  validate.fragments = FALSE
)

fragments220711006_root.counts <- FeatureMatrix(
  fragments = fragments220711006,
  features = combined.peaks,
  cells = cells.use220711006
)
# 提取 GRanges 对象的 seqnames

# 第三步：查看更新后的 FeatureMatrix 对象
#head(seqnames(A220322002_Silique_AT_0322.counts))

fragments220711005_root.counts <- FeatureMatrix(
  fragments = fragments220711005,
  features = combined.peaks,
  cells = cells.use220711005
)

#创建三个assay
fragments220711006_root_assay <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = fragments220711006_root.counts,
    fragments = fragments220711006,
    annotation = annot
  ),
  assay = "ATAC"
)

fragments220711005_root_assay <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = fragments220711005_root.counts,
    fragments = fragments220711005,
    annotation = annot
  ),
  assay = "ATAC"
)

# add information to identify dataset of origin
fragments220711006_root_assay$dataset <- 'A221008001_Shoot'
fragments220711005_root_assay$dataset <- 'A221008002_Shoot'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = fragments220711006_root_assay,
  y = fragments220711005_root_assay,
  add.cell.ids = c("A221008001_Shoot","A221008002_Shoot")
)
print(combined[["ATAC"]])
print(head(combined[["ATAC"]]))

saveRDS(object = combined, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/Shoot_AT_1008.fragments_combined.rds")

################################################################################################### QC ################################################################################################


combined <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/Shoot_AT_1008.fragments_combined.rds")

annotations <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangquanyin/Silique/gtf.rds")
seqlevelsStyle(annotations) <- 'UCSC'
seqlevels(annotations) <- sub("^chr", "Chr", seqlevels(annotations))
annotations <- annotations[!is.na(mcols(annotations)$gene_biotype)]
annotations <- annotations[annotations$gene_biotype == "protein_coding"]
annotations$type <- annotations$gene_biotype
Annotation(combined) <- annotations
#QC
# compute nucleosome signal score per cell并可视化
pbmc <- NucleosomeSignal(object = combined)

#在这里卡了好久的bug，原因可能是染色体号不对应，fragment里面是Chr，但是annotations里面是chr，或者是需要做上述的注释部分的处理，这一步非常地耗时，跑完保存
#这一步还是没有跑出来，跑着跑着就中断了，不知道什么原因


pbmc <- TSSEnrichment(object = pbmc, fast = TRUE)

# 当fast=T时，仅计算TSS分数，而不存储每个细胞Tn5插入频率的位置矩阵，这样后续将不能使用TSSPlot画图。
saveRDS(object = pbmc, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/Shoot_AT_1008.fragments_combined_QC.rds")



#pbmc <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/Signac/Shoot_AT_1008.fragments_combined_QC.rds")

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 0.5, 'NS > 0.5', 'NS < 0.5') #这里的阈值本来是4
figNS <- FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
# 基于TSS分数将不同的细胞进行分组以及画所有TSS位点的可及性信号图来检查TSS富集分数
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 1, 'High', 'Low')

print("step1 ok")
#figTSS <- TSSPlot(atac, group.by = 'high.tss')
#figTSS_group <- TSSPlot(atac, group.by = 'dataset') 

DS <- DensityScatter(
  pbmc,
  x = 'nCount_ATAC',
  y = 'TSS.enrichment',
  log_x = TRUE,
  quantiles = TRUE
)

# 修改颜色为蓝色渐变
DS <- DS + scale_fill_gradientn(
  colors = c("white", "lightblue", "blue", "darkblue"),
  name = "Density"  # 图例名称
)

# 添加标题和轴标签
DS <- DS + ggtitle("Density Scatter Plot: nCount_ATAC vs. TSS Enrichment") +
  xlab("Log-scaled nCount_ATAC") +
  ylab("TSS Enrichment")

# 美化图形主题
DS <- DS + theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )


QC <- VlnPlot(
  object = pbmc,
  features = c( 'nCount_ATAC','nucleosome_signal','TSS.enrichment'),
  pt.size = 0,
  ncol = 4
) 


ggsave(filename = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/shoot_NS.a1.png", plot = figNS, height = 6, width =12, dpi = 300)
#ggsave(filename = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/shoot_TSS.a1.png", plot = figTSS, height = 6, width = 12, dpi = 300)
ggsave(filename = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/shoot_QC.a1.png", plot = QC, height = 6, width = 12, dpi = 300)
ggsave(filename = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/shoot_DS.a1.png", plot = DS, height = 6, width = 12, dpi = 300)
#ggsave(filename = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/shoot_figTSS_group.a1.png", plot = figTSS_group, height = 6, width = 12, dpi = 300)
print("step2 ok")


#QC质量筛选，根据数据QC分布（小提琴图）设置数值
atac <- subset(
  x = atac,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 
)

#peaks <- granges(atac)
#peaks <- as.data.frame(peaks)
#write.table(x = peaks, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangquanyin/Silique/atac/process_data/silique_AT_combined_peaks.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#writeLines(text = colnames(x = atac), con = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangquanyin/Silique/atac/process_data/silique_AT_combined_cell.txt")

# cluster and make UMAP
atac <- FindTopFeatures(atac, min.cutoff = 10)
atac <- RunTFIDF(atac)
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30,min.dist = 0.1)
atac <- FindNeighbors(atac, reduction = "lsi", dims = 2:30)
atac <- FindClusters(atac, algorithm = 1, resolution = 1.2)
UMap <- DimPlot(object = atac, label = TRUE) + NoLegend()


#存基础分析RDS数据，后面可以与RNA数据联合分析

saveRDS(object = atac, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/Shoot_AT_1008.fragments_combined_QC_UMAP.rds")

print("step3 ok")
#图像
ggsave(filename = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/shoot_UMap.png", plot = UMap, height = 6, width = 12, dpi = 300)

# Load the pre-processed scRNA-seq data for ATs
atac_rna <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202407.At_PCA/2.Cluster2024/snRNA_202409/Shoot.rds")
atac_rna_Umap <- DimPlot(atac_rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()  #查看scRNAUmap结果并画图
print(atac_rna_Umap)

atac_rna <- UpdateSeuratObject(atac_rna)#将较旧版本的 Seurat 对象更新到当前版本所支持的格式，可选
atac_rna <- RenameIdents(
  object = atac_rna,
  '0' = 'unknown',
  '1' = 'Cortex',
  '2' = 'Trichoblast type I',
  '3' = 'Endodermis(QC)',
  '4' = 'Stressed cells',
  '5' = 'epidermis',
  '6' = 'Pericycle and vascular ( not xylem)', 
  '7' = 'vascular bundles',
  '8' = 'Cortex',
  '9' = 'unknown',
  '10' = 'Atrichoblast',
  '11' = 'Atrichoblast',
  '12' = 'Trichoblast type II (atrichoblast feature）',
  '13' = 'xylem(meristem feature)',
  '14' = 'phloem(meristem feature)',
  '15' = 'Atrichoblast',
  '16' = 'Dividing cells',
  '17' = 'Low quality cells',
  '18' = 'Endodermis',
  '19' = 'Cortex (QC)',
  '20' = 'Trichoblast type I',
  '21' = 'Endodermis',
  '22' = 'cortex',
  '23' = 'Endodermis',
  '24' = 'unknown',
  '25' = 'Root Cap',
  '26' = 'Trichoblast type I',
  '27' = 'Xylem  Pole Pericycle(meristem feature)',
  '28' = 'Phloem ',
  '29' = 'QC',
  '30' = 'Xylem type(death) III',
  '31' = 'Cortex (分化)'
)


atac_rna$celltype <- Idents(atac_rna)
#给RNA矩阵加上标签
atac_rna$tech <- 'scRNA'
DefaultAssay(atac_rna) <- 'RNA'
atac_rna_Umap <- DimPlot(atac_rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()  #查看scRNAUmap结果并画图
print(atac_rna_Umap)
#读入ATAC数据
#atac <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangquanyin/Silique/atac/process_data/silique_AT_combined_pre_analyze.rds")
#给ATAC加标签
atac$tech <- 'scATAC'


#画迁移前的图像
pRNA <- DimPlot(atac_rna, label=TRUE, repel=TRUE) + NoLegend() + ggtitle("SnATAC")
pATAC <- DimPlot(atac, reduction = "umap", label = TRUE)+ NoLegend() + ggtitle("snRNA")
#print(pRNA + pATAC)
# 检查参考数据中的特征
#print(rownames(atac_rna))

# 检查查询数据中的特征
#print(rownames(atac))
#common_features <- intersect(rownames(atac_rna), rownames(atac))
#print(common_features)
gene.activities <- GeneActivity(atac)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(
  object = atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_RNA)
)
#在标签迁移前必须做这一步，否则会出现没有共同特征的报错
DefaultAssay(atac) <- 'RNA'

#标签迁移
#通过找到参考数据和查询数据之间的对应关系（锚点），并使用这些锚点来预测查询数据的细胞类型。

#查询锚点
transfer.anchors <- FindTransferAnchors(
  reference = atac_rna,
  query = atac,
  reduction = 'cca'
)

#根据参考数据（atac_rna）中的细胞类型标签，推测查询数据中的细胞类型。
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = atac_rna$celltype,
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

atac <- AddMetaData(object = atac, metadata = predicted.labels)

plot1 <- DimPlot(
  object = atac_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

pz <- plot1 + plot2
ggsave(filename = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/shoot_anchor.a1.png", plot = pz, height = 6, width = 12, dpi = 300)
#绘制预测分数小提琴图，较高的预测分数通常表示细胞类型预测的可信度较高
vln <- VlnPlot(atac, 'prediction.score.max', group.by = 'predicted.id')
ggsave(filename = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/shoot_vln.score.a1.png", plot = vln, height = 6, width = 12, dpi = 300)

saveRDS(object = atac, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/Shoot_AT_1008.fragments_combined_QC_UMAP_anchor.rds")

######################################################################################################################################################################################################
# 添加黑名单比例和片段在峰值区域的比例
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# 可视化变量之间的关系
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# 添加黑名单比例和片段在峰值区域的比例
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# 可视化变量之间的关系
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

############################################################################################ 可以自己使用ggplot2画图 #################################################################################
#绘制质控信息和TSS富集得分的对比图
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

png("TSS-vs-Frags.png")
plot(p)
dev.off()
getwd()

#为ArchRProject每个样本绘制统计信息
p1 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

#根据TSS富集得分为每个样本绘制小提琴图
p2 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

#根据log10(unique nuclear fragments)为每个样本绘制山脊图
p3 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p3

#根据log10(unique nuclear fragments)为每个样本绘制小提琴图。

p4 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)


#绘制样本的TSS富集谱和Fragment大小分布
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p1

#TSS富集谱:

p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
p2

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

#保存和加载
saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = FALSE)

#使用ArchRProject过滤doublets
projHeme2 <- filterDoublets(projHeme1)
# Filtering 410 cells from ArchRProject!
#   scATAC_BMMC_R1 : 243 of 4932 (4.9%)
#   scATAC_CD34_BMMC_R1 : 107 of 3275 (3.3%)
#   scATAC_PBMC_R1 : 60 of 2454 (2.4%)

projHemeTmp <- filterDoublets(projHeme1, filterRatio = 1.5)
# Filtering 614 cells from ArchRProject!
#   scATAC_BMMC_R1 : 364 of 4932 (7.4%)
#   scATAC_CD34_BMMC_R1 : 160 of 3275 (4.9%)
#   scATAC_PBMC_R1 : 90 of 2454 (3.7%)





