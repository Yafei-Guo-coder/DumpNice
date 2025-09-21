library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(multtest)
library(gridExtra)
library(ComplexHeatmap)
library(ggalt)
library(clustree)
#####以Flower为例
pbmc.data <- Read10X(data.dir = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/arabidopsis/data/20250225.At_PCA/0.data/snRNA/Flower/R1063_Flower_NPSB_AT_0309/04.Matrix", gene.column = 1)
#####自动读取10X的数据，是一些tsv与mtx文件
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "biomamba")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200) #问题一：实际处理中这两个参数是怎么设置的？
##pbmc
##An object of class Seurat 
##34205 features across 6926 samples within 1 assay 
##Active assay: RNA (34205 features, 0 variable features)
##pbmc$
##pbmc$orig.ident    pbmc$nCount_RNA    pbmc$nFeature_RNA  

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-") #问题二：拟南芥中匹配线粒体基因的表达式是什么？
p <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("violin_plot_raw.pdf", plot = p, width = 6, height = 2)

p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave("combined_plot.png", combined_plot, width = 10, height = 5)

#过滤具有独特特征计数超过2500或少于200的细胞并且过滤线粒体计数为0.5%的细胞
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
p <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("violin_plot_filter.pdf", plot = p, width = 6, height = 2)

p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave("combined_plot_filter.png", combined_plot, width = 10, height = 5)

#对每一个细胞的基因表达进行标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#鉴定在细胞之间表达高变的基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #这里默认返回2000个高变feature

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
p1 <- VariableFeaturePlot(pbmc)
p2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave("combined_plot_highvari.png", combined_plot, width = 10, height = 5)

#使用线性变换（缩放），这是在PCA等降维技术之前的标准预处理步骤。ScaleData（）函数：调整每个基因的表达，使细胞间的平均表达为0缩放每个基因的表达，使细胞间的方差为1。
#这一步在下游分析中给予相同的权重，因此高表达的基因不占主导地位。
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

#根据基因表达进行降维
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#VizDimReduction(), DimPlot(), and DimHeatmap()
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

#PCA降维
p <- DimPlot(pbmc, reduction = "pca") + NoLegend()
ggsave("pca_plot.png", p, width = 10, height = 10)

p <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") #展示前三十个基因的loading
ggsave("pcaDimLoading_plot.png", p, width = 10, height = 10)

p <- DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) #问题三：这个热图一直画不出来
ggsave("pcaHeatmap_plot.png", p, width = 50, height = 50, limitsize=FALSE)

pdf("pcaheat_plot.pdf") #换了不同的画图方式还是画不出来
draw(p)
dev.off()

p <- DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) #这个也是一样的画不出来
pdf("pca15heat_plot.pdf")
print(p)
dev.off()

p <- ElbowPlot(pbmc) #类似于碎石图
ggsave("ElbowPlot.png", p, width = 10, height = 10)

#接下来是根据PC进行细胞聚类，上述过程可以选择出合适的PC数量，可以倾向于选择多个PC以便于囊括更多的信息
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#根据基因表达进行UMAP降维
pbmc <- RunUMAP(pbmc, dims = 1:10) #问题四：什么决定了聚类的数量，为什么有的是11，有的是22
p <- DimPlot(pbmc, reduction = "umap")
ggsave("umap_plot2.png", p, width = 10, height = 10)


#鉴定差异表达特征（聚类生物标志物）
#找cluster2中的所有marker
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

#找出区分聚类5与聚类0和聚类3的所有标记
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

#与所有剩余的细胞相比，找到每个集群的标记，只报告阳性的
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#可视化，建议探索ridgeplot(), cellScatter()和dotplot()作为查看数据集的其他方法
p <- VlnPlot(pbmc, features = c("AT3G06355", "AT1G35310"))
ggsave("exprediff_plot.png", p, width = 10, height = 5)

#画原始数量分布
p <- VlnPlot(pbmc, features = c("AT3G06355", "AT1G35310"), slot = "counts", log = TRUE)
ggsave("exprediff_raw_plot.png", p, width = 10, height = 5)

#画不同的feature在细胞cluster中的分布
p <- FeaturePlot(pbmc, features = c("AT3G06355", "AT1G35310", "AT3G41768", "AT2G47040"))
ggsave("feature_plot.png", p, width = 10, height = 10)

#DoHeatmap()为给定的细胞和特征生成一个表达式热图。在本例中，我们为每个集群绘制前20个标记（如果少于20个，则绘制所有标记）。
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
p <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave("heatmap_plot.png", p, width = 10, height = 10)

#####仅有一个稀疏矩阵时的读取方法-----
matrix_data <- read.table("single_cell_datamatrix.txt", sep="\t", header=T, row.names=1)
dim(matrix_data)
seurat_obj <- CreateSeuratObject(counts = matrix_data)

#####读取RDS文件-----
rm(list = ls())
pbmc <- readRDS("panc8.rds")
saveRDS(pbmc,"pbmc.rds")

#可视化聚类图-----
seurat_object <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/arabidopsis/data/20250225.At_PCA/2.Cluster2024/snRNA_202502.summary/rds/Flower_final202502.rds")
pbmc <- seurat_object

p <- DimPlot(pbmc, reduction = "umap", group.by = "cell_class", label=TRUE)
ggsave("cell_class.pdf", plot = p, width = 10, height = 8)

pdf()

p <- DimPlot(pbmc, reduction = "umap", group.by = "cell_type", label=TRUE)
ggsave("cell_type.pdf", plot = p, width = 10, height = 8)

### Cluster number
# Ruiying Chen


## Flower ----------------------------------------------------------------------
data <- readRDS("00.Data/Resolu.rds/Flower_scRNA.harmony.resolu.rds")
### Plot clustree
## With filter probability
clus.tree.out <- clustree(data, prop_filter = 0.01) +
  theme(legend.position = "right") + scale_edge_color_continuous(low = "#D42825", high = "#D42825")
pdf("00.Data/Resolu.rds/00.Plots.fix/Flower_cluster_tree.fix.pdf",
    width = 12, height = 20)
print(clus.tree.out)
dev.off()

