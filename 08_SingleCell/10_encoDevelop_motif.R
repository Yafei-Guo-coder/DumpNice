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
library(universalmotif, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
library(tidydr)
library(ggalluvial)
library(ggplot2)
library(chromVAR, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
library(TFBSTools, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
library(SummarizedExperiment)
library(motifmatchr)
library(Matrix)
library(rhdf5)
library(readxl)
library(ggrepel)
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
library(GenomeInfoDb)

projAll <- readRDS("rds_Endosperm/ArchR_EndospermD_rmTorpedo_reAnchor2_motif2.rds")
projAll <- addImputeWeights(projAll)

#====================================================================================================================#
###                  plot trajectory heatmap - use new motif matrix based on pseudotime trajectory                 ###
#====================================================================================================================#

projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage", threads = 1,force = TRUE)
pathToMacs2 <- findMacs2()
# 使用w/MACS2鉴定peak
projAll <- addReproduciblePeakSet(ArchRProj = projAll, groupBy = "Stage", pathToMacs2 = pathToMacs2, genomeSize = "1.25e+08")
projAll <- addPeakMatrix(projAll)
# 添加motif
meme_file <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/Ath_TF_binding_motifs_add.meme" 
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

projAll <- addMotifAnnotations(ArchRProj = projAll, motifPWMs = pwmlist, name = "Motif", force = TRUE)

# Motif偏离，不基于细胞分类情况
projAll <- addBgdPeaks(projAll,force = TRUE)
# 根据所有的motif注释计算每个细胞的偏离值
projAll <- addDeviationsMatrix(ArchRProj = projAll, peakAnnotation = "Motif", force = TRUE)

# 保存
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor2_motif3.rds")

gene <- c("AT3G26744","AT1G71250","AT1G21970","AT1G08560","AT1G30690","AT1G48930","AT3G62830","AT5G13170","AT3G08900","AT1G68560","AT2G47650","AT2G25540","AT4G33600","AT3G19820","AT4G04460","AT1G62340","AT5G50750","AT3G48740","AT5G50260","AT2G19500","AT1G23320","AT2G35670","AT5G60440","AT2G35230","AT1G65330","AT1G67775","AT1G50650","AT1G62340","AT3G08900","AT2G47650","AT1G49770","AT1G68560","AT3G11520","AT1G48930","AT1G75080","AT5G07280","AT1G71250","AT3G48740","AT1G11190","AT5G18270")

# plot motifMatrix
trajMM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, varCutOff = 0,pal = paletteContinuous(set = "horizonExtra"),rowOrder = names(trajMM))

# plot geneScoreMatrix 
trajGSM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM, varCutOff = 0,  pal = paletteContinuous(set = "horizonExtra"))

# plot geneIntegrationMatrix 
trajGIM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM, varCutOff = 0, pal = paletteContinuous(set = "horizonExtra"))

# trajPM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
# p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "horizonExtra")) 

# select gene to plot

pattern <- paste(gene, collapse = "|")
motifs <- grep(pattern, rownames(trajGSM), value = TRUE)


# motifs_z <- motifs[grep("^z:", motifs)]
trajGSM_sub <- trajGSM[motifs, ]
trajGIM_sub <- trajGIM[motifs, ]

mat <- t(scale(t(assay(trajGSM_sub))))

p2_2 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)    
mat <- t(scale(t(assay(trajGIM_sub))))      
p3_2 <- pheatmap(mat, color = paletteContinuous(set = "horizonExtra"),
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)       

pdf("EndospermTrajectory_umap_reanchor_add_plot1.pdf", width = 6, height = 8)
print(p1)
print(p2)
print(p2_2)
print(p3)
print(p3_2)
dev.off()

# combine motifMatrix and geneScore matrix

corGSM_MM <- correlateTrajectories(trajGSM, trajMM, corCutOff = -1,  varCutOff1 = 0, varCutOff2 = 0)

idxToRemove <- grep(pattern = "deviations", x = corGSM_MM[["correlatedMappings"]]$name2)
if(length(idxToRemove > 1)){
  corGSM_MM[["correlatedMappings"]] <- corGSM_MM[["correlatedMappings"]][-idxToRemove,]
}

# data <- as.data.frame(corGSM_MM[["correlatedMappings"]])
trajGSM2 <- trajGSM[corGSM_MM[["correlatedMappings"]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[["correlatedMappings"]]$name2, ]

trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)


rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajGSM2, labelTop = 1, pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, labelTop = 1, pal = paletteContinuous(set = "horizonExtra"), varCutOff = 0, rowOrder = rowOrder)

# ComplexHeatmap::draw(ht1 + ht2)
pdf("EndospermTrajectory_umap_reanchor_add_plot2.pdf", width = 6, height = 8)
print(ht1)
print(ht2)
dev.off()

# combine motifMatrix and geneInteration matrix
corGIM_MM <- correlateTrajectories(trajGSM, trajMM,  corCutOff = -1,  varCutOff1 = 0, varCutOff2 = 0)
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
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "horizonExtra"), varCutOff = 0, rowOrder = rowOrder)
# ComplexHeatmap::draw(ht1 + ht2)

pdf("EndospermTrajectory_umap_reanchor_add_plot3.pdf", width = 6, height = 8)
print(ht1)
print(ht2)
dev.off()

#====================================================================================================================#
###                                                  正向TF-调控因子鉴定                                             ###
#====================================================================================================================#

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
  
pdf("Endosperm-Tracks-Marker-Genes-with-Peak2Gene-heatmap2_add_plot5-stage.pdf",width = 8,height = 8)

print(p_GSM)
print(p_GIM)
dev.off()


seGroupMotif <- getGroupSE(ArchRProj = projAll, useMatrix = "MotifMatrix", groupBy = "predictedScore_endo")
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
  
pdf("Endosperm-Tracks-Marker-Genes-with-Peak2Gene-heatmap2_add_plot5-endo.pdf",width = 8,height = 8)

print(p_GSM)
print(p_GIM)
dev.off()


$ projAll <- addImputeWeights(projAll)
p1 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "MotifMatrix", name = "z:AT1G75080", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))
p2 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "GeneScoreMatrix", name = "AT1G75080", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))
p3 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "GeneIntegrationMatrix", name = "AT1G75080", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))
p4 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "MotifMatrix", name = "z:AT3G26744", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))
p5 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "GeneScoreMatrix", name = "AT3G26744", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))
p6 <- plotEmbedding(ArchRProj = projAll, pal = paletteContinuous(set = "solarExtra"), colorBy = "GeneIntegrationMatrix", name = "AT3G26744", embedding = "UMAP",imputeWeights = getImputeWeights(projAll))

pdf("integration_add_plot6.pdf", width = 13, height = 13)
grid.arrange(p1, p2, p3,p4,p5,p6, ncol = 3, top = "combined")
dev.off()

#====================================================================================================================#
#                                        整合分析:peak关联gene调控                                                ####
#====================================================================================================================#
# 可以只用ATAC-seq数据进行分析，如识别peak之间的共开放性来预测调控相互作用，也可以整合scRNA-seq数据，如通过peak-基因的连锁分析预测增性子活性
# 查看peak之间的共开放相关性
projAll <- addCoAccessibility(ArchRProj = projAll, reducedDims = "IterativeLSI")

# markerGenes  <- c("AT1G75080","AT3G26744")
markerGenes  <- c("AT3G19820","AT1G48930","AT3G48740","AT1G68560","AT1G21970","AT1G62340","AT2G47650","AT3G62830","AT3G08900","AT1G08560","AT1G30690","AT5G50750","AT1G71250","AT4G04460","AT2G25540","AT5G13170","AT4G33600","AT5G50260")
p1 <- plotBrowserTrack(ArchRProj = projAll, groupBy = "Stage", geneSymbol = markerGenes, upstream = 10000, downstream = 5000, loops = getCoAccessibility(projAll))
p2 <- plotBrowserTrack(ArchRProj = projAll, groupBy = "predictedGroup_endo", geneSymbol = markerGenes, upstream = 10000, downstream = 5000, loops = getCoAccessibility(projAll))

pdf("Endosperm-Tracks-Marker-Genes-with-CoAccessibility_add_plot7-stage.pdf",width = 35,height = 25)
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[9]],p1[[10]],p1[[11]],p1[[12]],p1[[13]],p1[[14]],p1[[15]],p1[[16]], p1[[17]],p1[[18]], nrow=5)
dev.off()

pdf("Endosperm-Tracks-Marker-Genes-with-CoAccessibility_add_plot7-endo.pdf",width = 35,height = 25)
grid.arrange(p2[[1]],p2[[2]],p2[[3]],p2[[4]],p2[[5]],p2[[6]],p2[[7]],p2[[8]],p2[[9]],p2[[10]],p2[[11]],p2[[12]],p2[[13]],p2[[14]],p2[[15]],p2[[16]], p2[[17]],p2[[18]], nrow=5)
dev.off()

# Peak2GeneLinkages peak关联基因分析
projAll <- addPeak2GeneLinks(ArchRProj = projAll, reducedDims = "IterativeLSI")
# p2g <- getPeak2GeneLinks(ArchRProj = projAll, corCutOff = 0.45, resolution = 1, returnLoops = FALSE)
# metadata(p2g)[[1]]
# p2g <- getPeak2GeneLinks(ArchRProj = projAll, corCutOff = 0.45, resolution = 10000, returnLoops = TRUE)
# 在browser track中绘制每种细胞类型的peak-to-gene连接

p1 <- plotBrowserTrack(ArchRProj = projAll, groupBy = "Stage", geneSymbol = markerGenes, upstream = 10000, downstream = 5000, loops = getPeak2GeneLinks(projAll))
p2 <- plotBrowserTrack(ArchRProj = projAll, groupBy = "predictedGroup_endo", geneSymbol = markerGenes, upstream = 10000, downstream = 5000, loops = getPeak2GeneLinks(projAll))

pdf("Endosperm-Tracks-Marker-Genes-with-CoAccessibility_add_plot8-stage.pdf",width = 35,height = 25)
grid.arrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],p1[[5]],p1[[6]],p1[[7]],p1[[8]],p1[[9]],p1[[10]],p1[[11]],p1[[12]],p1[[13]],p1[[14]],p1[[15]],p1[[16]], p1[[17]],p1[[18]], nrow=5)
dev.off()

pdf("Endosperm-Tracks-Marker-Genes-with-CoAccessibility_add_plot8-endo.pdf",width = 35,height = 25)
grid.arrange(p2[[1]],p2[[2]],p2[[3]],p2[[4]],p2[[5]],p2[[6]],p2[[7]],p2[[8]],p2[[9]],p2[[10]],p2[[11]],p2[[12]],p2[[13]],p2[[14]],p2[[15]],p2[[16]], p2[[17]],p2[[18]], nrow=5)
dev.off()

# 绘制Peak-to-gene连接热图
p1 <- plotPeak2GeneHeatmap(ArchRProj = projAll, groupBy = "Stage")
p2 <- plotPeak2GeneHeatmap(ArchRProj = projAll, groupBy = "predictedGroup_endo")

pdf("Endosperm-Tracks-Marker-Genes-with-CoAccessibility_add_plot9.pdf",width = 8,height = 10)
draw(p1)
draw(p2)
dev.off()

# 保存
saveRDS(projAll, "ArchR_EndospermD_rmTorpedo_reAnchor2_motif3_track.rds")

#====================================================================================================================#
###                                                    motif足迹分析                                               ####
#====================================================================================================================#

motifPositions <- getPositions(projAll)
# 提取部分感兴趣的TF motifs用于展示
# motifs <- c("AT1G75080","AT3G26744")
motifs <- c("AT1G75080")
projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "predictedScore_endo",force = TRUE)
seFoot_endo <- getFootprints(ArchRProj = projAll, positions = motifPositions[motifs], groupBy = "predictedScore_endo")
# 减去Tn5偏好
pdf("Endosperm_footprint_add_plot4-1.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot_endo, plotName = "Footprints-Subtract-Bias-Stage", ArchRProj = projAll, normMethod = "Subtract", addDOC = FALSE, smoothWindow = 5)
dev.off()

projAll <- addGroupCoverages(ArchRProj = projAll, groupBy = "Stage",force = TRUE)
seFoot_stage <- getFootprints(ArchRProj = projAll, positions = motifPositions[motifs], groupBy = "Stage")
# 减去Tn5偏好
pdf("Endosperm_footprint_add_plot4-2.pdf", width = 20, height = 20)
plotFootprints(seFoot = seFoot_stage, plotName = "Footprints-Subtract-Bias-Stage", ArchRProj = projAll, normMethod = "Subtract", addDOC = FALSE, smoothWindow = 5)
dev.off()


#====================================================================================================================#
###                           plot trajectory heatmap - combined GSM GIM and MM matrix                             ###
#====================================================================================================================#

projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor2_motif3.rds")
# plot motifMatrix
trajMM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
# plot geneScoreMatrix 
trajGSM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
# plot geneIntegrationMatrix 
trajGIM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)

# combine motifMatrix and geneInteration matrix
corGIM_MM <- correlateTrajectories(trajGSM, trajMM,  corCutOff = -1,  varCutOff1 = 0, varCutOff2 = 0)
idxToRemove2 <- grep(pattern = "deviations", x = corGIM_MM[["correlatedMappings"]]$name2)
if(length(idxToRemove2 > 1)){
  corGIM_MM[["correlatedMappings"]] <- corGIM_MM[["correlatedMappings"]][-idxToRemove2,]
}

trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajGSM2 <- trajGSM[corGIM_MM[[1]]$name1, ]
data <- as.data.frame(corGIM_MM[1])

# constant oppening peaks
sub1 <- data
corGSM_GIM_C <- correlateTrajectories(trajGSM2, trajGIM2, corCutOff = 0.3, varCutOff1 = 0, varCutOff2 = 0)
sub2 <- as.data.frame(corGSM_GIM_C[[1]])
name_1 <- union(sub1[which(sub1$VarAssay1 < 0.1 | sub1$VarAssay2 < 0.1),]$name1, sub2[which(sub2$VarAssay1 < 0.1 | sub2$VarAssay2 < 0.1),]$name1)


trajGIM_C <- trajGIM2[name_1, ]
mat <- t(scale(t(assay(trajGIM_C))))
mat <- mat[!rowSums(is.na(mat)) > 0, ]

name_2 <- data[match(rownames(mat), data$name1), 8]
trajGIM_C_2 <- trajGIM_C[rownames(mat),]

trajGSM_C <- trajGSM2[rownames(mat), ]
trajMM_C <- trajMM[name_2, ]

trajCombined2_C <- trajGIM_C_2
combinedMat2_C <- plotTrajectoryHeatmap(trajCombined2_C,limits = c(-10, 10), returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat2_C), rownames(trajGIM_C_2))

ht1 <- plotTrajectoryHeatmap(trajGSM_C, labelTop = 50, pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajGIM_C_2, labelTop = 40, pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht3 <- plotTrajectoryHeatmap(trajMM_C, labelTop = 50, pal = paletteContinuous(set = "horizonExtra"), varCutOff = 0, rowOrder = rowOrder)



pdf("combined_constant_V0.1.pdf", width = 20, height = 8)
ComplexHeatmap::draw(ht1 + ht2 + ht3)
dev.off()

# all oppening peaks
p1 <- plotTrajectoryHeatmap(trajMM, varCutOff = 0,limits = c(-10, 10), pal = paletteContinuous(set = "horizonExtra"))
p2 <- plotTrajectoryHeatmap(trajGSM2, varCutOff = 0, limits = c(-10, 10),  pal = paletteContinuous(set = "horizonExtra"))
mat <- t(scale(t(assay(trajGIM2))))
mat <- mat[!rowSums(is.na(mat)) > 0, ]
trajGIM3 <- trajGIM2[rownames(mat),]
p3 <- plotTrajectoryHeatmap(trajGIM3, varCutOff = 0, limits = c(-10, 10),  pal = paletteContinuous(set = "horizonExtra"))
pdf("combined_constant_all.pdf",width = 8,height = 8)
print(p1)
print(p2)
print(p3)
dev.off()

# combine geneScoreMatrix and geneInteration matrix
corGSM_GIM <- correlateTrajectories(trajGSM2, trajGIM2, corCutOff = 0.3, varCutOff1 = 0, varCutOff2 = 0)

# positive regulation
select <- data[which(data$Correlation > 0.3 & (data$name1 %in% corGSM_GIM[[1]]$name1)),]
trajGSM3 <- trajGSM2[select$name1, ]
trajMM3 <- trajMM[select$name2, ]
trajGIM3 <- trajGIM2[select$name1, ]
trajCombined2 <- trajGSM3
combinedMat2 <- plotTrajectoryHeatmap(trajCombined2, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat2), rownames(trajGSM3))
ht1 <- plotTrajectoryHeatmap(trajGSM3, labelTop = 1, pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajGIM3, labelTop = 1, pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht3 <- plotTrajectoryHeatmap(trajMM3, labelTop = 1, pal = paletteContinuous(set = "horizonExtra"), varCutOff = 0, rowOrder = rowOrder)
pdf("combined_R0.5.pdf", width = 20, height = 8)
ComplexHeatmap::draw(ht1 + ht2 + ht3)
dev.off()

# negative regulation
select <- data[which(data$Correlation < (-0.3) & (data$name1 %in% corGSM_GIM[[1]]$name1)),]
trajGSM3 <- trajGSM2[select$name1, ]
trajMM3 <- trajMM[select$name2, ]
trajGIM3 <- trajGIM2[select$name1, ]
trajCombined2 <- trajGSM3
combinedMat2 <- plotTrajectoryHeatmap(trajCombined2, returnMat = T, varCutOff = 0)

rowOrder <- match(rownames(combinedMat2), rownames(trajGSM3))
ht1 <- plotTrajectoryHeatmap(trajGSM3, labelTop = 1, pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajGIM3, labelTop = 1, pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht3 <- plotTrajectoryHeatmap(trajMM3, labelTop = 1, pal = paletteContinuous(set = "horizonExtra"), varCutOff = 0, rowOrder = rowOrder)
pdf("combined_R-0.5.pdf", width = 20, height = 8)
ComplexHeatmap::draw(ht1 + ht2 + ht3)
dev.off()

#################### plot 2000 RNA HVGs peak accessibility ###########################

data <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor2_motif3.rds")

# plot motifMatrix
trajMM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
# plot geneScoreMatrix 
trajGSM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
# plot geneIntegrationMatrix 
trajGIM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)

gene <- read.table("target_cluster6.txt", header = F, stringsAsFactors=F)
gene1 <- gene[,1]
pattern <- paste(gene1, collapse = "|")
motifs <- grep(pattern, rownames(trajGSM), value = TRUE)
ordered_motifs <- unlist(sapply(gene1, function(g) grep(g, motifs, value = TRUE)))

# motifs_z <- motifs[grep("^z:", motifs)]
trajGSM_sub <- trajGSM[motifs, ]
trajGIM_sub <- trajGIM[motifs, ]
trajGSM_sub <- trajGSM_sub[ordered_motifs, ]
trajGIM_sub <- trajGIM_sub[ordered_motifs, ]

mat <- t(scale(t(assay(trajGSM_sub))))
p2_2 <- pheatmap(mat, color = paletteContinuous(set = "solarExtra"),
                 cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE)    
mat <- t(scale(t(assay(trajGIM_sub))))      
p3_2 <- pheatmap(mat, color = paletteContinuous(set = "solarExtra"),
                 cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE)     

pdf("target_cluster6.pdf", width = 6, height = 8)
print(p2_2)
print(p3_2)
dev.off()

# combine geneScoreMatrix and geneInteration matrix
gene <- read.table("target_cluster6.txt", header = F, stringsAsFactors=F)
gene1 <- gene[,1]
pattern <- paste(gene1, collapse = "|")
motifs <- grep(pattern, rownames(trajGSM), value = TRUE)
ordered_motifs <- unlist(sapply(gene1, function(g) grep(g, motifs, value = TRUE)))

# motifs_z <- motifs[grep("^z:", motifs)]
trajGSM2 <- trajGSM[motifs, ]
trajGIM2 <- trajGIM[motifs, ]

corGSM_GIM <- correlateTrajectories(trajGSM2, trajGIM2, corCutOff = -1, varCutOff1 = 0, varCutOff2 = 0)
a <- as.data.frame(corGSM_GIM[[1]][,c(3,7)])

#====================================================================================================================#
###                               plot ICE1 and BZR1 target genes peak accessibility                               ###
#====================================================================================================================#

data <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor2_motif3.rds")
gene <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/ICE_BZR_target.gene.txt", header = F, stringsAsFactors=F)

projAll <- readRDS("ArchR_EndospermD_rmTorpedo_reAnchor2_motif3.rds")
# plot motifMatrix
trajMM  <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
# plot geneScoreMatrix 
trajGSM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
# plot geneIntegrationMatrix 
trajGIM <- getTrajectory(ArchRProj = projAll, name = "EndospermTrajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)



