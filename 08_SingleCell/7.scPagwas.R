library(Seurat)
library(scPagwas, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib")

# 单细胞数据处理
# input单核转录组
silique <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/silique_final202502.rds")
shoot <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/shoot_final202502.rds")
root <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/root_final202502.rds")
FM <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/FM_IM_final202502.rds")
flower <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/Flower_final202502.rds")
leaf <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/leaf_final202502.rds")
sam <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/sam_final202502.rds")

silique <- RenameCells(silique, add.cell.id = "silique")
shoot <- RenameCells(shoot, add.cell.id = "shoot")
root <- RenameCells(root, add.cell.id = "root")
FM <- RenameCells(FM, add.cell.id = "FM")
flower <- RenameCells(flower, add.cell.id = "flower")
leaf <- RenameCells(leaf, add.cell.id = "leaf")
sam <- RenameCells(sam, add.cell.id = "sam")

object.list <- list(silique = silique, shoot = shoot, root = root, FM = FM, flower = flower, leaf = leaf, sam = sam)
# 标准化 + 选择特征
object.list <- lapply(object.list, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  return(obj)})
# 查找锚点 + 整合
anchors <- FindIntegrationAnchors(object.list = object.list)
combined <- IntegrateData(anchorset = anchors)
# 后处理
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)

combined <- FindVariableFeatures(combined, nfeatures = 8000, selection.method = "vst")
length(VariableFeatures(combined))
saveRDS(combined, "snRNA_combinedAll.rds")


library(Seurat)
library(ff,lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
library(scPagwas, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
# 准备单细胞输入文件
#snRNA <- readRDS("snRNA_combined.rds")
#DefaultAssay(snRNA) <- "RNA"
#Idents(snRNA) <- snRNA$cell_type
#snRNA$celltypes <- snRNA$cell_type


# 准备GWAS文件
# input GWAS数据
# /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/1001genome/vcf/all/maf003_filter_ld_prunin
# 找离遗传变异位点最近的基因
# /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/01_ref/Araport11_GTF_genes.gtf
# /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/1001genome/vcf/all/maf005_filter_ld_prunin.bed
# /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/1-scPagwas/gwas/output/elev.gcta_with_gene.txt 

# 准备KEGG的输入文件
# wget -O ath_pathways.txt http://rest.kegg.jp/list/pathway/ath
# library(KEGGREST)
# library(dplyr)
# pathways <- read.delim("ath_pathways.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# pathway_ids <- pathways$V1
# gene2pathway <- lapply(pathway_ids, function(pid) {
#   res <- tryCatch({
#     df <- keggLink("ath", pid)
#     data.frame(
#       gene_id = sub("ath:", "", df),
#       pathway_id = pid,
#       stringsAsFactors = FALSE
#     )}, error = function(e) NULL)
#   return(res)})
# gene2pathway_df <- bind_rows(gene2pathway)
# write.csv(gene2pathway_df, "ath_gene2kegg_pathways.csv", row.names = FALSE)
# kegg_df <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/1-scPagwas/kegg/ath_gene2kegg_pathways.txt", header = TRUE, sep = "\t",stringsAsFactors = F) 
kegg_df <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/1-scPagwas/At.KEGG/ath_gene2kegg_pathways.txt", header = TRUE, sep = "\t",stringsAsFactors = F)
Genes_by_pathway_kegg <- split(kegg_df$gene_id, kegg_df$pathway_id)


# 准备基因注释信息
gtf_df <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/01_ref/Araport11_GTF_genes_sort_addchr.gtf")
gtf_df <- as.data.frame(gtf_df)
block_annotation <- gtf_df[,c(1,4,5,13)]
colnames(block_annotation)<-c("chrom", "start","end","label")
#save(block_annotation,file="block_annotation.RData")

# 准备LD文件
# plink1.9 --bfile maf005_filter_ld_prunin2_ld --allow-no-sex --autosome --r2 --ld-window-kb 1000 --ld-window-r2 0.2 --out maf005_filter_ld_prunin_ld
# covid_ld <-read.delim("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/1001genome/vcf/all/maf005_filter_ld_prunin2_ld.ld")
# colnames(covid_ld)[7]<-"R"
# print out the result in chrom number
# lapply(unique(covid_ld$CHR_A), function(i){
#   a<- covid_ld[covid_ld$CHR_A == i, ]
#   file_name <- paste0("./ld/",i,".Rds")
#   saveRDS(a, file = file_name)
# })
# integrate the data
ld_folder <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/1-scPagwas/ld/"
r2_threshold <- 0.2 
chrom_ld <- lapply(as.character(1:5), function(chrom){
  chrom_ld_file_path <- paste0(ld_folder, '/', chrom, '.Rds')
  ld_data <- readRDS(chrom_ld_file_path)
  ld_data_filtered <- ld_data[ld_data$R^2 > r2_threshold, c("SNP_A", "SNP_B", "R")]
  return(ld_data_filtered)})
# 运行
#source("/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/scPagwas/R/Singlecell_heritability_contributions.R")

#source("/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/scPagwas/R/Pathway_annotation_input.R")
#source("/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/scPagwas/R/Link_pathway_blocks_gwas.R")

Pagwas <- scPagwas_main(Pagwas =NULL,
                        gwas_data="/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/1-scPagwas/gwas/output/temp11.mf005.gcta_with_gene_mf005.txt",
                        Single_data ="/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/1-scPagwas/snRNA_combined_out.8kHV.rds",
                        output.prefix="temp11",
                        output.dirs="temp11",
                        Pathway_list=Genes_by_pathway_kegg,
                        n.cores=2,
                        assay="RNA",
                        singlecell=T, 
                        iters_singlecell = 100,
                        celltype=T,
                        block_annotation = block_annotation,
                        chrom_ld = chrom_ld)

save(Pagwas,file="./temp11_scPagwas.RData")
print(Pagwas)
names(Pagwas@misc)

# 绘图
# load("./temp11_scPagwas.RData")
library(ff,lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
library(scPagwas, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
library(grImport2, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
require("RColorBrewer")
require("Seurat")
require("SeuratObject")
require("ggsci")
require("dplyr")
require("ggplot2")
require("ggpubr")

Pagwas$positivecell <- ifelse(Pagwas$Random_Correct_BG_adjp < 0.05, 1, 0)

color29 <- c("#3CCF4E","#4D96FF","#8338EC","#D9DD6B","#D54C4C","#FDD2BF","#492F10","#00A19D","#5D8233","#284E78","#3E215D","#835151",
             "#F08FC0","#C6B4CE","#BB8760","#558776","#E99497","#FFBD9B","#0A1D37","#7FB285","#BDB2FF","#FF9B54","#FF006E","#FFD93D",
             "#247BA0","#E63946","#8AC926","#06D6A0","#F4A261")


pdf("DimPlot_singlecell_data.pdf",width = 10, height = 8)
Seurat::DimPlot(Pagwas,pt.size=1,reduction="umap",label = T, repel=TRUE)+
  scale_colour_manual(name = "cell_type", values = color29)+
  umap_theme()+
  ggtitle("root_shoot_silique")+
  labs(x="UMAP",y="")+
  theme(aspect.ratio=1)
dev.off()

pdf("DimPlot_scPagwas_gPAS_score.pdf",width = 10, height = 8)
scPagwas_Visualization(Single_data=Pagwas,
                       p_thre = 0.05,
                       FigureType = "umap",
                       width = 7,
                       height = 7,
                       lowColor = "white", 
                       highColor = "red",
                       output.dirs="figure",
                       size = 0.5,
                       do_plot = T)
dev.off()

pdf("BarPlot_proportion_positiveCell.pdf",width = 10, height = 8)
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                         var_ident="cell_type",#you should change this to you specific columns
                         var_group="positiveCells",
                         vec_group_colors=c("#E8D0B3","#7EB5A6"),
                         do_plot = T)
dev.off()

pdf("BarPlot_proportion_positiveNagtiveCell.pdf")
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                         var_ident="positiveCells",
                         var_group="cell_type", #you should change this to you specific columns
                         p_thre = 0.01,
                         vec_group_colors=NULL,
                         f_color=colorRampPalette(brewer.pal(n=10, name="RdYlBu")),
                         do_plot = T)
dev.off()

pdf("DotPlot_heritability_correlated_genes.pdf")
heritability_cor_scatterplot(gene_heri_cor=Pagwas@misc$PCC,
                             topn_genes_label=10,
                             color_low="#035397",
                             color_high ="#F32424",
                             color_mid = "white",
                             text_size=2,
                             do_plot=T,
                             max.overlaps =20,
                             width = 7,
                             height = 7)
dev.off()



Index <- estimate <- lower <- upper <- label <- Pvalue <- NULL
bootstrap_results <- Pagwas@misc$bootstrap_results[-1, c(
  "bp_value",
  "bias_corrected_estimate",
  "CI_lo",
  "CI_hi"
)]
bootstrap_results <- bootstrap_results[order(
  bootstrap_results$bias_corrected_estimate,
  decreasing = T
), ]

dat <- data.frame(
  Index = seq_len(nrow(bootstrap_results)),
  label = rownames(bootstrap_results),
  estimate = bootstrap_results$bias_corrected_estimate,
  lower = bootstrap_results$CI_lo,
  upper = bootstrap_results$CI_hi,
  #转化为科学计数法
  Pvalue = bootstrap_results$bp_value
)

plot1 <- ggplot2::ggplot(dat, aes(
  y = Index,
  x = estimate
)) +
  geom_errorbarh(aes(
    xmin = lower,
    xmax = upper
  ),
  color = "#6D8299",
  height = 0.25
  ) +
  geom_point(
    shape = 18,
    size = 5,
    color = "#D57E7E"
  ) +
  geom_vline(
    xintercept = 0,
    color = "#444444",
    linetype = "dashed",
    cex = 1,
    alpha = 0.5
  ) +
  scale_y_continuous(
    name = "", breaks = seq_len(nrow(dat)),
    labels = dat$label,
    trans = "reverse"
  ) +
  xlab("bias_corrected_estimate (95% CI)") +
  ylab(" ") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x.bottom = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12, colour = "black")
  )

table_base <- ggplot2::ggplot(dat, aes(y = label)) +
  ylab(NULL) +
  xlab("  ") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(
      color = "white",
      hjust = -3,
      size = 25
    ), ## This is used to help with alignment
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank()
  )

tab1 <- table_base +
  labs(title = "space") +
  geom_text(aes(
    y = rev(Index), x = 1,
    label = sprintf("%.3f", round(Pvalue, digits = 3))
  ),
  size = 4
  ) + ## decimal places
  ggtitle("Pvalue")


pdf("Forest_estimateplot_1.pdf",width = 9,height = 7)
print(plot1)
dev.off()
pdf("Forest_estimateplot_2.pdf",width = 3,height = 7)
print(tab1)
dev.off()


source(system.file("extdata", "plot_scpathway_contri_dot.R", package = "scPagwas"))
library(tidyverse)
library("rhdf5")
library(ggplot2)
library(grDevices)
library(stats)
library(FactoMineR)
library(scales)
library(reshape2)
library(ggdendro)
library(grImport2, lib.loc = "/jdfsbjcas1/ST_BJ/PUB/User/guoyafei1/Tools/R/lib2")
library(gridExtra)
library(grid)
#library(sisal)

type <- names(table(Pagwas$cell_type))
pdf("plot_scpathway_dot.pdf", width = 25, height = 6)
plot_scpathway_dot(Pagwas=Pagwas,
                   celltypes=type,
                   topn_path_celltype=5,
                   filter_p=0.05,
                   max_logp=15,
                   display_max_sizes=F,
                   size_var ="logrankPvalue" ,
                   col_var="proportion",
                   shape.scale = 8,
                   cols.use=c("lightgrey", "#E45826"),
                   dend_x_var = "logrankPvalue",
                   dist_method="euclidean",
                   hclust_method="ward.D",
                   do_plot = T,
                   width = 7,
                   height = 7)
dev.off()                   



