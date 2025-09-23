
import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import scglue
import seaborn as sns
from IPython import display
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from networkx.drawing.nx_agraph import graphviz_layout

rna = ad.read_h5ad("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/integration.spectral/endo_rna_embeddings.h5ad")
atac = ad.read_h5ad("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/integration.spectral/endo_atac_embeddings.h5ad")
guidance_hvf = nx.read_graphml("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/integration.spectral/endo_guidance.graphml.gz")

rna.var["name"] = rna.var_names
atac.var["name"] = atac.var_names

genes = rna.var.query("highly_variable").index
peaks = atac.var.query("highly_variable").index
# 用GLUE特征嵌入进行顺式调控推断
features = pd.Index(np.concatenate([rna.var_names, atac.var_names]))
feature_embeddings = np.concatenate([rna.varm["X_glue"], atac.varm["X_glue"]])

skeleton = guidance_hvf.edge_subgraph(
    e for e, attr in dict(guidance_hvf.edges).items()
    if attr["type"] == "fwd"
).copy()

reginf = scglue.genomics.regulatory_inference(
    features, feature_embeddings,
    skeleton=skeleton, random_state=0
)

gene2peak = reginf.edge_subgraph(
    e for e, attr in dict(reginf.edges).items()
    if attr["pval"] < 0.05
)
print(gene2peak)

# 基于推断的顺式调控区域构建TF-靶基因调控网络

motif_bed = scglue.genomics.read_bed("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/db.for_GRN/ath_cisbp_all.GeneID.fixed.bed")
motif_bed.head()
tfs = pd.Index(motif_bed["name"]).intersection(rna.var_names)
tfs.size

rna[:, np.union1d(genes, tfs)].write_loom("rna.loom")
np.savetxt("tfs.txt", tfs, fmt="%s")

!pyscenic grn rna.loom tfs.txt \
    -o draft_grn.csv --seed 0 --num_workers 20 \
    --cell_id_attribute cellID --gene_attribute gene_id

# 经由ATAC峰生成TF顺式调控排序
peak_bed = scglue.genomics.Bed(atac.var.loc[peaks])
peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)

gene2tf_rank_glue = scglue.genomics.cis_regulatory_ranking(
    gene2peak, peak2tf, genes, peaks, tfs,
    region_lens=atac.var.loc[peaks, "chromEnd"] - atac.var.loc[peaks, "chromStart"],
    random_state=0
)
gene2tf_rank_glue.iloc[:5, :5]

flank_bed = scglue.genomics.Bed(rna.var.loc[genes]).strand_specific_start_site().expand(500, 500)
flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)
gene2flank = nx.Graph([(g, g) for g in genes])
gene2tf_rank_supp = scglue.genomics.cis_regulatory_ranking(
    gene2flank, flank2tf, genes, genes, tfs,
    n_samples=0
)
gene2tf_rank_supp.iloc[:5, :5]


# 使用顺式调控排序修剪共表达网络
gene2tf_rank_glue.columns = gene2tf_rank_glue.columns + "_glue"
gene2tf_rank_supp.columns = gene2tf_rank_supp.columns + "_supp"
scglue.genomics.write_scenic_feather(gene2tf_rank_glue, "glue.genes_vs_tracks.rankings.feather")
scglue.genomics.write_scenic_feather(gene2tf_rank_supp, "supp.genes_vs_tracks.rankings.feather")
pd.concat([
    pd.DataFrame({
        "#motif_id": tfs + "_glue",
        "gene_name": tfs
    }),
    pd.DataFrame({
        "#motif_id": tfs + "_supp",
        "gene_name": tfs
    })
]).assign(
    motif_similarity_qvalue=0.0,
    orthologous_identity=1.0,
    description="placeholder"
).to_csv("ctx_annotation.tsv", sep="\t", index=False)

!pyscenic ctx draft_grn.csv \
    glue.genes_vs_tracks.rankings.feather \
    supp.genes_vs_tracks.rankings.feather \
    --annotations_fname ctx_annotation.tsv \
    --expression_mtx_fname rna.loom \
    --output pruned_grn.csv \
    --rank_threshold 500 --min_genes 1 \
    --num_workers 20 \
    --cell_id_attribute cells --gene_attribute genes 2> /dev/null

















cd /jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/GRN

python 04.GRN_step3.GRN_scheme1.py \
    --rna_h5ad ../endo_scrna.h5ad \
    --tf_file /jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/db.for_GRN/TF_list.txt \
    --integration_dir ../integration.spectral \
        --feather_dir ./02.scglue_grn \
    --output_dir ./03.GRN_scheme1 \
    --output_prefix GRN_scheme1 \
    --num_workers 20 \
    --rank_threshold 500 \
    --min_genes 10 \
    --top_genes 2000

### 1. 快速测试模式（最快）
适用于验证流程是否正常，快速获得结果

python 04.GRN_step3.GRN_scheme1.py \
     --rna_h5ad ../endo_scrna.h5ad \
     --tf_file /jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/db.for_GRN/TF_list.txt \
     --integration_dir ../integration.spectral \
     --feather_dir ./02.scglue_grn \
     --output_dir ./03.GRN_scheme1_test \
     --output_prefix GRN_scheme1_test \
     --num_workers 10 \
     --rank_threshold 1000 \
     --min_genes 5 \
     --top_genes 1000

# 超快速测试模式（针对紧急验证）
python 04.GRN_step3.GRN_scheme1.py \
     --rna_h5ad ../endo_scrna.h5ad \
     --tf_file /jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/db.for_GRN/TF_list.txt \
     --integration_dir ../integration.spectral \
     --feather_dir ./02.scglue_grn \
     --output_dir ./03.GRN_scheme1_ultra_fast \
     --output_prefix GRN_scheme1_ultra_fast \
     --num_workers 10 \
     --rank_threshold 2000
     --min_genes 10 
     --top_genes 300 