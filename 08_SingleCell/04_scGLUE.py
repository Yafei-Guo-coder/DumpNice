import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import scglue
import seaborn as sns
# from IPython import display
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
from glob import glob

rna = ad.read_h5ad("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/integration.spectral/endo_rna_embeddings.h5ad")
atac = ad.read_h5ad("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/integration.spectral/endo_atac_embeddings.h5ad")
guidance_hvf = nx.read_graphml("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/integration.spectral/endo_guidance.graphml.gz")

atac.var["name"] = atac.var_names
rna.var["name"] = rna.var_names
# 这里也是 更新一下peak(V3)
peaks = atac.var.query("highly_variable").index #V1

# 找到所有 *_specific_peaks.txt 文件
# files = glob("*_specific_peaks.txt")
# all_peaks = []
# 
# for f in files:
#   df = pd.read_csv(f, sep="\t")
# 
# peaks = df["seqnames"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
# all_peaks.extend(peaks)
# all_peaks = sorted(set(all_peaks))
# peaks = pd.Index(all_peaks)

# 这里基因的高变基因没有关注的随时序变化的基因，所以直接使用随时序变化的基因(V2)
# genes = rna.var.query("highly_variable").index #V1
# new_genes = pd.read_csv("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/Figs/Fig.Endo_dev/03.dtw_cluster/0917.HVG2000_2nd.txt", header=None, names=["gene"] )
new_genes = new_genes["gene"]
valid_new_genes = new_genes[new_genes.isin(rna.var_names)]
genes = valid_new_genes
genes = genes.values 
genes = pd.Index(genes) 


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
  if attr["pval"] < 0.071
)
print(gene2peak)

# genes_to_check = ["AT3G26744", "AT1G75080", "AT5G67580"]
# for g in genes_to_check:
#     if g in gene2peak.nodes:
#         print(f"{g} 在 gene2peak 里")
#     else:
#         print(f"{g} 不在 gene2peak 里")
# 
# edges_list = []
# for g in ["AT3G26744", "AT1G75080", "AT5G67580,'AT1G01150', 'AT1G01160', 'AT1G04013', 'AT1G04007', 'AT1G01180', 'AT1G01170']:
#   for _, tgt, attr in gene2peak.out_edges(g, data=True):
#     edges_list.append({"TF": g, "Peak": tgt, "qval": attr["qval"]})

# df = pd.DataFrame(edges_list)
# print(df)
# 基于推断的顺式调控区域构建TF-靶基因调控网络

# motif_bed = scglue.genomics.read_bed("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/db.for_GRN/ath_cisbp_all.GeneID.fixed.bed")
# 上面这个bed文件缺少预测的ICE1的motif基因组位置，所以重新预测一下ICE1的bed然后添加到这个bed文件里面
# fimo --o output --motif AT3G26744 /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/motif/Ath_TF_binding_motifs_addV2.meme /jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/00_ref/Arabidopsis/TAIR10_chr_all.fa

motif_bed = scglue.genomics.read_bed("/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/seed/01_fragment/03_ATACseed_run2/anchor_ED/network/all_fimo.V2.bed")

motif_bed.head()
tfs = pd.Index(motif_bed["name"]).intersection(rna.var_names)
# tfs.size
# "AT3G26744" in tfs
# "AT1G75080" in tfs

rna[:, np.union1d(genes, tfs)].write_loom("rna.loom")
np.savetxt("tfs.txt", tfs, fmt="%s")

!pyscenic grn rna.loom tfs.txt -o draft_grn.csv --seed 0 --num_workers 20 --cell_id_attribute cellID --gene_attribute gene_id

# 经由ATAC峰生成TF顺式调控排序
peak_bed = scglue.genomics.Bed(atac.var.loc[peaks])
peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)

# print("节点数:", len(peak2tf.nodes))
# print("边数:", len(peak2tf.edges))
# print(list(peak2tf.nodes)[:5])
# 
# # 随机查看 5 条边及属性
# for u, v, attr in list(peak2tf.edges(data=True))[:5]:
#   print(u, "→", v, attr)
# 
# "Chr1:38668-39168" in peak2tf.nodes  # 检查某个 peak
# "AT3G26744" in peak2tf.nodes

peak2tf_filtered = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)
print("节点数:", len(peak2tf_filtered.nodes))
print("边数:", len(peak2tf_filtered.edges))
print(list(peak2tf_filtered.nodes)[:5])
"Chr1:38668-39168" in peak2tf_filtered.nodes  # 检查某个 peak
"AT3G26744" in peak2tf_filtered.nodes 
# 随机查看 5 条边及属性
# for u, v, attr in list(peak2tf_filtered.edges(data=True))[:5]:
#   print(u, "→", v, attr)
# 
# "Chr1:38668-39168" in peak2tf_filtered.nodes  # 检查某个 peak
# "AT3G26744" in peak2tf_filtered.nodes

# peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)

gene2tf_rank_glue = scglue.genomics.cis_regulatory_ranking(
  gene2peak, peak2tf_filtered, genes, peaks, tfs,
  region_lens=atac.var.loc[peaks, "chromEnd"] - atac.var.loc[peaks, "chromStart"],
  random_state=0
)
gene2tf_rank_glue.iloc[:5, :5]

flank_bed = scglue.genomics.Bed(rna.var.loc[genes]).strand_specific_start_site().expand(1000, 500)
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

# 检查文件里面我们关注的基因都在不在
import pandas as pd
glue_df = pd.read_feather("glue.genes_vs_tracks.rankings.feather")
supp_df = pd.read_feather("supp.genes_vs_tracks.rankings.feather")
print(glue_df.head())
print(supp_df.head())
print(glue_df.columns)
print(glue_df.shape)

genes = ['AT1G75080', 'AT3G26744']

# 检查 tracks 列是否包含这些基因（假设 tracks 值为基因ID + '_glue'）
for gene in genes:
  gene_with_glue = f"{gene}_glue"  # 构造类似 AT1G75080_glue 的值
if gene_with_glue in glue_df['tracks'].values:
  print(f"{gene} (as {gene_with_glue}) is present in the 'tracks' column.")
else:
  print(f"{gene} (as {gene_with_glue}) is NOT present in the 'tracks' column.")



tfs_of_interest = ["AT3G26744", "AT1G75080"]

# 统计每个 TF 在 glue feather 中小于阈值的 target 数量
tfs_of_interest = ["AT3G26744", "AT1G75080"]
rank_threshold = 800

for tf in tfs_of_interest:
  count_glue = (glue_df[tf] <= rank_threshold).sum()
  count_supp = (supp_df[tf] <= rank_threshold).sum()
  print(f"{tf}: glue={count_glue}, supp={count_supp}, total={count_glue+count_supp}")



print("AT1G75080 在列中:", "AT1G75080" in glue_df.columns)
print(glue_df["AT1G75080"].sort_values().head(10))

ann = pd.read_csv("ctx_annotation.tsv", sep="\t")
print(ann[ann["gene_name"] == "AT1G75080"])



import pandas as pd
annotations = pd.read_csv("ctx_annotation.tsv", sep="\t")
for tf in ['AT3G26744', 'AT1G75080']:
  if tf in annotations['gene_name'].values:
    print(f"{tf} has {len(annotations[annotations['gene_name'] == tf])} motif annotations")
else:
  print(f"{tf} is missing motif annotations")



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

pyscenic ctx draft_grn.csv glue.genes_vs_tracks.rankings.feather supp.genes_vs_tracks.rankings.feather --annotations_fname ctx_annotation.tsv --expression_mtx_fname rna.loom --output pruned_grn_1000.csv  --rank_threshold  1000 --min_genes 1 --num_workers 20 --cell_id_attribute cellID --gene_attribute gene_id 2> /dev/null
pyscenic ctx draft_grn.csv glue.genes_vs_tracks.rankings.feather supp.genes_vs_tracks.rankings.feather --annotations_fname ctx_annotation.tsv --expression_mtx_fname rna.loom --output pruned_grn_1000_gene100.csv --rank_threshold 1000 --min_genes 100 --num_workers 20 --cell_id_attribute cellID --gene_attribute gene_id 2> /dev/null
pyscenic ctx draft_grn.csv glue.genes_vs_tracks.rankings.feather supp.genes_vs_tracks.rankings.feather --annotations_fname ctx_annotation.tsv --expression_mtx_fname rna.loom --output pruned_grn_1200.csv --rank_threshold 1200 --min_genes 1 --num_workers 20 --cell_id_attribute cellID --gene_attribute gene_id 2> /dev/null
pyscenic ctx draft_grn.csv glue.genes_vs_tracks.rankings.feather supp.genes_vs_tracks.rankings.feather --annotations_fname ctx_annotation.tsv --expression_mtx_fname rna.loom --output pruned_grn_1500.csv --rank_threshold 1500 --min_genes 1 --num_workers 20 --cell_id_attribute cellID --gene_attribute gene_id 2> /dev/null
pyscenic ctx draft_grn.csv glue.genes_vs_tracks.rankings.feather supp.genes_vs_tracks.rankings.feather --annotations_fname ctx_annotation.tsv --expression_mtx_fname rna.loom --output pruned_grn_1800.csv --rank_threshold 1800 --min_genes 1 --num_workers 20 --cell_id_attribute cellID --gene_attribute gene_id 2> /dev/null

grn = scglue.genomics.read_ctx_grn("pruned_grn_1000_gene100.csv")
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(grn, seed=42)
nx.draw(grn, pos=pos, with_labels=True, node_size=500, font_size=10)
plt.savefig("gene_TF_network_1000_gene100.pdf", format="pdf")
plt.close()
















# gemi 
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import feather
import loompy
import matplotlib.pyplot as plt

# --- 1. Configuration ---
# Target TFs to investigate (Arabidopsis AGI IDs)
TFS_OF_INTEREST = ['AT1G75080', 'AT3G26744'] 
# Maximum rank to consider for a target gene
RANK_THRESHOLD = 1000 
# Input File Paths
FEATHER_PATH = 'glue.genes_vs_tracks.rankings.feather'
LOOM_PATH = 'rna.loom'
ADATA_PATH = '/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250827.Endo.snATAC_snRNA.GRN/integration.spectral/endo_rna_embeddings.h5ad' 

# --- 2. Load and Fix Data ---

# Load AnnData object
adata = sc.read_h5ad(ADATA_PATH) 

# Load Feather file
feather_df = glue_df = pd.read_feather("glue.genes_vs_tracks.rankings.feather")

# Extract gene IDs from the Loom file (assuming gene IDs are under the 'gene_id' attribute)
with loompy.connect(LOOM_PATH) as ds:
  gene_ids = ds.ra['gene_id'] 

# Set the correct AGI gene IDs as the index of the Feather DataFrame
if len(gene_ids) == len(feather_df):
  feather_df.index = gene_ids
else:
  # Exit if row counts don't match, indicating a serious file mismatch
  print("FATAL ERROR: Gene count mismatch between Loom file and Feather file.")
exit()

# --- 3. Calculate and Visualize Regulon Scores ---

for tf in TFS_OF_INTEREST:
  if tf not in feather_df.columns:
  print(f"Skipping TF {tf}: Not found in the Feather file columns.")
continue

# 3.1. Extract Target Genes (now using AGI IDs)
# Find indices (AGI IDs) where rank is below the threshold
target_genes = feather_df.index[feather_df[tf] <= RANK_THRESHOLD].tolist()

# Filter for genes that exist in the AnnData object
valid_genes = [g for g in target_genes if g in adata.var_names]

print(f"\nProcessing TF: {tf}. Valid targets found: {len(valid_genes)}")

if len(valid_genes) >= 5: 
  
  # 3.2. Compute Module Score (Regulon Activity)
  score_key = f'{tf}_Regulon_Score'
sc.tl.score_genes(
  adata, 
  gene_list=valid_genes, 
  ctrl_size=50, 
  score_name=score_key,
  use_raw=True 
)

# 3.3. Visualize Regulon Activity on UMAP
sc.pl.umap(
  adata, 
  color=score_key, 
  title=f'Regulon Score for {tf} (Targets N={len(valid_genes)})',
  save=f'_{tf}_regulon_score_check.pdf',
  show=True
)

else:
  print(f"Skipping visualization for {tf}: Insufficient valid target genes (< 5).")

print("\n--- Regulon Quality Check Complete ---")
print("Check the generated PDF files for Regulon activity on the UMAP plot.")
















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
