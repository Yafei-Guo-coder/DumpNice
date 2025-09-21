import matplotlib.patches as mpatches
import numpy as np
import stereo as st
import os
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from stereo.plots import PlotCollection
import pandas as pd



# from stereopy.tools import PlotCollection
# read the GEF file
# data_path = '/home/xiaomanyi/1_soybean_spatial/merge/cd_07/cd_07.merged.gem.gz'
# data_path = '/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/03_soybean_nodule/01_stero/addUSDA110/merge/cd_07/cd_07.merged.gem.gz'
# data = st.io.read_gem(file_path=data_path, bin_size=50)
# data.tl.cal_qc()
# data.tl.raw_checkpoint()
# remember to set flavor as seurat
# adata = st.io.stereo_to_anndata(data,flavor='seurat',output='cd_07.merged.h5ad')
# read data
input_file = '/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/03_soybean_nodule/01_stero/01_stereopy/cd_07.merged.h5ad'
data = st.io.read_ann_h5ad(input_file, spatial_key='spatial')
# preprocessing
data.tl.cal_qc()
# data.tl.filter_cells(min_counts=100, pct_counts_mt=10)
data.tl.filter_genes(min_cells=3)
data.tl.raw_checkpoint() # save raw data
data.tl.normalize_total(target_sum=10000)
data.tl.log1p()
# PCA
data.tl.pca(use_highly_genes=False, n_pcs=50, svd_solver='arpack')
# Gaussian smoothing
data.tl.gaussian_smooth(n_neighbors=10, smooth_threshold=90)
data.tl.normalize_total(target_sum=10000)
data.tl.log1p()
# clustering
data.tl.pca(use_highly_genes=False, n_pcs=50, svd_solver='arpack')
data.tl.neighbors(pca_res_key='pca', n_pcs=30, res_key='neighbors')

# test the resolution parameter
resolutions = [0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.4, 0.5]
cluster_counts = []
for res in resolutions:
    data.tl.leiden(
        neighbors_res_key="neighbors", 
        res_key=f"leiden_{res}", 
        resolution=res
    )
    n_clusters = len(data.cells[f"leiden_{res}"].unique())
    cluster_counts.append(n_clusters)
    print(f"Resolution {res}: {n_clusters} clusters")

# Get coordinates
x = data.position[:, 0]
y = data.position[:, 1]

for res in resolutions:
    cluster_key = f"leiden_{res}"
    cluster = data.cells[cluster_key].astype(int)
    # Plot
    plt.figure(figsize=(8, 3))
    sc = plt.scatter(x, y, c=cluster, cmap="tab20", s=1, alpha=1.0)
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.title(f"Spatial Cluster: resolution = {res}")
    # Legend
    unique_clusters = np.unique(cluster)
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))
    patches = [
        mpatches.Patch(color=colors[i], label=f"Cluster {c}")
        for i, c in enumerate(unique_clusters)
    ]
    plt.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.15, 1.0),
               fontsize=6, title="Cluster", frameon=False)
    # Save plot
    plt.savefig(f"stereopy_clusters_res_{res}.png", dpi=300, bbox_inches='tight')
    plt.show()

plt.figure(figsize=(6, 4))
plt.plot(resolutions, cluster_counts, marker='o')
plt.xlabel("Resolution")
plt.ylabel("Number of clusters")
plt.title("Cluster number vs. Resolution")
plt.grid(True)
# Save the figure to a file, e.g., PNG format
plt.savefig("cluster_number_vs_resolution.png", dpi=300, bbox_inches='tight')
plt.show()

# write the result data to the output file (error)
st.io.write_h5ad(data, use_raw=True, use_result=True, key_record=None, output='./cd_07_paraTest.h5ad')

# use the resolution 0.1 to detect marker genes
data.tl.leiden(neighbors_res_key="neighbors", res_key="leiden", resolution=0.1)

data.tl.highly_variable_genes(
    method='seurat',          # 推荐：经典 Seurat 方法
    n_top_genes=2000,         # 保留前2000个变异最大的基因
    res_key='highly_variable_genes'
)

# Find marker genes (SVGs) per cluster
data.tl.find_marker_genes(
    cluster_res_key='leiden',
    method='t_test', # 't_test', 'wilcoxon_test', 'logreg'
    case_groups='all',
    control_groups='rest',
    corr_method='benjamini-hochberg',
    use_raw=True,
    use_highly_genes=False,
    hvg_res_key='highly_variable_genes',
    res_key='marker_genes',
    n_genes='all',
    n_jobs=4
)

marker_dict = data.tl.result['marker_genes']
sig_gene_list = []
for key, df in marker_dict.items():
    if not isinstance(df, pd.DataFrame):
        continue
    if not key.endswith('.vs.rest'):
        continue
    if 'pvalues_adj' not in df.columns:
        print(f"Warning: {key} missing 'pvalues_adj'")
        continue
    # 提取 cluster 名（如 '2.vs.rest' => '2'）
    cluster = key.split('.')[0]
    # 筛选显著 marker genes
    df_sig = df[df['pvalues_adj'] < 0.05].copy()
    df_sig['cluster'] = cluster
    sig_gene_list.append(df_sig)

# 合并所有结果
sig_marker_df = pd.concat(sig_gene_list, ignore_index=True)

# 可选：保存或查看
print(sig_marker_df.head())
sig_marker_df.to_csv("filtered_marker_genes2.csv", index=False)


# 实例化
pc = PlotCollection(data)

# 调用方法绘图
pc.marker_genes_heatmap(
    res_key='marker_genes',               # 与 find_marker_genes 时使用的参数一致
    cluster_res_key='leiden',             # 聚类的 key，一般是 'leiden'
    markers_num=50,                       # 每个 cluster 取前 10 个基因
    sort_key='scores',                    # 按 scores 排序
    ascend=False,                         # 从高到底
    show_labels=True,
    show_group=True,
    show_group_txt=True,
    do_log=True,                          # 如果你的数据是 log1p 归一化的，设为 True
    width=800,
    height=600,
    out_path='marker_gene_heatmap2.png',   # 输出文件名
    out_dpi=300                           # 分辨率
)

pc.marker_genes_heatmap(
    res_key='marker_genes',
    cluster_res_key='leiden',
    gene_list=['GmCH_02G027480', 'GmCH_10G262450', 'GmCH_20G561760', 'GmCH_02G053760', 'GmCH_07G175080', 'GmCH_06G162110', 'GmCH_08G229910', 'GmCH_17G456780', 'GmCH_06G154030', 'GmCH_12G344250', 'GmCH_13G354660', 'GmCH_18G503850', 'GmCH_10G284570', 'GmCH_02G035860', 'GmCH_11G315630', 'GmCH_13G354640', 'GmCH_07G196560', 'GmCH_17G455200'],  # replace with real gene names
    width=800,
    height=600,
    out_path='rna_marker_gene_heatmap.png',   # optional: save to file
    out_dpi=300  
)


pc.marker_genes_scatter(
    res_key='marker_genes',
    markers_num=20,
    values_to_plot='logfoldchanges',
    sort_by='scores',
    out_path='marker_genes_scatter.png',
    out_dpi=300,
    width=1200,
    height=800
)

# 调用实例方法
pc.marker_genes_scatter(
    res_key='marker_genes',
    genes=[
        'GmCH_02G027480', 'GmCH_10G262450', 'GmCH_20G561760', 'GmCH_02G053760',
        'GmCH_07G175080', 'GmCH_06G162110', 'GmCH_08G229910', 'GmCH_17G456780',
        'GmCH_06G154030', 'GmCH_12G344250', 'GmCH_13G354660', 'GmCH_18G503850',
        'GmCH_10G284570', 'GmCH_02G035860', 'GmCH_11G315630', 'GmCH_13G354640',
        'GmCH_07G196560', 'GmCH_17G455200'
    ],
    values_to_plot='scores',   # 或 'scores', 'logfoldchanges'
    out_path='marker_genes_rna_scatter.png',
    out_dpi=300,
    width=1000,
    height=600
)

# 生成 marker genes 散点图，保存到文件
pc.marker_genes_text(
    res_key="marker_genes",
    groups='all',          # 显示所有cluster
    markers_num=20,        # 每个cluster显示20个marker基因
    sort_key='scores',     # 以scores排序
    ascend=False,          # 降序排序，分数高的排前面
    fontsize=8,            # 字体大小
    ncols=4,               # 每行4列子图
    sharey=True,           # 共享y轴刻度
    width=2400,            # 图宽，像素
    height=2400,            # 图高，像素
    out_path="./marker_genes_scatter.png",  # 保存文件路径
    out_dpi=300            # 保存文件分辨率，300dpi
)


pc.spatial_scatter_by_gene(
    gene_name=['GmCH_02G027480', 'GmCH_10G262450', 'GmCH_20G561760', 'GmCH_02G053760', 'GmCH_07G175080', 'GmCH_06G162110', 'GmCH_08G229910', 'GmCH_17G456780', 'GmCH_06G154030', 'GmCH_12G344250', 'GmCH_13G354660', 'GmCH_18G503850', 'GmCH_10G284570', 'GmCH_02G035860', 'GmCH_11G315630', 'GmCH_13G354640', 'GmCH_07G196560', 'GmCH_17G455200'],  # 多基因
    dot_size=1,
    palette='CET_L4',
    color_bar_reverse=True,
    width=2000,
    height=2400,
    out_path='multi_gene_spatial.png',  # 多图会自动拼图排列
    out_dpi=300
)

