import scanpy as sc
import pandas as pd
import pyarrow.feather as feather
import os
import sys
import argparse
from multiprocessing import cpu_count
import numpy as np
if not hasattr(np, 'object'):
    np.object = object

from spagrn.regulatory_network import InferNetwork as irn

adata = sc.read_h5ad("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250707.3ST/C05912F1.Cell_type_updated.h5ad")
# 假设你的坐标列是这些（请按你的实际列名改）
adata.obs['x'] = adata.obs['bin_x']
adata.obs['y'] = adata.obs['bin_y']
# 转为 dense 表达矩阵
expr_df = pd.DataFrame(adata.X.toarray(), columns=adata.var_names)
# 提取坐标
coords = adata.obs[['x', 'y']].reset_index(drop=True)
# 合并表达矩阵 + 坐标
feather_df = pd.concat([coords, expr_df.reset_index(drop=True)], axis=1)
# 保存为 spaGRN 格式的 feather 文件
feather.write_feather(feather_df, "spaGRN_input_C05912F1.feather")
print("✅ Feather 文件已保存为 spaGRN_input_arabidopsis.feather")


fn = '/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250707.3ST/C05912F1.Cell_type_updated.h5ad'  # mouse_brain slide_seq data
data = irn.read_file(fn)
print(data)
data = irn.preprocess(data, min_genes=5, min_cells=1, min_counts=5, max_gene_num=8000)
print(data)
sc.tl.pca(data)
# sc.pp.neighbors(data)  # 先计算邻居图
# sc.tl.leiden(data) 
sc.pl.pca(data, color='Cell_type')
tfs_fn = '/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/Ath_TF.txt'
database_fn = 'spaGRN_input_C05912F1.feather'
motif_anno_fn = '/jdfsbjcas1/ST_BJ/P21Z28400N0234/guoyafei1/ATAC_out/ArchR/motif_AT/Ath_TF.tbl'
grn = irn(data)
grn.add_params({'prune_auc_threshold': 0.05, 'rank_threshold': 9000, 'auc_threshold': 0.05})
niche_human = pd.read_csv('resource/lr_network_human.csv')
niche_mouse = pd.read_csv('resource/lr_network_mouse.csv')
# concat the databases together if you are using multiple ligand-receptors databases
# niches = pd.concat([niche_mouse, niche_human])
# print(niches)


grn.infer(database_fn,
          motif_anno_fn,
          tfs_fn,
          # niche_df=niches,
          # gene_list=None,
          num_workers=20,
          cache=True,
          # output_dir='exp/ouput',
          save_tmp=True,
          layers='raw_counts',
          latent_obsm_key='spatial',
          model='bernoulli',
          n_neighbors=10,
          methods=['FDR_I','FDR_C','FDR_G'],
          operation='intersection',
          mode='moran',
          cluster_label='leiden')
          
/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250707.3ST/C05912F1.Cell_type_updated.h5ad
/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250707.3ST/C05912F3.Cell_type_updated.h5ad
/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202208.At_silique_dev/250707.3ST/C05913E5.Cell_type_updated.h5ad

C05912F1鱼雷 C05912F3拐棍 
adata = sc.read_h5ad('project_spagrn.h5ad')
adata
adata.uns['regulon_dict']

import spagrn.plot as prn
import matplotlib.pyplot as plt


#prn.plot_2d(adata, adata.obsm['auc_mtx'], pos_label="spatial", reg_name="AT3G26744", fn='AT3G26744.png')

reg_name = 'AT3G26744(+)'
# 取调控子活性值
values = adata.obsm['auc_mtx'][reg_name]

values_fixed = values.copy()
values_fixed[values_fixed < 0] = 0
values_shift = values_fixed - values_fixed.min() + 1e-6  # 防止log(0)
values_log = np.log(values_shift)

# 对log后的数据再做线性缩放到0-15范围
values_log_scaled = (values_log - values_log.min()) / (values_log.max() - values_log.min()) * 10

# 取空间坐标
x = adata.obs['x']
y = adata.obs['y']
plt.figure(figsize=(16,8))
sc = plt.scatter(x, y, c=values_log_scaled, cmap='plasma', s=0.4, marker='o')
plt.colorbar(sc, label=f'{reg_name} regulon activity')
plt.title('Spatial distribution of ICE1')
plt.xlabel('x')
plt.ylabel('y')
plt.gca().invert_yaxis()  # 空间坐标系中通常y轴反转更符合视觉习惯，可根据实际情况调整
plt.savefig(f'{reg_name}_spatial2.pdf')
plt.show()


reg_name = 'AT1G75080(+)'
# 取调控子活性值
values = adata.obsm['auc_mtx'][reg_name]
import numpy as np

values_fixed = values.copy()
values_fixed[values_fixed < 0] = 0
values_shift = values_fixed - values_fixed.min() + 1e-6  # 防止log(0)
values_log = np.log(values_shift)

# 对log后的数据再做线性缩放到0-15范围
values_log_scaled = (values_log - values_log.min()) / (values_log.max() - values_log.min()) * 10

# 取空间坐标
x = adata.obs['x']
y = adata.obs['y']
plt.figure(figsize=(16,8))
sc = plt.scatter(x, y, c=values_log_scaled, cmap='plasma', s=0.4, marker='o')
plt.colorbar(sc, label=f'{reg_name} regulon activity')
plt.title('Spatial distribution of BZR1')
plt.xlabel('x')
plt.ylabel('y')
plt.gca().invert_yaxis()  # 空间坐标系中通常y轴反转更符合视觉习惯，可根据实际情况调整
plt.savefig(f'{reg_name}_spatial2.pdf')
plt.show()


prn.auc_heatmap(adata,
                adata.obsm['auc_mtx'],
                cluster_label='Cell_type',
                topn=10,
                subset=False,
                save=True,
                fn='C05912F3_clusters_heatmap_top10.png',
                legend_fn="C05912F3_rss_celltype_legend_top10.png")





import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# 1. 选取你关心的10个基因（调控子）名字，比如：
genes = ['AT3G26744(+)', 'AT1G75080(+)', 'AT1G71250(+)', 'AT1G21970(+)', 'AT1G08560(+)', 
            'AT1G30690(+)', 'AT1G48930(+)', 'AT3G62830(+)', 'AT5G13170(+)', 'AT3G08900(+)',
            'AT1G68560(+)', 'AT2G47650(+)', 'AT2G25540(+)', 'AT4G33600(+)', 'AT3G19820(+)',
            'AT4G04460(+)', 'AT1G62340(+)', 'AT5G50750(+)']

genes_10 = ['AT3G26744(+)', 'AT1G75080(+)']
data_subset = adata.obsm['auc_mtx'][genes_10]


values_fixed = values.copy()
values_fixed[values_fixed < 0] = 0
values_shift = values_fixed - values_fixed.min() + 1e-6  # 防止log(0)
values_log = np.log(values_shift)

# 对log后的数据再做线性缩放到0-15范围
values = (values_log - values_log.min()) / (values_log.max() - values_log.min()) * 10


cell_types = adata.obs['Cell_type']

sorted_idx = np.argsort(cell_types.values)
data_sorted = data_subset.iloc[sorted_idx]
cell_types_sorted = cell_types.iloc[sorted_idx]

# 转成字符串
cell_types_sorted = cell_types_sorted.astype(str)

unique_types = cell_types_sorted.unique()
palette = sns.color_palette("tab20", len(unique_types))
lut = dict(zip(unique_types, palette))
row_colors = cell_types_sorted.map(lut)

import seaborn as sns
import matplotlib.pyplot as plt

sns.clustermap(data_sorted,
               row_cluster=True,
               col_cluster=True,
               row_colors=row_colors,
               cmap='viridis',
               standard_scale=1,
               yticklabels=False)

plt.savefig('auc_mtx_10genes_heatmap_by_Cell_type_clust.pdf')
plt.show()





import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

genes_10 = ['AT3G26744(+)', 'AT1G75080(+)']

x = adata.obs['x']
y = adata.obs['y']
cell_types = adata.obs['Cell_type']

# -------------------------------
# 1. 单独基因空间分布图，逐个保存pdf
for gene in genes_10:
    values = adata.obsm['auc_mtx'][gene].copy()
    
    # 负值设为0
    values[values < 0] = 0
    
    # 防止log(0)做平移
    values_shift = values - values.min() + 1e-6
    
    # log转换
    values_log = np.log(values_shift)
    
    # 线性缩放到0-10范围
    values_scaled = (values_log - values_log.min()) / (values_log.max() - values_log.min()) * 10
    
    plt.figure(figsize=(8,6))
    sc = plt.scatter(x, y, c=values_scaled, cmap='coolwarm', s=0.1)
    plt.colorbar(sc, label=f'Scaled activity of {gene}')
    plt.title(f'Spatial distribution of {gene}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(f'{gene}_spatial_activity.pdf')
    plt.close()

# -------------------------------
# 2. 聚类热图，两个基因活性，按细胞类型排序

data_subset = adata.obsm['auc_mtx'][genes_10].copy()

# 对矩阵整体做负值置0 + log + 缩放
data_subset[data_subset < 0] = 0
min_val = data_subset.min().min()
data_shifted = data_subset - min_val + 1e-6
data_log = np.log(data_shifted)
min_log = data_log.min().min()
max_log = data_log.max().max()
data_scaled = (data_log - min_log) / (max_log - min_log) * 10

# 根据细胞类型排序
sorted_idx = np.argsort(cell_types.values)
data_sorted = data_scaled.iloc[sorted_idx]
cell_types_sorted = cell_types.iloc[sorted_idx].astype(str)

# 颜色映射
unique_types = cell_types_sorted.unique()
palette = sns.color_palette("tab20", len(unique_types))
lut = dict(zip(unique_types, palette))
row_colors = cell_types_sorted.map(lut)

# 画聚类热图
g = sns.clustermap(data_sorted,
                   row_cluster=False,
                   col_cluster=True,
                   row_colors=row_colors,
                   cmap='viridis',
                   standard_scale=None,
                   yticklabels=False,
                   figsize=(8, 12))

g.fig.suptitle('Clustered regulon activities by Cell type', y=1.02)
g.fig.savefig('auc_mtx_10genes_heatmap_by_Cell_type.pdf')
plt.close()
















