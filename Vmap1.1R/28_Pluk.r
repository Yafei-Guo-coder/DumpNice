library(vegan)
setwd("/Users/guoyafei/Downloads/200328-普鲁克分析（Procrustes Analysis）评估物种-环境或功能关联度的一个示例")
#PCA降维及探索性分析
##样方-环境属性矩阵
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#环境变量的 PCA 需要标准化，详情 ?rda
env_pca <- rda(env, scale = TRUE)

##样方-物种丰度矩阵
otu <- read.delim('spe_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#物种数据 Hellinger 预转化
otu_hel <- decostand(otu, method = 'hellinger')

#对转化后的物种数据执行 PCA，无需标准化
otu_pca <- rda(otu_hel, scale = FALSE)

##排序图比较，以 PCA 的 I 型标尺为例
par(mfrow = c(1, 2))
biplot(env_pca, choices = c(1, 2), scaling = 1, 
       main = '环境组成的PCA', col = c('red', 'blue'))
biplot(otu_pca, choices = c(1, 2), scaling = 1, 
       main = '物种组成的PCA', col = c('red', 'blue'))
#Procrustes分析度量两数据集一致性
#Procrustes 分析
#提取两个 PCA 中的样方排序坐标，均以 I 型标尺为例
site_env <- summary(otu_pca, scaling = 1)$site
site_otu <- summary(otu_pca, scaling = 1)$site

#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(X = env_pca, Y = otu_pca, symmetric = TRUE)
summary(proc)

#旋转图
plot(proc, kind = 1, type = 'text')

#一些重要的结果提取
names(proc)
head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置
#残差图
plot(proc, kind = 2)
residuals(proc)  #残差值
#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(X = env_pca, Y = otu_pca, permutations = how(nperm = 999))
prot

#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量
library(ggplot2)

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')
#----------------------------------------------------
#ggplot2 作图
p <- ggplot(Y) +
  geom_point(aes(X1, X2, color = groups), size = 1.5, shape = 16) +
  geom_point(aes(PC1, PC2, color = groups), size = 1.5, shape = 1) +
  scale_color_manual(values = c('red2', 'purple2', 'green3'), limits = c('A', 'B', 'C')) +
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.1, 'cm')),
               color = 'blue', size = 0.3) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = sprintf('M^2 == 0.2178'),
           x = -0.21, y = 0.42, size = 3, parse = TRUE) +
  annotate('text', label = 'P < 0.001',
           x = -0.21, y = 0.38, size = 3, parse = TRUE)

p

#输出图片
ggsave('procrustes.pdf', p, width = 6, height = 5)
ggsave('procrustes.png', p, width = 6, height = 5)