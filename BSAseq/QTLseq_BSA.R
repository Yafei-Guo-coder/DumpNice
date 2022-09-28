setwd("C:/Users/Administrator/Desktop")
library(QTLseqr)
#读取GATK VariantsToTable导出的.table数据
df <- importFromGATK(file = "BSA.filter.table",
                     highBulk = "NHR",
                     lowBulk = "HR",
                     chromList = 1:5)

#针对导入的数据进行SNP过滤：
#绘制测序深度直方图：查看两个混池中reads的深度分布
library("ggplot2")
ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
  xlim(0,1000)

#查看参考基因组上的等位基因频率分布，检查F2群体结构
ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))

#查看高低混池的SNPindex情况
ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.HIGH))

ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.LOW))


#正式开始过滤
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 40,
    maxTotalDepth = 200,
    minSampleDepth = 20,
    depthDifference = 100,
    minGQ = 95,
    verbose = TRUE
  )


#分析F2群体,计算deltaSNP
df_filt <- runQTLseqAnalysis(df_filt,
                             windowSize = 1e6,
                             popStruc = "F2",
                             bulkSize = c(25, 25),
                             replications = 10000,
                             intervals = c(95, 99) )
df_filt <- subset(df_filt, !is.na(SNPindex.LOW) | !is.na(SNPindex.HIGH))


#绘制SNP数量值在染色体的折线图
p1 <- plotQTLStats(SNPset = df_filt, var = "nSNPs")
p1
#绘制deltaSNP在染色体的折线图
p2 <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
p2
#绘制deltaSNP在某条染色体（1号和3号染色体）的曼哈顿图
p4 <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE,line=FALSE)
p4


#计算G值
df_filt <- runGprimeAnalysis(df_filt,
                             windowSize = 1e6,
                             outlierFilter = "deltaSNP",
                             filterThreshold = 0.1)
df_filt <- subset(df_filt, !is.na(SNPindex.LOW) | !is.na(SNPindex.HIGH))
head(df_filt)

#查看G‘分布与count数量的分布
plotGprimeDist(SNPset =df_filt, outlierFilter = "deltaSNP", filterThreshold = 0.1)


#绘制G'在染色体的折线图
p3 <- plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
p3



#提取出全部显著的QTL位点，并且输出为csv文件
results <- getQTLTable(SNPset = df_filt, method = "Gprime",alpha = 0.01, export = FALSE)
results








#建议先跑完上述所有步骤之后得到大图。再根据大图标明的位置，拆解染色体片段。
#下面以一号染色体为例，要获得一号染色体deltaSNP显著的位置
#绘制两个混池的SNPindex散点图
library(ggplot2)
df2 <- importFromGATK(file = "BSA.filter.table",
                      highBulk = "NHR",
                      lowBulk = "HR",
                      chromList = paste0(rep("Chr", 1), 1))
df_filt2 <-filterSNPs(
  SNPset = df2,
  refAlleleFreq = 0.20,
  minTotalDepth = 100,
  maxTotalDepth = 400,
  depthDifference = 100,
  minSampleDepth = 40,
  minGQ = 99,
  verbose = TRUE
)

#分析F2群体
df_filt2 <- runQTLseqAnalysis(df_filt2,
                              windowSize = 1e6,
                              popStruc = "F2",
                              bulkSize = c(385, 430),
                              replications = 10000,
                              intervals = c(95, 99) )

#计算deltaSNP和G值
df_filt2 <- runGprimeAnalysis(df_filt2,
                              windowSize = 1e6,
                              outlierFilter = "deltaSNP",
                              filterThreshold = 0.1)
df_filt2 <- subset(df_filt2, !is.na(SNPindex.LOW) | !is.na(SNPindex.HIGH))


#write.csv(SNPindex.HIGH_all, "SNPindex.HIGH_all.csv")
SNPindex.HIGH_Chr1<- df_filt2[,grepl("CHROM|POS|SNPindex.HIGH", colnames(df_filt2))]
ggplot(SNPindex.HIGH_Chr1, aes(x=POS/1000000, y=SNPindex.HIGH),axes = T)+
  ylim(0,1)+
  geom_point(size=1.5, color="#0000FF")+
  labs(title="HIGH-bulk",x="Chromosome 1 position(Mb)", y = "SNP-index")+
  geom_hline(aes(yintercept=0.5),linetype="dashed")+
  theme_classic()
