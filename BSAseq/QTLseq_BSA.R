setwd("C:/Users/Administrator/Desktop")
library(QTLseqr)
#��ȡGATK VariantsToTable������.table����
df <- importFromGATK(file = "BSA.filter.table",
                     highBulk = "NHR",
                     lowBulk = "HR",
                     chromList = 1:5)

#��Ե�������ݽ���SNP���ˣ�
#���Ʋ������ֱ��ͼ���鿴���������reads����ȷֲ�
library("ggplot2")
ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
  xlim(0,1000)

#�鿴�ο��������ϵĵ�λ����Ƶ�ʷֲ������F2Ⱥ��ṹ
ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))

#�鿴�ߵͻ�ص�SNPindex���
ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.HIGH))

ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.LOW))


#��ʽ��ʼ����
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


#����F2Ⱥ��,����deltaSNP
df_filt <- runQTLseqAnalysis(df_filt,
                             windowSize = 1e6,
                             popStruc = "F2",
                             bulkSize = c(25, 25),
                             replications = 10000,
                             intervals = c(95, 99) )
df_filt <- subset(df_filt, !is.na(SNPindex.LOW) | !is.na(SNPindex.HIGH))


#����SNP����ֵ��Ⱦɫ�������ͼ
p1 <- plotQTLStats(SNPset = df_filt, var = "nSNPs")
p1
#����deltaSNP��Ⱦɫ�������ͼ
p2 <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
p2
#����deltaSNP��ĳ��Ⱦɫ�壨1�ź�3��Ⱦɫ�壩��������ͼ
p4 <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE,line=FALSE)
p4


#����Gֵ
df_filt <- runGprimeAnalysis(df_filt,
                             windowSize = 1e6,
                             outlierFilter = "deltaSNP",
                             filterThreshold = 0.1)
df_filt <- subset(df_filt, !is.na(SNPindex.LOW) | !is.na(SNPindex.HIGH))
head(df_filt)

#�鿴G���ֲ���count�����ķֲ�
plotGprimeDist(SNPset =df_filt, outlierFilter = "deltaSNP", filterThreshold = 0.1)


#����G'��Ⱦɫ�������ͼ
p3 <- plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
p3



#��ȡ��ȫ��������QTLλ�㣬�������Ϊcsv�ļ�
results <- getQTLTable(SNPset = df_filt, method = "Gprime",alpha = 0.01, export = FALSE)
results








#�����������������в���֮��õ���ͼ���ٸ��ݴ�ͼ������λ�ã����Ⱦɫ��Ƭ�Ρ�
#������һ��Ⱦɫ��Ϊ����Ҫ���һ��Ⱦɫ��deltaSNP������λ��
#����������ص�SNPindexɢ��ͼ
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

#����F2Ⱥ��
df_filt2 <- runQTLseqAnalysis(df_filt2,
                              windowSize = 1e6,
                              popStruc = "F2",
                              bulkSize = c(385, 430),
                              replications = 10000,
                              intervals = c(95, 99) )

#����deltaSNP��Gֵ
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