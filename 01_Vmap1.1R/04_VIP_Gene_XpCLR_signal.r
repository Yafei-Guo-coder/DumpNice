library(qqman)
#画信号位点基因的XP-CLR的染色体信号
#Working directory
#xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%/VIP_gene_XpCLR
#VIP_gene
#1. TraesCS4D02G364400	Vrn2-2	24	58277633	58279728	chr4D:509282253-509284348(-)	Transcription factor GHD7 [UniProtKB/Swiss-Prot:E5RQA1]	Os10g0560400	NA
#2. TraesCS2D02G079600	Ppd-1(PRR)	11	33952048	33956269	chr2D:33952048-33956269(-)	Two-component response regulator-like PRR37 [UniProtKB/Swiss-Prot:A2YQ93]	Os07g0695100	NA
#3. TraesCS2A02G081900	Ppd-A1	7	36933684	36938202	chr2A:36933684-36938202(-)	Two-component response regulator-like PRR37 [UniProtKB/Swiss-Prot:A2YQ93]	Os07g0695100	NA
#4. TraesCS5A02G473800	Q-5A	26	196896739	196900381	chr5A:650127258-650130900(-)	APETALA2-like protein 2 [UniProtKB/Swiss-Prot:Q84TB5]	Os03g0818800	NA
#5. TraesCS4B02G043100	Rht-B1	21	30861268	30863723	chr4B:30861268-30863723(+)	DELLA protein RHT-1 [UniProtKB/Swiss-Prot:Q9ST59]	Os03g0707600	AT1G66350
#6. TraesCS1D02G029100	Sr33	5	11451423	11459353	chr1D:11451423-11459353(+)	Disease resistance protein RGA5 [UniProtKB/Swiss-Prot:F7J0N2]	NA	AT3G46530
#7. TraesCS1D02G040400	Sr45	5	19341296	19346065	chr1D:19341296-19346065(+)	Putative disease resistance protein RGA4 [UniProtKB/Swiss-Prot:Q7XA39]	Os10g0130800	N

#File Gene_Name Chr Start Stop
#North2_South_smooth_A Ppd-A1 7	36933684	36938202
#North2_South_smooth_B Rht-B1 21	30861268	30863723
#North2_South_smooth_D Sr45 5	19341296	19346065
#WA_EU_smooth_A  Q-5A 26	196896739	196900381
#WA_EU_smooth_D  Sr45 5	19341296	19346065
#EU_South_smooth_B Rht-B1 21	30861268	30863723
#EU_South_smooth_D Flowering1 24	58277633	58279728
#EU_South_smooth_D Sr45 5	19341296	19346065
#WA_South_smooth_D Flowering1 24	58277633	58279728
#WA_South_smooth_D Sr33 5	11451423	11459353
#Tibet_South_smooth_D  Ppd-1(PRR) 11	33952048	33956269
cat $1 |while read file chr1 chr2
do
awk '{if($1== "'${chr1}'" || $1=="'${chr2}'") {print $0}}' ../../${file}| sort -k1,1n -k2,2g |sed '1i  Chr\tWindowStart\tWindowStop\tSNPcount\tMeanY\tWstat' > ${file::-4}.${chr1}.${chr2}.txt
done
#把42条染色体合并成21条
#转移到本地画图
#/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Xpclr/V3
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Xpclr/V3")
#批量读取Xp-clr结果文件
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Xpclr/V3/TXT" ##文件目录
fileNames <- dir(path)  ##获取该路径下的文件名
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   ##生成读取文件路径
data <- lapply(filePath, function(x){
  read.table(x, header=T,stringsAsFactors = F)})
#读取阈值文件
thresHold <- read.table("thresHold.txt",header=T,stringsAsFactors = F)
rownames(thresHold) <- thresHold$GEO_region
#读取着丝粒文件
centromer <- read.table("centromerer.txt",header=T,stringsAsFactors = F)
rownames(centromer) <- centromer$chr
#读取VIP基因文件
gene <- read.table("VIP_gene.txt",header=T,stringsAsFactors = F)
rownames(gene) <- paste(gene$File,gene$Position,sep=".")
out1 <- strsplit(sub('.chr',':chr', names(data)), ":")
out2 <- strsplit(sub('.txt',':txt', names(data)), ":")
pdf("Xp-clr.pdf",width = 10,height = 5)
for(i in c(1:length(data))){
  name <- out1[[i]][1]
  tit <- out2[[i]][1]
  chr <- strsplit(out[[i]][2], ".txt")[[1]][1]
  title <- paste(gene[tit,1],gene[tit,2],sep=":")
  cen <- centromer[chr,4]
  plot(data[[i]]$WindowStart/1000000, data[[i]]$MeanY, ylim = c(-2,50),axes = T,col = "skyblue",type="l",cex = 2,lwd = 3,ann = F)
  abline(v = gene[tit,7]/1000000, col = "red", cex=2,lwd = 1)
  abline(h=thresHold[name,3], col = "lightgrey", lty = 3,cex=2,lwd = 4)
  title(title,xlab= chr, ylab = 'Xp-clr', font.lab = 1, cex.lab = 1.5)
  points(cen/1000000,0, pch=20,cex=2,col=col)
}
dev.off()
