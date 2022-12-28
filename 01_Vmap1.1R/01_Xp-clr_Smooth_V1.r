#XP-CLR
#工作路径：xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V2
#XPCLR的流程：XPCLR -xpclr ../groupSouth/groupSouthChr${i}.geno ../groupSCA/groupSCAChr${i}.geno ../groupSCA/groupSCAChr${i}.snp South_SCA_10kchr${i} -w1 0.005 500 10000 $i -p1 0.95 &
#这是参数 -w1 0.005 500 10000 -p1 0.95
#下面的是做smooth标准化
#该目录下一共有16个结果文件夹，分别是
#XPCLRresult_CA_NW_A
#XPCLRresult_CA_SA
#XPCLRresult_NW_A_NE_A
#XPCLRresult_SA_SW_A
#XPCLRresult_SA_Tibet
#XPCLRresult_SE_A_NE_A
#XPCLRresult_SW_A_NE_A
#XPCLRresult_Strang_WA
#XPCLRresult_WA_CA
#XPCLRresult_WA_EU
#XPCLRresult_WA_NE_A
#XPCLRresult_WA_NW_A
#XPCLRresult_WA_SA
#XPCLRresult_WA_SE_A
#XPCLRresult_WA_SW_A
#XPCLRresult_WA_Tibet
#!/usr/bin/Rscript.R
library(GenWin)
library(dplyr)
Folder <- c("CA_NW_A","CA_SA","NW_A_NE_A","SA_SW_A","SA_Tibet","SE_A_NE_A","SW_A_NE_A","Strang_WA","WA_CA","WA_EU","WA_NE_A","WA_NW_A","WA_SA","WA_SE_A","WA_SW_A","WA_Tibet")
for (num in c(1:16)){
  filePath <- paste("/data2/xuebo/Projects/Speciation/xpclr/Selection_V2/XPCLRresult_",Folder[num],sep="")
  setwd(filePath)
  for(i in c(1:42)){
    file=paste(Folder[num],"_10kchr",i,".xpclr.txt",sep="")
    Data1=read.table(file,sep=" ")
    Data2<-na.omit(Data1)
    Data <- filter(Data2, Data2[,6] != Inf)
    Z=matrix(, nrow = nrow(Data),ncol=1)
    figure=paste("pdfnom_c",i,".pdf",sep="")
    pdf(figure,width=12,height=8)
    for(j in 1:nrow(Data)){
      Z[j,1]=(Data[j,6]-mean(Data$V6))/sd(Data$V6)
    }
    NORM=splineAnalyze(Y=Z,map=Data$V4,smoothness = 2000,plotRaw=T,plotWindows = T,method = 4)
    dev.off()
    normScore=NORM$windowData
    outFile=paste("/data2/xuebo/Projects/Speciation/xpclr/Selection_V2/smooth/",Folder[num],"_smooth",i,".txt",sep="")
    write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
  }
}
#Rscript smooth.R

#GenWin(smooth)的结果是这样的"WindowStart" "WindowStop" "SNPcount" "MeanY" "Wstat"，没有染色体，要加上染色体的信息
#工作路径：xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/North_South_SCA/XPCLRresultSouth_SCA/smooth
#接下来给smooth结果添加染色体号并进行排序，并且合并成A，B，D lineage.

library(GenWin)
library(dplyr)
Folder <- c("SH_EA", "SH_IA")
for (num in c(1:9)){
  filePath <- paste("/data1/home/yafei/004_Vmap3/xpclr/ResultXPCLR/",Folder[num],"/",sep="")
  setwd(filePath)
  for(i in c(1:9)){
    file=paste(Folder[num],"_10kchr00",i,".xpclr.txt",sep="")
    Data1=read.table(file,sep=" ")
    Data2<-na.omit(Data1)
    Data <- filter(Data2, Data2[,6] != Inf)
    Z=matrix(, nrow = nrow(Data),ncol=1)
    figure=paste("pdfnom_c",i,".pdf",sep="")
    pdf(figure,width=12,height=8)
    m <- mean(Data$V6)
    s <- sd(Data$V6)
    for(j in 1:nrow(Data)){
      Z[j,1]=(Data[j,6]-m)/s
    }
    NORM=splineAnalyze(Y=Z,map=Data$V4,smoothness = 2000,plotRaw=T,plotWindows = T,method = 4)
    dev.off()
    normScore=NORM$windowData
    outFile=paste("/data1/home/yafei/004_Vmap3/xpclr/ResultXPCLR/",Folder[num],"/smooth/",Folder[num],"_smooth",i,".txt",sep="")
    write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
  }
  for(i in c(10:42)){
    file=paste(Folder[num],"_10kchr0",i,".xpclr.txt",sep="")
    Data1=read.table(file,sep=" ")
    Data2<-na.omit(Data1)
    Data <- filter(Data2, Data2[,6] != Inf)
    Z=matrix(, nrow = nrow(Data),ncol=1)
    figure=paste("pdfnom_c",i,".pdf",sep="")
    pdf(figure,width=12,height=8)
    m <- mean(Data$V6)
    s <- sd(Data$V6)
    for(j in 1:nrow(Data)){
      Z[j,1]=(Data[j,6]-m)/s
    }
    NORM=splineAnalyze(Y=Z,map=Data$V4,smoothness = 2000,plotRaw=T,plotWindows = T,method = 4)
    dev.off()
    normScore=NORM$windowData
    outFile=paste("/data1/home/yafei/004_Vmap3/xpclr/ResultXPCLR/",Folder[num],"/smooth/",Folder[num],"_smooth",i,".txt",sep="")
    write.table(normScore,outFile,sep="\t",col.names = T,row.names = F)
  }
}
}
