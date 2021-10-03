library(CMplot)
setwd("/Users/guoyafei/Documents/个人项目/Project-4-VmapIII/Fastcall2/测试数据/")

mydata<-read.table("/Users/guoyafei/Documents/个人项目/Project-4-VmapIII/Fastcall2/测试数据/fastcall2_001_pos.txt",header=TRUE,sep="\t")
head(mydata)
# snp         chr       pos
# snp1_1    1        2041
# snp1_2    1        2062
# snp1_3    1        2190
CMplot(mydata,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300) 
mydata<-read.table("/Users/guoyafei/Documents/个人项目/Project-4-VmapIII/Fastcall2/测试数据/fastcall_001_pos.txt",header=TRUE,sep="\t")
head(mydata)
# snp         chr       pos
# snp1_1    1        2041
# snp1_2    1        2062
# snp1_3    1        2190
CMplot(mydata,plot.type="d",bin.size=1e4,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300) 
