install.packages("gdata")
library(cinaR)
library(ggplot2)
library(RColorBrewer)
library(aplot)
library(grid)
library(gridExtra)
display.brewer.all()
col <- brewer.pal(n = 8, name = "Spectral")
setwd("/Users/guoyafei/Documents/02_VmapIII/04_Statistics/04_Structure")
eigvec<- read.table("/Users/guoyafei/Documents/02_VmapIII/04_Statistics/04_Structure/ABlineage.maf0.05.5k.eigenvec",header=T,stringsAsFactors = F)
eigval <- read.table("/Users/guoyafei/Documents/02_VmapIII/04_Statistics/04_Structure/ABlineage.maf0.05.5k.eigenval",header=F,stringsAsFactors = F)
mds <- read.table("/Users/guoyafei/Documents/02_VmapIII/04_Statistics/04_Structure/ABlineage.maf0.05.60k.mds",header=T,stringsAsFactors = F)
popinfo <- "/Users/guoyafei/Documents/02_VmapIII/01_表格/V.2/Lulab_germplasm_Info.xlsx"
poptable <- read.xls(popinfo, sheet=1, verbose=FALSE, fill=NA)
popAB <- poptable[which(poptable$GenomeNAtype!="DD"),c(1,2,3,12,13,14,18,20)]
colnames(popAB)[1] <- "FID"
all <- merge(mds,popAB,by="FID")
#PC1:19.03% PC2:11.18%
#PC1:19.03% PC2:11.18%
#colnames(popAB) <- c("order","exp_id","vcf_id","group","color","pch","group2")
#source("/Users/guoyafei/Documents/01_个人项目/04_VmapIII/16_分群/pca.plot2d.r")
#pca_plot(eigenvector = eigvec, eigenvalue = eigval,
#         group = popinfo, key = key, outdir = od,
#         shape = T, shapes = pop$group2, border = T, border_size = 2.5,
#         line0 = T, line0_size = 1)
p<-ggplot(all,aes(x=C1,y=C2,shape=as.factor(all$PCA_shape),color=as.factor(all$PCA_color)))
p<-p+geom_point()+
  scale_color_brewer(palette = "Accent")+
  scale_shape_manual(values = c(1:11))+
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
  labs(x = 'PCA 1 (19.03%)', y = 'PCA 2 (11.18%)') +  
  theme_bw() +
  theme(
    #legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.title.y=element_text(size = 15),
    axis.title.x=element_text(size = 15),
)
site <- read.table("/Users/guoyafei/Documents/02_VmapIII/04_Statistics/04_Structure/AB.maf0.05_siteQCfile.txt",header=T,stringsAsFactors = F)
taxa <- read.table("/Users/guoyafei/Documents/02_VmapIII/04_Statistics/04_Structure/AB.maf0.05_taxaQCfile.txt",header=T,stringsAsFactors = F)
p0 <- ggplot(data=site,aes(x=MissingRate)) + 
  geom_histogram(binwidth=0.01,fill="#69b3a2", 
                 color="#e9ecef", alpha=0.9)+
  theme_bw()+
  labs(x="MissingRate",y="SNP Count")+
  xlim(0,0.2)+
  theme(
  #legend.position="none",
  #panel.border = element_blank(),
  axis.line.y = element_line(),
  axis.line.x = element_line(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.text.x=element_text(size=15),
  axis.text.y=element_text(size=15),
  axis.title.y=element_text(size = 15),
  axis.title.x=element_text(size = 15),

)
p0
