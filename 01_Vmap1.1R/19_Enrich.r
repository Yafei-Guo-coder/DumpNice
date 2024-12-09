#Go analysis
#Working directory
#xuebo@204:/data2/xuebo/Projects/Speciation/xpclr/Selection_V3/smooth/lineage_V2/Top5%/Go
#以下invalid
for i in A B D
do
for j in `cat names`
do
grep -v -f neg_WA_North2_smooth_${i}.top5.txt ${j}_smooth_${i}.top5.txt > ${j}.go.gene.txt
done
done

#由于clusterProfiler安装的原因，转移到66上继续
#66@yafei:/data1/home/yafei/008_Software/wheatGO-v1.1
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/North2_South_smooth_A.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/Tibet_South_smooth_A.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/WA_EU_smooth_A.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/WA_South_smooth_A.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All

./wheatGO-v1.1-Ontologizer -g All_VIP_gene/North2_South_smooth_B.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/Tibet_South_smooth_B.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/WA_EU_smooth_B.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/WA_South_smooth_B.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All

./wheatGO-v1.1-Ontologizer -g All_VIP_gene/North2_South_smooth_D.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/Strang_WA_smooth_D.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/WA_South_smooth_D.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/Tibet_South_smooth_D.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All
./wheatGO-v1.1-Ontologizer -g All_VIP_gene/WA_EU_smooth_D.top5.txt -m GOSLIM -c Parent-Child-Intersection -p Benjamini-Hochberg -r 100 -s All

#文件目录
#ABD_GOMAP
#ABD_GOMAP_All
#ABD_GOSLIM
#ABD_GOSLIM_All
#A_B_D_GOMAP
#A_B_D_GOMAP_All
#A_B_D_GOSLIM
#A_B_D_GOSLIM_All

for i in `cat all_names`
do
sort -k11,11g A_B_D_GOSLIM_All/table-${i}.top5-Parent-Child-Intersection-Benjamini-Hochberg.txt | awk '{if($11<0.05) {print $0}}' | sed '1i ID\tPop.total\tPop.term\tStudy.total\tStudy.term\tPop.family\tStudy.family\tnparents\tis.trivial\tp\tp.adjusted\tp.min\tname' > A_B_D_GOSLIM_All/${i}.go.txt   
done


library(ggplot2)
library(RColorBrewer)
display.brewer.all()
col <- brewer.pal(n = 8, name = "Blues")[c(4,7)]
col <- brewer.pal(n = 8, name = "Oranges")[c(4,7)]
col <- brewer.pal(n = 8, name = "Greens")[c(4,7)]
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Go/V3/GOMAP")


#画整体的Go富集图----
#path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Go/V3/GOMAP/D"
setwd("/Users/guoyafei/Desktop/baypass/GO")
lat <- read.table("type4.go.txt", header=T,stringsAsFactors = F,sep="\t")
ggplot(data=lat)+
  geom_bar(aes(x=Description, y=number, fill=(p.adjust)), stat='identity') +
  coord_flip() +
  scale_fill_gradient(expression(p.adjust),low="#FDAE6B", high = "#D94801") +
  ylab("Gene number") +
  xlab("GO term description") +
  theme(
    axis.text.x=element_text(color="black",size=rel(0.3)),
    axis.text.y=element_text(color="black", size=rel(0.3)),
    axis.title.x = element_text(color="black", size=rel(5)),
    axis.title.y = element_blank(),
    legend.text=element_text(color="black",size=rel(0.2)),
    legend.title = element_text(color="black",size=rel(0.7))
  )+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))




fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T, sep = "\t", stringsAsFactors=F)}) 
plot_list = list() 
for(i in 1:length(data)){
  p <- ggplot(data=data[[i]])+
    geom_bar(aes(x= Description, y=Count, fill=(p.adjust)), stat='identity') +
    coord_flip() +
    #facet_grid(.~Name,scales="free") +
    #facet_wrap(~Name,ncol = 1)+
    #A
    scale_fill_gradient(expression(p.adjust),low="#FDAE6B", high = "#D94801") +
    #B
    #scale_fill_gradient(expression(p.adjust),low="#9ECAE1", high = "#2171B5") +
    #D
    #scale_fill_gradient(expression(p.adjust),low="#A1D99B", high = "#238B45") +
    ylab("Gene number") +
    xlab("GO term description") +
    #expand_limits(y=c(0,8))+
    #ylim(0,100)+
    theme(
      axis.text.x=element_text(color="black",size=rel(0.3)),
      axis.text.y=element_text(color="black", size=rel(0.3)),
      axis.title.x = element_text(color="black", size=rel(5)),
      axis.title.y = element_blank(),
      legend.text=element_text(color="black",size=rel(0.2)),
      legend.title = element_text(color="black",size=rel(0.7))
      #legend.position=c(0,1),legend.justification=c(-1,0)
      #legend.position="top",
    )+
    #  scale_x_discrete(limits= c("Domesticated_einkorn", "Domesticated_emmer","Durum","EA","EU","Indian_dwarf","Khorasan_wheat","Landrace","Macha","Persian_wheat","Polish_wheat","Rivet_wheat","SCA","Spelt","Tibetan_semi_wild","Urartu","Vavilovii","WA","Wild_Einkorn","Wild_emmer","Xinjiang_wheat"))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
  
  plot_list[[i]] = p
}

name <- strsplit(names(data), ".txt")
for (i in 1:length(data)) { 
  file_name = paste(name[[i]],".pdf",sep="") 
  pdf(file_name,height = 12,width = 9) 
  print(plot_list[[i]]) 
  dev.off() 
} 
