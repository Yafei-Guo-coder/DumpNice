library(ggplot2)
library(RColorBrewer)
display.brewer.all()
col <- brewer.pal(n = 8, name = "Blues")[c(4,7)]
col <- brewer.pal(n = 8, name = "Oranges")[c(4,7)]
col <- brewer.pal(n = 8, name = "Greens")[c(4,7)]
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Go/V3/clusterProfiler/GOMAP")

#画整体的Go富集图----
path <- "/Users/guoyafei/Documents/01_Migration/02_Environment/02_XP-CLR/Go/V3/clusterProfiler/GOMAP/lineage/B"
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
    #scale_fill_gradient(expression(p.adjust),low="#9ECAE1", high = "#2171B5") +
    #B
    scale_fill_gradient(expression(p.adjust),low="#FDAE6B", high = "#D94801") +
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

