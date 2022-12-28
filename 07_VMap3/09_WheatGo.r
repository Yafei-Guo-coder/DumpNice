#wheatGo
library(ggplot2)
library(RColorBrewer)
display.brewer.all()
col <- brewer.pal(n = 8, name = "Blues")[c(4,7)]
col <- brewer.pal(n = 8, name = "Oranges")[c(4,7)]
col <- brewer.pal(n = 8, name = "Greens")[c(4,7)]
col <- brewer.pal(n = 8, name = "Red")[c(4,7)]

setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/shufgene")
setwd("/Users/guoyafei/Documents/02_VmapIII/08_Network/PPI")
data <- read.table("go.txt",header=T,stringsAsFactors = F,sep="\t")
colnames(data)[c(3,4,6)] <- c("Description", "Count", "p.adjust")
ggplot(data=data)+
  geom_bar(aes(x= Description, y=Count, fill=(-log10(p.adjust))), stat='identity') +
  coord_flip() +
  scale_fill_gradient(expression(-log(p.adjust)),low="#FDAE6B", high = "#D94801") +
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
  ylim(0,130)+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))


#画整体的Go富集图----
path <- "/Users/guoyafei/Documents/02_VmapIII/08_Network/shufgene/Go"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=F, sep = "\t", stringsAsFactors=F)}) 
plot_list = list() 
for(i in 1:length(data)){
  colnames(data[[i]])[c(3,4,6)] <- c("Description", "Count", "p.adjust")
  p <- ggplot(data=data[[i]])+
    geom_bar(aes(x= Description, y=Count, fill=(-log(p.adjust))), stat='identity') +
    coord_flip() +
    #facet_grid(.~Name,scales="free") +
    #facet_wrap(~Name,ncol = 1)+
    #A
    scale_fill_gradient(expression(-log(p.adjust)),low="#FDAE6B", high = "#D94801") +
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
    )+
    theme_bw()+
    theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
  
  plot_list[[i]] = p
}

name <- strsplit(names(data), ".txt")
for (i in 1:length(data)) { 
  file_name = paste(name[[i]],".pdf",sep="") 
  pdf(file_name,height = 12,width = 9) 
  print(plot_list[[i]]) 
  dev.off() 
} 
