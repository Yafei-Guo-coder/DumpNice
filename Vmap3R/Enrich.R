library(ggplot2)
library(RColorBrewer)
display.brewer.all()
col1 <- brewer.pal(n = 8, name = "YlOrBr")[c(2,5)]
col2 <- brewer.pal(n = 8, name = "YlGnBu")[c(2,5)]
col3 <- brewer.pal(n = 8, name = "YlGn")[c(2,5)]
col4 <- brewer.pal(n = 8, name = "Purples")[c(2,5)]
col5 <- brewer.pal(n = 8, name = "Greys")[c(2,5)]

#setwd("/Users/guoyafei/Documents/Lulab/Project-2-Migration/anno_gene")
#nameA <- c("Domesticated_einkorn", "Domesticated_emmer","Durum","EA","EU","Indian_dwarf","Khorasan_wheat","Landrace","Macha","Persian_wheat","Polish_wheat","Rivet_wheat","SCA","Spelt","Tibetan_semi_wild","Urartu","Vavilovii","WA","Wild_Einkorn","Wild_emmer","Xinjiang_wheat")
#nameB <- c("Domesticated_emmer","Durum","EA","EU","Khorasan_wheat", "Landrace","Macha","Persian_wheat","Polish_wheat","Rivet_wheat","SCA","Spelt","Speltoides","Vavilovii","WA","Wild_emmer","Xinjiang_wheat","Yunan_wheat")
#nameD <- c("Anathera","Club_wheat","EA","EU","Indian_dwarf","Landrace","Macha","Meyeri","SCA","Spelt","Strangulata","Tibetan_semi_wild","Vavilovii","WA","Xinjiang_wheat","Yunan_wheat")

path <- "//Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/XP-CLR/Go/TXT"
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x, header=T, sep = "\t", stringsAsFactors=F)}) 
#for(i in 1:length(data)){
#  data[[i]]$Name <- nameB[i]
#}
#all <- data.frame()
#for(i in 1:length(data)){
#  all <- rbind(all, data[[i]])
#}
#plot_list = list() 
#for(i in 1:16){
p <- ggplot(data=data[[6]])+
  geom_bar(aes(x= Description, y=Number.in.input.list, fill=(FDR)), stat='identity') +
  coord_flip() +
  #facet_grid(.~Name,scales="free") +
  #facet_wrap(~Name,ncol = 1)+
  scale_fill_gradient(expression(FDR),low=col1[2], high = col1[1]) +
  ylab("Gene number") +
  xlab("GO term description") +
  #expand_limits(y=c(0,8))+
  #ylim(0,100)+
  theme(
    axis.text.x=element_text(color="black",size=rel(0.8)),
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
#plot_list[[i]] = p
#}
for (i in 1:16) { 
  file_name = paste("D",i,".pdf",sep="") 
  pdf(file_name,height = 12,width = 9) 
  print(plot_list[[i]]) 
  dev.off() 
} 

