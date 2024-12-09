library(RColorBrewer)
#color <- brewer.pal(5, "RdBu")[c(2,5)]
#color <- append(color,brewer.pal(5, "PuOr")[c(2,5)])
#brewer.pal(5, "PiYG")[c(2,5)]
#brewer.pal(11, "RdYlBu")[1]
color <- c(brewer.pal(5, "RdBu")[c(2,5)],brewer.pal(5, "PuOr")[c(2,5)],brewer.pal(5, "PiYG")[c(2,5)],brewer.pal(11, "RdYlBu")[c(10,11)])
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/Ne")

setwd("~/Desktop/")
####改格式
library(RColorBrewer)
library(cowplot)
library(viridis)
display.brewer.all()
brewer.pal(12,'Set3')
brewer.pal(7,'PuRd')
brewer.pal(7,'Blues')
brewer.pal(7,'OrRd')
brewer.pal(7,'Greens')
a = seq(1:10)
#cols = viridis(8)
cols <- color
plot(NULL,xlim = c(2,3),ylim = c(0,8),axes=F,xlab="",ylab="")
#set.seed(2)
#rom = sample(1:14,14,replace = F)
j = 0
for (i in (1:8)){
  j=j+1
  lines(a,rep(j,10),type="l",col=cols[i],lwd = 19)
}
#col1 = c("#66C2A5", "#FC8D62", "#8DA0CB","#FFD92F","#E78AC3","#07dede","#E5C494")
#col2 = c("#f15a24","#779970","#8dd3c7","#fcee21","#ba9bc9","#619cff","#c69c6d")
col1 = c("#f5e1fa", "#C6DBEF", "#FDD49E","#9FDA3AFF","#8dd3c7","#277F8EFF","#FC8D59")
#setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/smcpp/getpoints")
dataf = read.table("/Volumes/My Passport/01_VMap1.1/02_Environment/03_Ne/count2.txt",header = T)
unique(dataf$Subspecies)
datafd <- dataf[which(dataf$Subspecies != "Urartu" & dataf$Subspecies != "Speltoilds" & dataf$Time <=10000),]
#datafd$Subspecies = factor(datafd$Subspecies,levels=c("Urartu","Speltoilds","Strangulata","Wild_emmer",
#                                                      "Domesticated_emmer","Free_threshing","Bread_wheat"))

datafd$Subspecies = factor(datafd$Subspecies,levels=c("Strangulata","Wild_emmer",
                                                      "Domesticated_emmer","Free_threshing","Bread_wheat"))


#col1 <- c(brewer.pal(10, "RdBu")[c(3,4,5,7,8,9)],brewer.pal(8, "Set1")[6])

col1 <- c(brewer.pal(10, "RdBu")[c(3,4,7,8)],brewer.pal(8, "Set1")[6])

col1 <- c(brewer.pal(5, "RdBu")[5],brewer.pal(5, "PuOr")[5],brewer.pal(5, "PiYG")[5],brewer.pal(5, "RdBu")[2],brewer.pal(5, "PuOr")[2],brewer.pal(5, "PiYG")[2],brewer.pal(11, "RdYlBu")[11])

p1 = ggplot(datafd, aes(x = Time,y = Contribution,fill = Subspecies))+
  ####position="stack"堆叠状
  ####position="stack" 改为"postion_stack(reverse=T) 反向堆叠
  geom_bar(stat ="identity",width = 800,position ="stack")+
  scale_fill_manual(values = col1)+
  labs(x = "",y = "All effective population size",size = rel(1)) +
  #geom_text(aes(label = datafd$Num),position=position_stack(vjust =0.5),size = 5)+ 这行代码是显示个数的数字
  guides(fill = guide_legend(reverse = F))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,size = rel(0.8))) +
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))

ggplot(datafd, aes(x = Time,y = Contribution,fill = Subspecies))+
  ####position="stack"堆叠状
  ####position="stack" 改为"postion_stack(reverse=T) 反向堆叠
  geom_bar(stat ="identity",width = 800,position ="fill")+
  scale_fill_manual(values = col1)+
  labs(x = "",y = "All effective population size",size = rel(1)) +
  #geom_text(aes(label = datafd$Num),position=position_stack(vjust =0.5),size = 5)+ 这行代码是显示个数的数字
  guides(fill = guide_legend(reverse = F))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,size = rel(0.8))) +
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
