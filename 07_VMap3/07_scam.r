require(scam)
testx=data1$physicalPosition
length(testx)
testy=data1$geneticPosition
length(testy)
fit = scam(testy~s(testx,k=30,bs="cr"), family=gaussian(link="identity"))
predict(fit,data.frame(testx))
#plot(testx,testy)
#lines(testx,predict(fit),col="red")

plot(testx,predict(fit,data.frame(testx)))
#head(data.frame(testx))
#z<-diff(predict(fit,data.frame(testx)))
#sort(z)
#zz<-data.frame(testx,testy,predict(fit,data.frame(testx)))
#summary(zz)
#zz%>%
#  gather(key="group",value="y",c(2,3)) %>%
#  ggplot(aes(x=testx,y=y,color=group))+geom_point(size=0.1)
#zz%>%ggplot(aes(x=predict.fit..data.frame.testx..,y=testy))+geom_point()
#test<-c(1157133,2410604)
#approx(df1$physicalPosition,df1$geneticPosition,xout = test,method = "linear",rule = 1,ties = "ordered")$y
#head(df1$physicalPosition)
#zlist<- unique(df$chromosome)
#for(p in zlist){
#  print(p)
#}
#head(df)
require(scam)
setwd("/Users/guoyafei/Documents/02_VmapIII/06_Selection/selscan")
map <- read.table("/Users/guoyafei/公共资源/01_数据资源/iwgsc_refseqv1.0_mapping_data_42chr.txt",header=T,stringsAsFactors = F)
filenames <- paste(c(paste("chr00",1:9,sep=""),paste("chr0",10:42,sep="")),".pos.txt",sep="")
for(i in c(6)){
  sub <- map[which(map$chromosome == i),]
  pos <- read.table(filenames[i],header=F,stringsAsFactors = F)
  colnames(pos) <- c("chr","physicalPosition")
  physicalPosition=sub$physicalPosition
  geneticPosition=sub$geneticPosition
  fit = scam(geneticPosition~s(physicalPosition,k=20,bs="mpi"), family=gaussian(link="identity"))
  #plot(physicalPosition,predict(fit,data.frame(geneticPosition)))
  ID <- paste(pos$chr,pos$physicalPosition,sep="-")
  physicalPosition <- pos$physicalPosition
  out <- data.frame(pos$chr,ID,predict(fit,data.frame(physicalPosition)),pos$physicalPosition)
  outfile <- paste("/Users/guoyafei/Documents/02_VmapIII/06_Selection/selscan/map/",filenames[i],sep="")
  write.table(out,outfile,col.names = F,quote = F,sep="\t", row.names = F)
}

#chr036.ihs.out
data <- read.table("/Users/guoyafei/Desktop/chr036.ihs.out",header=F,stringsAsFactors = F)
colnames(data)<-c("locus","pos","1_freq","ihh_1","ihh_0","ihs")
data2 <- data[sort(sample(c(1:52691),10000)),]
p<- ggplot(data2, aes(x=pos/1000000, y=ihs)) +
      # Show all points
      geom_point(alpha=0.8, size=1.3) +
      #scale_color_manual(values = rep(c("grey","grey", "skyblue", "skyblue"), 7)) +
      # custom X axis:
      #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
      # Add highlighted points
      #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
      # Add label using ggrepel to avoid overlapping
      #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size = 15),
        axis.title.x=element_text(size = 15),
      )
      #scale_y_continuous(limits = c(0,7))+
      #geom_point(data=point,aes(x=BPcum,y=-log10(P)),color="red")
    #geom_vline(xintercept = 528377322, colour="red",linetype=2, size=1)
    #geom_hline(yintercept = -log10(thresh[a,1]), colour="red",linetype=2, size=1)+
    #geom_hline(yintercept = -log10(thresh[a,2]), colour="blue",linetype=2, size=1)














