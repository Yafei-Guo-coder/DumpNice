#working directory
#203:yafei:/data1/home/yafei/008_Software/snpEff/
#input file: Xp-clr_9VIP/* from 204:yafei:/data2/yafei/003_Project3/Vmap1.1/E6/Xp-CLR_V2/VIP_genes
#output file: data/Xp-clr_9VIP_out/* 
#shell: 
for i in `cat /data1/home/yafei/008_Software/snpEff/Xp-clr_9VIP/names.txt`
do
  java -Xmx10G -jar ../snpEff.jar eff -c ../snpEff.config AT_10 /data1/home/yafei/008_Software/snpEff/Xp-clr_9VIP/${i}.pos.recode.vcf > ${i}.snp.eff.vcf -csvStats ${i}.csv -stats ${i}.html &
done
for i in `cat /data1/home/yafei/008_Software/snpEff/Xp-clr_9VIP/names.txt`
do
  sed '1,45d' ${i}.snp.eff.vcf | awk '{if(NF > 100) print $1"\t"$2"\t"$4"\t"$5}' > ${i}.snpEff1
  sed '1,45d' ${i}.snp.eff.vcf | awk '{if(NF>100) print $8}' |awk -F"ANN=" '{print $2}' | awk -F"|" '{print $2"\t"$3"\t"$11}' > ${i}.snpEff2
  paste -d "\t" ${i}.snpEff1 ${i}.snpEff2 > ${i}.snpEff
  rm ${i}.snpEff1 ${i}.snpEff2
done

#画变异序列分布
#将VCF通过tassel转换成table格式的table.txt
for i in `cat /data1/home/yafei/008_Software/snpEff/Xp-clr_9VIP/names.txt`
do
  sed '1d' table/${i}.txt | sed 's/\t/:/g' | sed 's/:/\t/1' | sed 's/://g' > ${i}.logo.seq
done

library(ggplot2)
library(ggseqlogo)
setwd("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment/02_XP-CLR/Gene/VIP_gene/Plot")
substrRight <- function(x){
  num = nchar(x)-2
  gsub('[0-9.]', '', substr(x, 6, nchar(x)))
}
names <- read.table("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment/02_XP-CLR/Gene/VIP_gene/names.txt",header=F,stringsAsFactors = F)
names <- names[,1]
#names <- c("","","","","","","","","","","","","","","","","GA2ox3-A1","","GA2ox3-D1","","","","","","","","","","","","","","","","","","","","","","","","","","","","NGR5_1","NGR5_2","NGR5_3","PIF_1","PIF_2","PIF_3")

pdf("plot1.pdf", width = 60, height = 8)
for (i in c(1:9)){
  file1 <- paste("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment/02_XP-CLR/Gene/VIP_gene/SeqLog/",names[i],".logo.seq",sep="")
  fasta_file = read.table(file1,header=F,stringsAsFactors = F)
  fasta = fasta_file[,2]
  file2 <- paste("/Users/guoyafei/Documents/01_个人项目/01_Migration/02_Add_ZNdata/02_Environment/02_XP-CLR/Gene/VIP_gene/SnpEff/",names[i],".snpEff",sep="")
  snpEff <- read.table(file2,header=F,stringsAsFactors = F,fill=TRUE,sep="\t")
  p1 <- ggseqlogo(fasta,method="prob")+theme(axis.text.x = element_blank())+labs(title = names[i])+theme(plot.title = element_text(hjust = 0.5,size = 25))
  Ref <- as.data.frame(snpEff$V3)
  colnames(Ref) <- "letter"
  Alt <- as.data.frame(snpEff$V4)
  colnames(Alt) <- "letter"
  if(dim(snpEff)[2] >6){
    Old <- as.data.frame(substring(snpEff$V7,3,5)) 
    colnames(Old) <- "letter"
    Pos <- as.data.frame(gsub('[a-zA-Z.*]', '', snpEff$V7))
    colnames(Pos) <- "letter"
    New <- as.data.frame(substrRight(snpEff$V7))
    colnames(New) <- "letter"
    all <- rbind(Ref,Alt,Old,Pos,New)
    num <- dim(Ref)[1]
    aln <- data.frame(
      letter = all,
      Type=rep(c("Ref","Alt","Old","Pos","New"), each=num),
      x=rep(1:num,5)
    )
    p2 <- ggplot(aln, aes(x, Type)) +
      geom_text(aes(label=letter,)) +
      scale_y_discrete(limits=c("New","Pos","Old","Alt","Ref"))+
      scale_x_continuous(breaks=1:10, expand = c(0.105, 0)) + xlab('') +
      theme_logo() +
      theme(legend.position = 'none', axis.text.x = element_blank())
  } else{
    all <- rbind(Ref,Alt)
    num <- dim(Ref)[1]
    aln <- data.frame(
      letter = all,
      species=rep(c("Ref","Alt"), each=num),
      x=rep(1:num,2)
    )
    p2 <- ggplot(aln, aes(x, Type)) +
      geom_text(aes(label=letter,),check_overlap = TRUE) +
      scale_y_discrete(limits=c("Alt","Ref"))+
      scale_x_continuous(breaks=1:10, expand = c(0.105, 0)) + xlab('') +
      theme_logo() +
      theme(legend.position = 'none', axis.text.x = element_blank())  
  }
  snpEff$impact <- NA
  snpEff[which(snpEff$V6=="HIGH"),8] <- 4
  snpEff[which(snpEff$V6=="MODERATE"),8] <- 3
  snpEff[which(snpEff$V6=="LOW"),8] <- 2
  snpEff[which(snpEff$V6=="MODIFILER"),8] <- 1
  bp_data <- data.frame(
    x=snpEff$V2,
    Impact=snpEff$impact
  )
  p3 <- ggplot(bp_data, aes(factor(x), Impact))+
    geom_bar(stat = "identity", fill="grey")+
    theme_logo()+
    xlab("")+
    theme(axis.text.x = element_text(angle = 90))
  suppressMessages(require(cowplot))
  p <- plot_grid(p1,p2,p3,ncol = 1, align = "v")
  print(p)
  
}
dev.off()

