library(ggplot2)
library(ggseqlogo)
setwd("/Users/guoyafei/Documents/个人项目/傅老师/20210311/GWAS_sign_gene/")
substrRight <- function(x){
  num = nchar(x)-2
  gsub('[0-9.]', '', substr(x, 6, nchar(x)))
}
names <- c("","","","","","","","","","","","","","","","","GA2ox3-A1","","GA2ox3-D1","","","","","","","","","","","","","","","","","","","","","","","","","","","","NGR5_1","NGR5_2","NGR5_3","PIF_1","PIF_2","PIF_3")
pdf("plot1.pdf", width = 60, height = 8)
for (i in c(17,19,47,48,49,50,51,52)){
  file1 <- paste("/Users/guoyafei/Documents/个人项目/傅老师/20210311/GWAS_sign_gene/TXT/noheader/",i,".txt.log",sep="")
  fasta_file = read.table(file1,header=F,stringsAsFactors = F)
  fasta = fasta_file[,2]
  file2 <- paste("/Users/guoyafei/Documents/个人项目/傅老师/20210311/GWAS_sign_gene/snpEff/",i,".pos.snpEff",sep="")
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

