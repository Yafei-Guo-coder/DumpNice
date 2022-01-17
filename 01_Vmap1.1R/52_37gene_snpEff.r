#取这些基因上下游5k的snp
#位点
#yafei@204: /data1/home/yafei/003_Project3/bayenv/Merge/Gene
#snp
#yafei@204: /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225
#chr1.E6_Landrace_locate.vcf.gz
cat $1 |while read chr from to gene
do
vcftools --gzvcf /data2/yafei/003_Project3/Vmap1.1/E6/Landrace_locate_225/chr${chr}.E6_Landrace_locate.vcf.gz --chr $chr --from-bp $from --to-bp $to  --recode --recode-INFO-all --out 5k.$gene.$chr.$from-$to
done
#将VCF转移到 yafei@203: /data1/home/yafei/008_Software/snpEff/37Gene
for i in `cat /data1/home/yafei/008_Software/snpEff/37Gene/vcfName.txt`
do
grep -v "#" ${i}.snp.eff.vcf | awk '{if(NF > 100) print $1"\t"$2"\t"$4"\t"$5}' > ${i}.snpEff1
grep -v "#" ${i}.snp.eff.vcf | awk '{if(NF > 100) print $8}' |awk -F"ANN=" '{print $2}' | awk -F"|" '{print $2"\t"$3"\t"$11}' > ${i}.snpEff2
paste -d "\t" ${i}.snpEff1 ${i}.snpEff2 > ${i}.snpEff
rm ${i}.snpEff1 ${i}.snpEff2
done
#将VCF文件和snpEff文件转移到本地: /Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff

#画变异序列分布
for i in `cat /data1/home/yafei/008_Software/snpEff/37Gene/vcfName.txt`
do
run_pipeline.pl -Xms512m -Xmx5g -vcf ${i}.snp.eff.vcf  -export ${i}.tassel -exportType Table
sed '1d' ${i}.tassel.txt | sed 's/\t/:/g' | sed 's/:/\t/1' | sed 's/://g' > ${i}.logo.seq
done

library(ggplot2)
library(ggseqlogo)
library(cowplot)
setwd("/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff")
substrRight <- function(x){
  num = nchar(x)-2
  gsub('[0-9.]', '', substr(x, 6, nchar(x)))
}
names <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/name.txt",header=F,stringsAsFactors = F,sep="\t")
pdf_tit <- read.table("/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/pdf_title.txt", header=F,stringsAsFactors = F)
names <- names[,1]
#names <- c("","","","","","","","","","","","","","","","","GA2ox3-A1","","GA2ox3-D1","","","","","","","","","","","","","","","","","","","","","","","","","","","","NGR5_1","NGR5_2","NGR5_3","PIF_1","PIF_2","PIF_3")
pdf("Glu-1A.pdf", width = 13, height = 8)
#for (i in c(1:length(names))){
for (i in c(2)){
    file1 <- paste("/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/Logo_seq/",names[i],".logo.seq",sep="")
    fasta_file = read.table(file1,header=F,stringsAsFactors = F)
    fasta = fasta_file[,2]
    file2 <- paste("/Users/guoyafei/Documents/01_Migration/02_Environment/08_snpEff/snpEff/",names[i],".snpEff",sep="")
    snpEff <- read.table(file2,header=F,stringsAsFactors = F,fill=TRUE,sep="\t")
    p1 <- ggseqlogo(fasta,method="prob")+
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position="none")+
      labs(title = pdf_tit[i,1])+
      theme(plot.title = element_text(hjust = 0.5,size = 15))
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
        scale_x_continuous(breaks=1:10, expand = c(0.07, 0)) + xlab('') +
        theme_logo() +
        theme(legend.position = 'none', axis.text.x = element_blank(),)
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
        scale_x_continuous(breaks=1:10, expand = c(0.07, 0)) + xlab('') +
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
    #bp_data$x <- c(1:dim(bp_data)[1])
    p3 <- ggplot(bp_data, aes(x, Impact))+
      geom_bar(stat = "identity", fill="grey")+
      theme_logo()+
      #scale_x_discrete( expand = c(0.055, 0))+
      xlab("")+
      theme(axis.text.x = element_text(angle = 90,size=15))
    suppressMessages(require(cowplot))
    p <- plot_grid(p1,p2,p3,ncol = 1, align = "v")
    print(p)
}
dev.off()
