#将hapmap转换成vcf的过程中，使用
#run_pipeline.pl -SortGenotypeFilePlugin -inputFile test.hmp.txt -outputFile test.sort.hmp.txt -fileType Hapmap
#排序时，部分等位基因的前后位置发生了转换
#由下列代码将World.hmp.txt转换为vcf文件，output_World.vcf。

setwd("/Users/guoyafei/Downloads/Results/ChangeFormat/")
T <- read.table("alleles", header =FALSE, stringsAsFactors = FALSE,sep=" ")
data <- read.table("output.txt", header =TRUE, stringsAsFactors = F)
colnames(T) <- c("psId","alleA","alleB")
World <- merge(T, data, by="psId", all.x=TRUE)
AA <- paste(World$alleA,World$alleA,sep="")
AB <- paste(World$alleA,World$alleB,sep="")
BB <- paste(World$alleB,World$alleB,sep="")
World$AA <- AA
World$AB <- AB
World$BB <- BB

write.table(World,"vcf_input.txt", quote = FALSE, row.names = FALSE, sep = " ")
#使用CreateVCF.java来创建vcf_output.txt文件
input <- read.table("vcf_output.txt", header =TRUE, stringsAsFactors = FALSE, sep=" ")
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
#out <- read.table("vcf_input.txt", header =TRUE, stringsAsFactors = FALSE, sep=" ")
names <- colnames(input)
colnames(input) <- c(names[-1],"QUAL")
input[,756] <- "."
input$CHROM <- input[,7]
input$POS <- input[,8]
input$ID <- rownames(input)
input$REF <- input[,1]
input$ALT <- input[,2]
input$FILTER <- "PASS"
input$INFO <- "."
input$FORMAT <- "GT"
new <- input[,c(757:761,756,762:764,10:755)]
new2 <- new[order(new[,2]),]
new3 <- new2[order(new2[,1]),]
write.table(new3,"output_World.vcf", quote = FALSE, row.names = FALSE, sep = "\t")


