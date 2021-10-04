#1. 提取需要的样本
EC <- read.table("EC.vcf", header=T,stringsAsFactors=F)
WW <- read.table("WW.vcf", header=T,stringsAsFactors=F)
TR <- read.table("TR.vcf", header = T,stringsAsFactors=F)
names <- read.table("855names", header=F, stringsAsFactors=F)
#dim(TR[1,colnames(TR) %in% names[,1]])
#TR[1,colnames(TR) %in% names[,1]]
#WW[1,colnames(WW) %in% names[,1]]
#dim(WW[1,colnames(WW) %in% names[,1]])

new_WW <- WW[,colnames(WW) %in% names[,1]]
new_TR <- TR[,colnames(TR) %in% names[,1]]
new_EC <- EC[,colnames(EC) %in% names[,1]]
#new_WW[1:3,1:10]
#new_TR[1:3,1:10]
#new_EC[1:3,1:10]

#2. 修改样本名字
WWn <- read.table("WW-new-name", header=F, stringsAsFactors=F)
#WWn[1,1] <- "#CHROM"
colnames(new_WW) <- WWn[,1]
#new_WW[1:3,1:10]
TRn <- read.table("TR-new-name", header=F, stringsAsFactors=F)
#TRn[1,1] <- "#CHROM"
colnames(new_TR) <- TRn[,1]
#new_TR[1:3,1:10]
ECn <- read.table("EC-new-name", header=F, stringsAsFactors=F)
#ECn[1,1] <- "#CHROM"
colnames(new_EC) <- ECn[,1]
#new_EC[1:3,1:10]

#3. 修改染色体号（21 -> 42）& POS 
new_WW[which(new_WW$CHROM == "chr1A"),1] <- 1
new_WW[which(new_WW$CHROM == "chr1B"),1] <- 3
new_WW[which(new_WW$CHROM == "chr1D"),1] <- 5
new_WW[which(new_WW$CHROM == "chr2A"),1] <- 7
new_WW[which(new_WW$CHROM == "chr2B"),1] <- 9
new_WW[which(new_WW$CHROM == "chr2D"),1] <- 11
new_WW[which(new_WW$CHROM == "chr3A"),1] <- 13
new_WW[which(new_WW$CHROM == "chr3B"),1] <- 15
new_WW[which(new_WW$CHROM == "chr3D"),1] <- 17
new_WW[which(new_WW$CHROM == "chr4A"),1] <- 19
new_WW[which(new_WW$CHROM == "chr4B"),1] <- 21
new_WW[which(new_WW$CHROM == "chr4D"),1] <- 23
new_WW[which(new_WW$CHROM == "chr5A"),1] <- 25
new_WW[which(new_WW$CHROM == "chr5B"),1] <- 27
new_WW[which(new_WW$CHROM == "chr5D"),1] <- 29
new_WW[which(new_WW$CHROM == "chr6A"),1] <- 31
new_WW[which(new_WW$CHROM == "chr6B"),1] <- 33
new_WW[which(new_WW$CHROM == "chr6D"),1] <- 35
new_WW[which(new_WW$CHROM == "chr7A"),1] <- 37
new_WW[which(new_WW$CHROM == "chr7B"),1] <- 39
new_WW[which(new_WW$CHROM == "chr7D"),1] <- 41

new_WW[which(new_WW$CHROM == "1" & new_WW$POS > 471304005),1] <- 2
pos <- new_WW[which(new_WW$CHROM == "2"),2]
new_WW[which(new_WW$CHROM == "2"),2] <- pos - 471304005

new_WW[which(new_WW$CHROM == "3" & new_WW$POS > 438720154),1] <- 4
pos <- new_WW[which(new_WW$CHROM == "4"),2]
new_WW[which(new_WW$CHROM == "4"),2] <- pos - 438720154

new_WW[which(new_WW$CHROM == "5" & new_WW$POS > 452179604),1] <- 6
pos <- new_WW[which(new_WW$CHROM == "6"),2]
new_WW[which(new_WW$CHROM == "6"),2] <- pos - 452179604

new_WW[which(new_WW$CHROM == "7" & new_WW$POS > 462376173),1] <- 8
pos <- new_WW[which(new_WW$CHROM == "8"),2]
new_WW[which(new_WW$CHROM == "8"),2] <- pos - 462376173

new_WW[which(new_WW$CHROM == "9" & new_WW$POS > 453218924),1] <- 10
pos <- new_WW[which(new_WW$CHROM == "10"),2]
new_WW[which(new_WW$CHROM == "10"),2] <- pos - 453218924

new_WW[which(new_WW$CHROM == "11" & new_WW$POS > 462216879),1] <- 12
pos <- new_WW[which(new_WW$CHROM == "12"),2]
new_WW[which(new_WW$CHROM == "12"),2] <- pos - 462216879

new_WW[which(new_WW$CHROM == "13" & new_WW$POS > 454103970),1] <- 14
pos <- new_WW[which(new_WW$CHROM == "14"),2]
new_WW[which(new_WW$CHROM == "14"),2] <- pos - 454103970

new_WW[which(new_WW$CHROM == "15" & new_WW$POS > 448155269),1] <- 16
pos <- new_WW[which(new_WW$CHROM == "16"),2]
new_WW[which(new_WW$CHROM == "16"),2] <- pos - 448155269

new_WW[which(new_WW$CHROM == "17" & new_WW$POS > 476235359),1] <- 18
pos <- new_WW[which(new_WW$CHROM == "18"),2]
new_WW[which(new_WW$CHROM == "18"),2] <- pos - 476235359

new_WW[which(new_WW$CHROM == "19" & new_WW$POS > 452555092),1] <- 20
pos <- new_WW[which(new_WW$CHROM == "20"),2]
new_WW[which(new_WW$CHROM == "20"),2] <- pos - 452555092

new_WW[which(new_WW$CHROM == "21" & new_WW$POS > 451014251),1] <- 22
pos <- new_WW[which(new_WW$CHROM == "22"),2]
new_WW[which(new_WW$CHROM == "22"),2] <- pos - 451014251

new_WW[which(new_WW$CHROM == "23" & new_WW$POS > 451004620),1] <- 24
pos <- new_WW[which(new_WW$CHROM == "24"),2]
new_WW[which(new_WW$CHROM == "24"),2] <- pos - 451004620

new_WW[which(new_WW$CHROM == "25" & new_WW$POS > 453230519),1] <- 26
pos <- new_WW[which(new_WW$CHROM == "26"),2]
new_WW[which(new_WW$CHROM == "26"),2] <- pos - 453230519

new_WW[which(new_WW$CHROM == "27" & new_WW$POS > 451372872),1] <- 28
pos <- new_WW[which(new_WW$CHROM == "28"),2]
new_WW[which(new_WW$CHROM == "28"),2] <- pos - 451372872

new_WW[which(new_WW$CHROM == "29" & new_WW$POS > 451901030),1] <- 30
pos <- new_WW[which(new_WW$CHROM == "30"),2]
new_WW[which(new_WW$CHROM == "30"),2] <- pos - 451901030

new_WW[which(new_WW$CHROM == "31" & new_WW$POS > 452440856),1] <- 32
pos <- new_WW[which(new_WW$CHROM == "32"),2]
new_WW[which(new_WW$CHROM == "32"),2] <- pos - 452440856

new_WW[which(new_WW$CHROM == "33" & new_WW$POS > 452077197),1] <- 34
pos <- new_WW[which(new_WW$CHROM == "34"),2]
new_WW[which(new_WW$CHROM == "34"),2] <- pos - 452077197

new_WW[which(new_WW$CHROM == "35" & new_WW$POS > 450509124),1] <- 36
pos <- new_WW[which(new_WW$CHROM == "36"),2]
new_WW[which(new_WW$CHROM == "36"),2] <- pos - 450509124

new_WW[which(new_WW$CHROM == "37" & new_WW$POS > 450046986),1] <- 38
pos <- new_WW[which(new_WW$CHROM == "38"),2]
new_WW[which(new_WW$CHROM == "38"),2] <- pos - 450046986

new_WW[which(new_WW$CHROM == "39" & new_WW$POS > 453822637),1] <- 40
pos <- new_WW[which(new_WW$CHROM == "40"),2]
new_WW[which(new_WW$CHROM == "40"),2] <- pos - 453822637

new_WW[which(new_WW$CHROM == "41" & new_WW$POS > 453812268),1] <- 42
pos <- new_WW[which(new_WW$CHROM == "42"),2]
new_WW[which(new_WW$CHROM == "42"),2] <- pos - 453812268

new_EC[which(new_EC$CHROM == "chr1A"),1] <- 1
new_EC[which(new_EC$CHROM == "chr1B"),1] <- 3
new_EC[which(new_EC$CHROM == "chr1D"),1] <- 5
new_EC[which(new_EC$CHROM == "chr2A"),1] <- 7
new_EC[which(new_EC$CHROM == "chr2B"),1] <- 9
new_EC[which(new_EC$CHROM == "chr2D"),1] <- 11
new_EC[which(new_EC$CHROM == "chr3A"),1] <- 13
new_EC[which(new_EC$CHROM == "chr3B"),1] <- 15
new_EC[which(new_EC$CHROM == "chr3D"),1] <- 17
new_EC[which(new_EC$CHROM == "chr4A"),1] <- 19
new_EC[which(new_EC$CHROM == "chr4B"),1] <- 21
new_EC[which(new_EC$CHROM == "chr4D"),1] <- 23
new_EC[which(new_EC$CHROM == "chr5A"),1] <- 25
new_EC[which(new_EC$CHROM == "chr5B"),1] <- 27
new_EC[which(new_EC$CHROM == "chr5D"),1] <- 29
new_EC[which(new_EC$CHROM == "chr6A"),1] <- 31
new_EC[which(new_EC$CHROM == "chr6B"),1] <- 33
new_EC[which(new_EC$CHROM == "chr6D"),1] <- 35
new_EC[which(new_EC$CHROM == "chr7A"),1] <- 37
new_EC[which(new_EC$CHROM == "chr7B"),1] <- 39
new_EC[which(new_EC$CHROM == "chr7D"),1] <- 41

new_EC[which(new_EC$CHROM == "1" & new_EC$POS > 471304005),1] <- 2
pos <- new_EC[which(new_EC$CHROM == "2"),2]
new_EC[which(new_EC$CHROM == "2"),2] <- pos - 471304005

new_EC[which(new_EC$CHROM == "3" & new_EC$POS > 438720154),1] <- 4
pos <- new_EC[which(new_EC$CHROM == "4"),2]
new_EC[which(new_EC$CHROM == "4"),2] <- pos - 438720154

new_EC[which(new_EC$CHROM == "5" & new_EC$POS > 452179604),1] <- 6
pos <- new_EC[which(new_EC$CHROM == "6"),2]
new_EC[which(new_EC$CHROM == "6"),2] <- pos - 452179604

new_EC[which(new_EC$CHROM == "7" & new_EC$POS > 462376173),1] <- 8
pos <- new_EC[which(new_EC$CHROM == "8"),2]
new_EC[which(new_EC$CHROM == "8"),2] <- pos - 462376173

new_EC[which(new_EC$CHROM == "9" & new_EC$POS > 453218924),1] <- 10
pos <- new_EC[which(new_EC$CHROM == "10"),2]
new_EC[which(new_EC$CHROM == "10"),2] <- pos - 453218924

new_EC[which(new_EC$CHROM == "11" & new_EC$POS > 462216879),1] <- 12
pos <- new_EC[which(new_EC$CHROM == "12"),2]
new_EC[which(new_EC$CHROM == "12"),2] <- pos - 462216879

new_EC[which(new_EC$CHROM == "13" & new_EC$POS > 454103970),1] <- 14
pos <- new_EC[which(new_EC$CHROM == "14"),2]
new_EC[which(new_EC$CHROM == "14"),2] <- pos - 454103970

new_EC[which(new_EC$CHROM == "15" & new_EC$POS > 448155269),1] <- 16
pos <- new_EC[which(new_EC$CHROM == "16"),2]
new_EC[which(new_EC$CHROM == "16"),2] <- pos - 448155269

new_EC[which(new_EC$CHROM == "17" & new_EC$POS > 476235359),1] <- 18
pos <- new_EC[which(new_EC$CHROM == "18"),2]
new_EC[which(new_EC$CHROM == "18"),2] <- pos - 476235359

new_EC[which(new_EC$CHROM == "19" & new_EC$POS > 452555092),1] <- 20
pos <- new_EC[which(new_EC$CHROM == "20"),2]
new_EC[which(new_EC$CHROM == "20"),2] <- pos - 452555092

new_EC[which(new_EC$CHROM == "21" & new_EC$POS > 451014251),1] <- 22
pos <- new_EC[which(new_EC$CHROM == "22"),2]
new_EC[which(new_EC$CHROM == "22"),2] <- pos - 451014251

new_EC[which(new_EC$CHROM == "23" & new_EC$POS > 451004620),1] <- 24
pos <- new_EC[which(new_EC$CHROM == "24"),2]
new_EC[which(new_EC$CHROM == "24"),2] <- pos - 451004620

new_EC[which(new_EC$CHROM == "25" & new_EC$POS > 453230519),1] <- 26
pos <- new_EC[which(new_EC$CHROM == "26"),2]
new_EC[which(new_EC$CHROM == "26"),2] <- pos - 453230519

new_EC[which(new_EC$CHROM == "27" & new_EC$POS > 451372872),1] <- 28
pos <- new_EC[which(new_EC$CHROM == "28"),2]
new_EC[which(new_EC$CHROM == "28"),2] <- pos - 451372872

new_EC[which(new_EC$CHROM == "29" & new_EC$POS > 451901030),1] <- 30
pos <- new_EC[which(new_EC$CHROM == "30"),2]
new_EC[which(new_EC$CHROM == "30"),2] <- pos - 451901030

new_EC[which(new_EC$CHROM == "31" & new_EC$POS > 452440856),1] <- 32
pos <- new_EC[which(new_EC$CHROM == "32"),2]
new_EC[which(new_EC$CHROM == "32"),2] <- pos - 452440856

new_EC[which(new_EC$CHROM == "33" & new_EC$POS > 452077197),1] <- 34
pos <- new_EC[which(new_EC$CHROM == "34"),2]
new_EC[which(new_EC$CHROM == "34"),2] <- pos - 452077197

new_EC[which(new_EC$CHROM == "35" & new_EC$POS > 450509124),1] <- 36
pos <- new_EC[which(new_EC$CHROM == "36"),2]
new_EC[which(new_EC$CHROM == "36"),2] <- pos - 450509124

new_EC[which(new_EC$CHROM == "37" & new_EC$POS > 450046986),1] <- 38
pos <- new_EC[which(new_EC$CHROM == "38"),2]
new_EC[which(new_EC$CHROM == "38"),2] <- pos - 450046986

new_EC[which(new_EC$CHROM == "39" & new_EC$POS > 453822637),1] <- 40
pos <- new_EC[which(new_EC$CHROM == "40"),2]
new_EC[which(new_EC$CHROM == "40"),2] <- pos - 453822637

new_EC[which(new_EC$CHROM == "41" & new_EC$POS > 453812268),1] <- 42
pos <- new_EC[which(new_EC$CHROM == "42"),2]
new_EC[which(new_EC$CHROM == "42"),2] <- pos - 453812268

new_TR[which(new_TR$CHROM == "chr1A"),1] <- 1
new_TR[which(new_TR$CHROM == "chr1B"),1] <- 3
new_TR[which(new_TR$CHROM == "chr1D"),1] <- 5
new_TR[which(new_TR$CHROM == "chr2A"),1] <- 7
new_TR[which(new_TR$CHROM == "chr2B"),1] <- 9
new_TR[which(new_TR$CHROM == "chr2D"),1] <- 11
new_TR[which(new_TR$CHROM == "chr3A"),1] <- 13
new_TR[which(new_TR$CHROM == "chr3B"),1] <- 15
new_TR[which(new_TR$CHROM == "chr3D"),1] <- 17
new_TR[which(new_TR$CHROM == "chr4A"),1] <- 19
new_TR[which(new_TR$CHROM == "chr4B"),1] <- 21
new_TR[which(new_TR$CHROM == "chr4D"),1] <- 23
new_TR[which(new_TR$CHROM == "chr5A"),1] <- 25
new_TR[which(new_TR$CHROM == "chr5B"),1] <- 27
new_TR[which(new_TR$CHROM == "chr5D"),1] <- 29
new_TR[which(new_TR$CHROM == "chr6A"),1] <- 31
new_TR[which(new_TR$CHROM == "chr6B"),1] <- 33
new_TR[which(new_TR$CHROM == "chr6D"),1] <- 35
new_TR[which(new_TR$CHROM == "chr7A"),1] <- 37
new_TR[which(new_TR$CHROM == "chr7B"),1] <- 39
new_TR[which(new_TR$CHROM == "chr7D"),1] <- 41

new_TR[which(new_TR$CHROM == "1" & new_TR$POS > 471304005),1] <- 2
pos <- new_TR[which(new_TR$CHROM == "2"),2]
new_TR[which(new_TR$CHROM == "2"),2] <- pos - 471304005

new_TR[which(new_TR$CHROM == "3" & new_TR$POS > 438720154),1] <- 4
pos <- new_TR[which(new_TR$CHROM == "4"),2]
new_TR[which(new_TR$CHROM == "4"),2] <- pos - 438720154

new_TR[which(new_TR$CHROM == "5" & new_TR$POS > 452179604),1] <- 6
pos <- new_TR[which(new_TR$CHROM == "6"),2]
new_TR[which(new_TR$CHROM == "6"),2] <- pos - 452179604

new_TR[which(new_TR$CHROM == "7" & new_TR$POS > 462376173),1] <- 8
pos <- new_TR[which(new_TR$CHROM == "8"),2]
new_TR[which(new_TR$CHROM == "8"),2] <- pos - 462376173

new_TR[which(new_TR$CHROM == "9" & new_TR$POS > 453218924),1] <- 10
pos <- new_TR[which(new_TR$CHROM == "10"),2]
new_TR[which(new_TR$CHROM == "10"),2] <- pos - 453218924

new_TR[which(new_TR$CHROM == "11" & new_TR$POS > 462216879),1] <- 12
pos <- new_TR[which(new_TR$CHROM == "12"),2]
new_TR[which(new_TR$CHROM == "12"),2] <- pos - 462216879

new_TR[which(new_TR$CHROM == "13" & new_TR$POS > 454103970),1] <- 14
pos <- new_TR[which(new_TR$CHROM == "14"),2]
new_TR[which(new_TR$CHROM == "14"),2] <- pos - 454103970

new_TR[which(new_TR$CHROM == "15" & new_TR$POS > 448155269),1] <- 16
pos <- new_TR[which(new_TR$CHROM == "16"),2]
new_TR[which(new_TR$CHROM == "16"),2] <- pos - 448155269

new_TR[which(new_TR$CHROM == "17" & new_TR$POS > 476235359),1] <- 18
pos <- new_TR[which(new_TR$CHROM == "18"),2]
new_TR[which(new_TR$CHROM == "18"),2] <- pos - 476235359

new_TR[which(new_TR$CHROM == "19" & new_TR$POS > 452555092),1] <- 20
pos <- new_TR[which(new_TR$CHROM == "20"),2]
new_TR[which(new_TR$CHROM == "20"),2] <- pos - 452555092

new_TR[which(new_TR$CHROM == "21" & new_TR$POS > 451014251),1] <- 22
pos <- new_TR[which(new_TR$CHROM == "22"),2]
new_TR[which(new_TR$CHROM == "22"),2] <- pos - 451014251

new_TR[which(new_TR$CHROM == "23" & new_TR$POS > 451004620),1] <- 24
pos <- new_TR[which(new_TR$CHROM == "24"),2]
new_TR[which(new_TR$CHROM == "24"),2] <- pos - 451004620

new_TR[which(new_TR$CHROM == "25" & new_TR$POS > 453230519),1] <- 26
pos <- new_TR[which(new_TR$CHROM == "26"),2]
new_TR[which(new_TR$CHROM == "26"),2] <- pos - 453230519

new_TR[which(new_TR$CHROM == "27" & new_TR$POS > 451372872),1] <- 28
pos <- new_TR[which(new_TR$CHROM == "28"),2]
new_TR[which(new_TR$CHROM == "28"),2] <- pos - 451372872

new_TR[which(new_TR$CHROM == "29" & new_TR$POS > 451901030),1] <- 30
pos <- new_TR[which(new_TR$CHROM == "30"),2]
new_TR[which(new_TR$CHROM == "30"),2] <- pos - 451901030

new_TR[which(new_TR$CHROM == "31" & new_TR$POS > 452440856),1] <- 32
pos <- new_TR[which(new_TR$CHROM == "32"),2]
new_TR[which(new_TR$CHROM == "32"),2] <- pos - 452440856

new_TR[which(new_TR$CHROM == "33" & new_TR$POS > 452077197),1] <- 34
pos <- new_TR[which(new_TR$CHROM == "34"),2]
new_TR[which(new_TR$CHROM == "34"),2] <- pos - 452077197

new_TR[which(new_TR$CHROM == "35" & new_TR$POS > 450509124),1] <- 36
pos <- new_TR[which(new_TR$CHROM == "36"),2]
new_TR[which(new_TR$CHROM == "36"),2] <- pos - 450509124

new_TR[which(new_TR$CHROM == "37" & new_TR$POS > 450046986),1] <- 38
pos <- new_TR[which(new_TR$CHROM == "38"),2]
new_TR[which(new_TR$CHROM == "38"),2] <- pos - 450046986

new_TR[which(new_TR$CHROM == "39" & new_TR$POS > 453822637),1] <- 40
pos <- new_TR[which(new_TR$CHROM == "40"),2]
new_TR[which(new_TR$CHROM == "40"),2] <- pos - 453822637

new_TR[which(new_TR$CHROM == "41" & new_TR$POS > 453812268),1] <- 42
pos <- new_TR[which(new_TR$CHROM == "42"),2]
new_TR[which(new_TR$CHROM == "42"),2] <- pos - 453812268

#4. 修改 ID & 去掉 chrUn
new_TR$ID <- paste(new_TR$CHROM, new_TR$POS, sep = "-")
new_WW$ID <- paste(new_WW$CHROM, new_WW$POS, sep = "-")
new_EC$ID <- paste(new_EC$CHROM, new_EC$POS, sep = "-")

new_TR <- new_TR[-which(new_TR$CHROM == "chrUn"), ]
new_WW <- new_WW[-which(new_WW$CHROM == "chrUn"), ]
new_EC <- new_EC[-which(new_EC$CHROM == "chrUn"), ]

new_TR$INFO <- "."
new_WW$INFO <- "."
new_EC$INFO <- "."

new_TR$FILTER <- "PASS"
new_WW$FILTER <- "PASS"
new_EC$FILTER <- "PASS"
write.table(new_WW,"new_WW", row.names=F,sep="\t",quote=F) 
write.table(new_TR,"new_TR", row.names=F,sep="\t",quote=F) 
write.table(new_EC,"new_EC", row.names=F,sep="\t",quote=F) 

#5. 拆分染色体
new_WW <- read.table("new_WW",header=T,stringsAsFactors=F)
new_TR <- read.table("new_TR",header=T,stringsAsFactors=F)
new_EC <- read.table("new_EC",header=T,stringsAsFactors=F)

WW_1 <- new_WW[which(new_WW$CHROM == "1"),]
WW_2 <- new_WW[which(new_WW$CHROM == "2"),]
WW_3 <- new_WW[which(new_WW$CHROM == "3"),]
WW_4 <- new_WW[which(new_WW$CHROM == "4"),]
WW_5 <- new_WW[which(new_WW$CHROM == "5"),]
WW_6 <- new_WW[which(new_WW$CHROM == "6"),]
WW_7 <- new_WW[which(new_WW$CHROM == "7"),]
WW_8 <- new_WW[which(new_WW$CHROM == "8"),]
WW_9 <- new_WW[which(new_WW$CHROM == "9"),]
WW_10 <- new_WW[which(new_WW$CHROM == "10"),]
WW_11 <- new_WW[which(new_WW$CHROM == "11"),]
WW_12 <- new_WW[which(new_WW$CHROM == "12"),]
WW_13 <- new_WW[which(new_WW$CHROM == "13"),]
WW_14 <- new_WW[which(new_WW$CHROM == "14"),]
WW_15 <- new_WW[which(new_WW$CHROM == "15"),]
WW_16 <- new_WW[which(new_WW$CHROM == "16"),]
WW_17 <- new_WW[which(new_WW$CHROM == "17"),]
WW_18 <- new_WW[which(new_WW$CHROM == "18"),]
WW_19 <- new_WW[which(new_WW$CHROM == "19"),]
WW_20 <- new_WW[which(new_WW$CHROM == "20"),]
WW_21 <- new_WW[which(new_WW$CHROM == "21"),]
WW_22 <- new_WW[which(new_WW$CHROM == "22"),]
WW_23 <- new_WW[which(new_WW$CHROM == "23"),]
WW_24 <- new_WW[which(new_WW$CHROM == "24"),]
WW_25 <- new_WW[which(new_WW$CHROM == "25"),]
WW_26 <- new_WW[which(new_WW$CHROM == "26"),]
WW_27 <- new_WW[which(new_WW$CHROM == "27"),]
WW_28 <- new_WW[which(new_WW$CHROM == "28"),]
WW_29 <- new_WW[which(new_WW$CHROM == "29"),]
WW_30 <- new_WW[which(new_WW$CHROM == "30"),]
WW_31 <- new_WW[which(new_WW$CHROM == "31"),]
WW_32 <- new_WW[which(new_WW$CHROM == "32"),]
WW_33 <- new_WW[which(new_WW$CHROM == "33"),]
WW_34 <- new_WW[which(new_WW$CHROM == "34"),]
WW_35 <- new_WW[which(new_WW$CHROM == "35"),]
WW_36 <- new_WW[which(new_WW$CHROM == "36"),]
WW_37 <- new_WW[which(new_WW$CHROM == "37"),]
WW_38 <- new_WW[which(new_WW$CHROM == "38"),]
WW_39 <- new_WW[which(new_WW$CHROM == "39"),]
WW_40 <- new_WW[which(new_WW$CHROM == "40"),]
WW_41 <- new_WW[which(new_WW$CHROM == "41"),]
WW_42 <- new_WW[which(new_WW$CHROM == "42"),]

TR_1 <- new_TR[which(new_TR$CHROM == "1"),]
TR_2 <- new_TR[which(new_TR$CHROM == "2"),]
TR_3 <- new_TR[which(new_TR$CHROM == "3"),]
TR_4 <- new_TR[which(new_TR$CHROM == "4"),]
TR_5 <- new_TR[which(new_TR$CHROM == "5"),]
TR_6 <- new_TR[which(new_TR$CHROM == "6"),]
TR_7 <- new_TR[which(new_TR$CHROM == "7"),]
TR_8 <- new_TR[which(new_TR$CHROM == "8"),]
TR_9 <- new_TR[which(new_TR$CHROM == "9"),]
TR_10 <- new_TR[which(new_TR$CHROM == "10"),]
TR_11 <- new_TR[which(new_TR$CHROM == "11"),]
TR_12 <- new_TR[which(new_TR$CHROM == "12"),]
TR_13 <- new_TR[which(new_TR$CHROM == "13"),]
TR_14 <- new_TR[which(new_TR$CHROM == "14"),]
TR_15 <- new_TR[which(new_TR$CHROM == "15"),]
TR_16 <- new_TR[which(new_TR$CHROM == "16"),]
TR_17 <- new_TR[which(new_TR$CHROM == "17"),]
TR_18 <- new_TR[which(new_TR$CHROM == "18"),]
TR_19 <- new_TR[which(new_TR$CHROM == "19"),]
TR_20 <- new_TR[which(new_TR$CHROM == "20"),]
TR_21 <- new_TR[which(new_TR$CHROM == "21"),]
TR_22 <- new_TR[which(new_TR$CHROM == "22"),]
TR_23 <- new_TR[which(new_TR$CHROM == "23"),]
TR_24 <- new_TR[which(new_TR$CHROM == "24"),]
TR_25 <- new_TR[which(new_TR$CHROM == "25"),]
TR_26 <- new_TR[which(new_TR$CHROM == "26"),]
TR_27 <- new_TR[which(new_TR$CHROM == "27"),]
TR_28 <- new_TR[which(new_TR$CHROM == "28"),]
TR_29 <- new_TR[which(new_TR$CHROM == "29"),]
TR_30 <- new_TR[which(new_TR$CHROM == "30"),]
TR_31 <- new_TR[which(new_TR$CHROM == "31"),]
TR_32 <- new_TR[which(new_TR$CHROM == "32"),]
TR_33 <- new_TR[which(new_TR$CHROM == "33"),]
TR_34 <- new_TR[which(new_TR$CHROM == "34"),]
TR_35 <- new_TR[which(new_TR$CHROM == "35"),]
TR_36 <- new_TR[which(new_TR$CHROM == "36"),]
TR_37 <- new_TR[which(new_TR$CHROM == "37"),]
TR_38 <- new_TR[which(new_TR$CHROM == "38"),]
TR_39 <- new_TR[which(new_TR$CHROM == "39"),]
TR_40 <- new_TR[which(new_TR$CHROM == "40"),]
TR_41 <- new_TR[which(new_TR$CHROM == "41"),]
TR_42 <- new_TR[which(new_TR$CHROM == "42"),]

EC_1 <- new_EC[which(new_EC$CHROM == "1"),]
EC_2 <- new_EC[which(new_EC$CHROM == "2"),]
EC_3 <- new_EC[which(new_EC$CHROM == "3"),]
EC_4 <- new_EC[which(new_EC$CHROM == "4"),]
EC_5 <- new_EC[which(new_EC$CHROM == "5"),]
EC_6 <- new_EC[which(new_EC$CHROM == "6"),]
EC_7 <- new_EC[which(new_EC$CHROM == "7"),]
EC_8 <- new_EC[which(new_EC$CHROM == "8"),]
EC_9 <- new_EC[which(new_EC$CHROM == "9"),]
EC_10 <- new_EC[which(new_EC$CHROM == "10"),]
EC_11 <- new_EC[which(new_EC$CHROM == "11"),]
EC_12 <- new_EC[which(new_EC$CHROM == "12"),]
EC_13 <- new_EC[which(new_EC$CHROM == "13"),]
EC_14 <- new_EC[which(new_EC$CHROM == "14"),]
EC_15 <- new_EC[which(new_EC$CHROM == "15"),]
EC_16 <- new_EC[which(new_EC$CHROM == "16"),]
EC_17 <- new_EC[which(new_EC$CHROM == "17"),]
EC_18 <- new_EC[which(new_EC$CHROM == "18"),]
EC_19 <- new_EC[which(new_EC$CHROM == "19"),]
EC_20 <- new_EC[which(new_EC$CHROM == "20"),]
EC_21 <- new_EC[which(new_EC$CHROM == "21"),]
EC_22 <- new_EC[which(new_EC$CHROM == "22"),]
EC_23 <- new_EC[which(new_EC$CHROM == "23"),]
EC_24 <- new_EC[which(new_EC$CHROM == "24"),]
EC_25 <- new_EC[which(new_EC$CHROM == "25"),]
EC_26 <- new_EC[which(new_EC$CHROM == "26"),]
EC_27 <- new_EC[which(new_EC$CHROM == "27"),]
EC_28 <- new_EC[which(new_EC$CHROM == "28"),]
EC_29 <- new_EC[which(new_EC$CHROM == "29"),]
EC_30 <- new_EC[which(new_EC$CHROM == "30"),]
EC_31 <- new_EC[which(new_EC$CHROM == "31"),]
EC_32 <- new_EC[which(new_EC$CHROM == "32"),]
EC_33 <- new_EC[which(new_EC$CHROM == "33"),]
EC_34 <- new_EC[which(new_EC$CHROM == "34"),]
EC_35 <- new_EC[which(new_EC$CHROM == "35"),]
EC_36 <- new_EC[which(new_EC$CHROM == "36"),]
EC_37 <- new_EC[which(new_EC$CHROM == "37"),]
EC_38 <- new_EC[which(new_EC$CHROM == "38"),]
EC_39 <- new_EC[which(new_EC$CHROM == "39"),]
EC_40 <- new_EC[which(new_EC$CHROM == "40"),]
EC_41 <- new_EC[which(new_EC$CHROM == "41"),]
EC_42 <- new_EC[which(new_EC$CHROM == "42"),]

#6. 合并
TR_EC_1 <- merge(TR_1, EC_1, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_1 <- merge(TR_EC_1,WW_1,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_1[is.na(TR_EC_WW_1)] <- "./."
write.table(TR_EC_WW_1,"TR_EC_WW_1.txt",row.names=F,sep="\t",quote=F)

TR_EC_2 <- merge(TR_2, EC_2, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_2 <- merge(TR_EC_2,WW_2,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_2[is.na(TR_EC_WW_2)] <- "./."
write.table(TR_EC_WW_2,"TR_EC_WW_2.txt",row.names=F,sep="\t",quote=F)

TR_EC_3 <- merge(TR_3, EC_3, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_3 <- merge(TR_EC_3,WW_3,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_3[is.na(TR_EC_WW_3)] <- "./."
write.table(TR_EC_WW_3,"TR_EC_WW_3.txt",row.names=F,sep="\t",quote=F)

TR_EC_4 <- merge(TR_4, EC_4, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_4 <- merge(TR_EC_4,WW_4,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_4[is.na(TR_EC_WW_4)] <- "./."
write.table(TR_EC_WW_4,"TR_EC_WW_4.txt",row.names=F,sep="\t",quote=F)

TR_EC_5 <- merge(TR_5, EC_5, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_5 <- merge(TR_EC_5,WW_5,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_5[is.na(TR_EC_WW_5)] <- "./."
write.table(TR_EC_WW_5,"TR_EC_WW_5.txt",row.names=F,sep="\t",quote=F)

TR_EC_6 <- merge(TR_6, EC_6, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_6 <- merge(TR_EC_6,WW_6,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_6[is.na(TR_EC_WW_6)] <- "./."
write.table(TR_EC_WW_6,"TR_EC_WW_6.txt",row.names=F,sep="\t",quote=F)

TR_EC_7 <- merge(TR_7, EC_7, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_7 <- merge(TR_EC_7,WW_7,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_7[is.na(TR_EC_WW_7)] <- "./."
write.table(TR_EC_WW_7,"TR_EC_WW_7.txt",row.names=F,sep="\t",quote=F)

TR_EC_8 <- merge(TR_8, EC_8, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_8 <- merge(TR_EC_8,WW_8,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_8[is.na(TR_EC_WW_8)] <- "./."
write.table(TR_EC_WW_8,"TR_EC_WW_8.txt",row.names=F,sep="\t",quote=F)

TR_EC_9 <- merge(TR_9, EC_9, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_9 <- merge(TR_EC_9,WW_9,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_9[is.na(TR_EC_WW_9)] <- "./."
write.table(TR_EC_WW_9,"TR_EC_WW_9.txt",row.names=F,sep="\t",quote=F)

TR_EC_10 <- merge(TR_10, EC_10, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_10 <- merge(TR_EC_10,WW_10,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_10[is.na(TR_EC_WW_10)] <- "./."
write.table(TR_EC_WW_10,"TR_EC_WW_10.txt",row.names=F,sep="\t",quote=F)

TR_EC_11 <- merge(TR_11, EC_11, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_11 <- merge(TR_EC_11,WW_11,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_11[is.na(TR_EC_WW_11)] <- "./."
write.table(TR_EC_WW_11,"TR_EC_WW_11.txt",row.names=F,sep="\t",quote=F)

TR_EC_12 <- merge(TR_12, EC_12, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_12 <- merge(TR_EC_12,WW_12,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_12[is.na(TR_EC_WW_12)] <- "./."
write.table(TR_EC_WW_12,"TR_EC_WW_12.txt",row.names=F,sep="\t",quote=F)

TR_EC_13 <- merge(TR_13, EC_13, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_13 <- merge(TR_EC_13,WW_13,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_13[is.na(TR_EC_WW_13)] <- "./."
write.table(TR_EC_WW_13,"TR_EC_WW_13.txt",row.names=F,sep="\t",quote=F)

TR_EC_14 <- merge(TR_14, EC_14, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_14 <- merge(TR_EC_14,WW_14,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_14[is.na(TR_EC_WW_14)] <- "./."
write.table(TR_EC_WW_14,"TR_EC_WW_14.txt",row.names=F,sep="\t",quote=F)

TR_EC_15 <- merge(TR_15, EC_15, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_15 <- merge(TR_EC_15,WW_15,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_15[is.na(TR_EC_WW_15)] <- "./."
write.table(TR_EC_WW_15,"TR_EC_WW_15.txt",row.names=F,sep="\t",quote=F)

TR_EC_16 <- merge(TR_16, EC_16, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_16 <- merge(TR_EC_16,WW_16,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_16[is.na(TR_EC_WW_16)] <- "./."
write.table(TR_EC_WW_16,"TR_EC_WW_16.txt",row.names=F,sep="\t",quote=F)

TR_EC_17 <- merge(TR_17, EC_17, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_17 <- merge(TR_EC_17,WW_17,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_17[is.na(TR_EC_WW_17)] <- "./."
write.table(TR_EC_WW_17,"TR_EC_WW_17.txt",row.names=F,sep="\t",quote=F)

TR_EC_18 <- merge(TR_18, EC_18, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_18 <- merge(TR_EC_18,WW_18,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_18[is.na(TR_EC_WW_18)] <- "./."
write.table(TR_EC_WW_18,"TR_EC_WW_18.txt",row.names=F,sep="\t",quote=F)

TR_EC_19 <- merge(TR_19, EC_19, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_19 <- merge(TR_EC_19,WW_19,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_19[is.na(TR_EC_WW_19)] <- "./."
write.table(TR_EC_WW_19,"TR_EC_WW_19.txt",row.names=F,sep="\t",quote=F)

TR_EC_20 <- merge(TR_20, EC_20, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_20 <- merge(TR_EC_20,WW_20,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_20[is.na(TR_EC_WW_20)] <- "./."
write.table(TR_EC_WW_20,"TR_EC_WW_20.txt",row.names=F,sep="\t",quote=F)

TR_EC_21 <- merge(TR_21, EC_21, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_21 <- merge(TR_EC_21,WW_21,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_21[is.na(TR_EC_WW_21)] <- "./."
write.table(TR_EC_WW_21,"TR_EC_WW_21.txt",row.names=F,sep="\t",quote=F)

TR_EC_22 <- merge(TR_22, EC_22, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_22 <- merge(TR_EC_22,WW_22,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_22[is.na(TR_EC_WW_22)] <- "./."
write.table(TR_EC_WW_22,"TR_EC_WW_22.txt",row.names=F,sep="\t",quote=F)

TR_EC_23 <- merge(TR_23, EC_23, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_23 <- merge(TR_EC_23,WW_23,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_23[is.na(TR_EC_WW_23)] <- "./."
write.table(TR_EC_WW_23,"TR_EC_WW_23.txt",row.names=F,sep="\t",quote=F)

TR_EC_24 <- merge(TR_24, EC_24, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_24 <- merge(TR_EC_24,WW_24,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_24[is.na(TR_EC_WW_24)] <- "./."
write.table(TR_EC_WW_24,"TR_EC_WW_24.txt",row.names=F,sep="\t",quote=F)

TR_EC_25 <- merge(TR_25, EC_25, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_25 <- merge(TR_EC_25,WW_25,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_25[is.na(TR_EC_WW_25)] <- "./."
write.table(TR_EC_WW_25,"TR_EC_WW_25.txt",row.names=F,sep="\t",quote=F)

TR_EC_26 <- merge(TR_26, EC_26, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_26 <- merge(TR_EC_26,WW_26,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_26[is.na(TR_EC_WW_26)] <- "./."
write.table(TR_EC_WW_26,"TR_EC_WW_26.txt",row.names=F,sep="\t",quote=F)

TR_EC_27 <- merge(TR_27, EC_27, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_27 <- merge(TR_EC_27,WW_27,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_27[is.na(TR_EC_WW_27)] <- "./."
write.table(TR_EC_WW_27,"TR_EC_WW_27.txt",row.names=F,sep="\t",quote=F)

TR_EC_28 <- merge(TR_28, EC_28, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_28 <- merge(TR_EC_28,WW_28,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_28[is.na(TR_EC_WW_28)] <- "./."
write.table(TR_EC_WW_28,"TR_EC_WW_28.txt",row.names=F,sep="\t",quote=F)

TR_EC_29 <- merge(TR_29, EC_29, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_29 <- merge(TR_EC_29,WW_29,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_29[is.na(TR_EC_WW_29)] <- "./."
write.table(TR_EC_WW_29,"TR_EC_WW_29.txt",row.names=F,sep="\t",quote=F)

TR_EC_30 <- merge(TR_30, EC_30, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_30 <- merge(TR_EC_30,WW_30,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_30[is.na(TR_EC_WW_30)] <- "./."
write.table(TR_EC_WW_30,"TR_EC_WW_30.txt",row.names=F,sep="\t",quote=F)

TR_EC_31 <- merge(TR_31, EC_31, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_31 <- merge(TR_EC_31,WW_31,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_31[is.na(TR_EC_WW_31)] <- "./."
write.table(TR_EC_WW_31,"TR_EC_WW_31.txt",row.names=F,sep="\t",quote=F)

TR_EC_32 <- merge(TR_32, EC_32, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_32 <- merge(TR_EC_32,WW_32,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_32[is.na(TR_EC_WW_32)] <- "./."
write.table(TR_EC_WW_32,"TR_EC_WW_32.txt",row.names=F,sep="\t",quote=F)

TR_EC_33 <- merge(TR_33, EC_33, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_33 <- merge(TR_EC_33,WW_33,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_33[is.na(TR_EC_WW_33)] <- "./."
write.table(TR_EC_WW_33,"TR_EC_WW_33.txt",row.names=F,sep="\t",quote=F)

TR_EC_34 <- merge(TR_34, EC_34, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_34 <- merge(TR_EC_34,WW_34,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_34[is.na(TR_EC_WW_34)] <- "./."
write.table(TR_EC_WW_34,"TR_EC_WW_34.txt",row.names=F,sep="\t",quote=F)

TR_EC_35 <- merge(TR_35, EC_35, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_35 <- merge(TR_EC_35,WW_35,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_35[is.na(TR_EC_WW_35)] <- "./."
write.table(TR_EC_WW_35,"TR_EC_WW_35.txt",row.names=F,sep="\t",quote=F)

TR_EC_36 <- merge(TR_36, EC_36, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_36 <- merge(TR_EC_36,WW_36,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_36[is.na(TR_EC_WW_36)] <- "./."
write.table(TR_EC_WW_36,"TR_EC_WW_36.txt",row.names=F,sep="\t",quote=F)

TR_EC_37 <- merge(TR_37, EC_37, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_37 <- merge(TR_EC_37,WW_37,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_37[is.na(TR_EC_WW_37)] <- "./."
write.table(TR_EC_WW_37,"TR_EC_WW_37.txt",row.names=F,sep="\t",quote=F)

TR_EC_38 <- merge(TR_38, EC_38, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_38 <- merge(TR_EC_38,WW_38,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_38[is.na(TR_EC_WW_38)] <- "./."
write.table(TR_EC_WW_38,"TR_EC_WW_38.txt",row.names=F,sep="\t",quote=F)

TR_EC_39 <- merge(TR_39, EC_39, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_39 <- merge(TR_EC_39,WW_39,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_39[is.na(TR_EC_WW_39)] <- "./."
write.table(TR_EC_WW_39,"TR_EC_WW_39.txt",row.names=F,sep="\t",quote=F)

TR_EC_40 <- merge(TR_40, EC_40, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_40 <- merge(TR_EC_40,WW_40,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_40[is.na(TR_EC_WW_40)] <- "./."
write.table(TR_EC_WW_40,"TR_EC_WW_40.txt",row.names=F,sep="\t",quote=F)

TR_EC_41 <- merge(TR_41, EC_41, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_41 <- merge(TR_EC_41,WW_41,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_41[is.na(TR_EC_WW_41)] <- "./."
write.table(TR_EC_WW_41,"TR_EC_WW_41.txt",row.names=F,sep="\t",quote=F)

TR_EC_42 <- merge(TR_42, EC_42, by="ID", all=TRUE)[,-c(2:9,38:45)]
TR_EC_WW_42 <- merge(TR_EC_42,WW_42,by="ID",all=TRUE)[,-c(116:123)]
TR_EC_WW_42[is.na(TR_EC_WW_42)] <- "./."
write.table(TR_EC_WW_42,"TR_EC_WW_42.txt",row.names=F,sep="\t",quote=F)