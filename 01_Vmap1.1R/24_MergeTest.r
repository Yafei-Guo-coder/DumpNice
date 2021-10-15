setwd("/Users/guoyafei/Downloads/Results/ChangeFormat/test_merge")
test.Exome <- read.table("test.Exome.vcf", head = FALSE, stringsAsFactors = FALSE)
test.Trancing <- read.table("test.Trancing.vcf", head = FALSE, stringsAsFactors = FALSE)
test.World <- read.table("test.World.vcf", head = FALSE, stringsAsFactors = FALSE)
test.merge <- read.table("test.merge.vcf", head = FALSE, stringsAsFactors = FALSE)

a <- paste(test.Trancing[,1],test.Trancing[,2],sep="_")
b <- paste(test.World[,1],test.World[,2],sep="_")
c  <- paste(test.merge[,1],test.merge[,2],sep="_")
d <- paste(test.Exome[,1],test.Exome[,2],sep="_")
venn.diagram(x= list(FileT = a, FileE= d,FileW = b,Merge = c), 
             filename = "merge.png", height = 450, width = 450,resolution =300, 
             imagetype="png", col ="transparent", 
             fill =c("cornflowerblue","green","yellow","darkorchid1"),
             alpha = 0.5, 
             label.col = c("orange", "white","darkorchid4", "white", "white", "white",    "white", "white","darkblue", "white", "white", "white","white", "darkgreen", "white"), 
             cex = 0.45,fontfamily = "serif", fontface = "bold",
             cat.col =c("darkblue", "darkgreen", "orange","darkorchid4"), 
             cat.cex = 0.45,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif", 
             rotation.degree = 270)

ab <- a[a%in%b]
ad <- a[a%in%d]
abd <- d[d%in%ab]
a[-which(a%in%c)]