setwd("/Users/guoyafei/Desktop/")
data <- read.table("Barcode",header=F,stringsAsFactors = F)
data <- as.data.frame.matrix(data)
mylist <- vector( mode = "list" )
for(i in 1:108){ 
  if(length(strsplit(data[i,1],"")[[1]])==4){
    for(j in 1:2){
      print(j)
      print(substr(data[i,1], j, j+2))
    }}else if(length(strsplit(data[i,1],"")[[1]])==5){
      for(j in 1:2){
        print(j)
        print(substr(data[i,1], j, j+3))
      }
    } else if(length(strsplit(data[i,1],"")[[1]])==6){
      for(j in 1:2){
        print(j)
        print(substr(data[i,1], j, j+4))
      }} else {
        for(j in 1:2){
          print(j)
          print(substr(data[i,1], j, j+5))
        }
      }
}
