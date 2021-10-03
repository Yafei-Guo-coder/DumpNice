setwd("/Users/guoyafei/Downloads/Results/")
W.pos <- read.table("World.pos",header=T,stringsAsFactors = F)
B.pos <- read.table("sorted.pos",header=F,stringsAsFactors = F)
colnames(B.pos) <- c("psid", "alle1", "alle2")
World <- merge(W.pos, B.pos, by="psId",all.x=TRUE)
rm(W.pos)

Trancing <- read.table("Trancing.SNPs.pos",header = T,stringsAsFactors = F)
Exome <- read.table("Exome.pos",header=F,stringsAsFactors = F)

W <- paste(World[,2],World[,3],sep="")
T <- paste(Trancing[,1],Trancing[,2],sep="")
E <- paste(Exome[,1],Exome[,2],sep="")

newW <- cbind(W,World[,1],World[,4],World[,5])
newW <- as.data.frame(newW)
colnames(newW) <- c("chr","psid","alleWA","alleWB")
newE <- cbind(E,Exome[,3],Exome[,4])
newE <- as.data.frame(newE)
colnames(newE) <- c("chr","alleEA","alleEB")
newT <- cbind(T,Trancing[,3],Trancing[,4])
newT <- as.data.frame(newT)
colnames(newT) <- c("chr","alleTA","alleTB")
rm(Exome)
rm(Trancing)

newWT <- merge(newW,newT,by="chr",all.x=TRUE)
newWTE <- merge(newWT, newE,by="chr",all.x=TRUE)

index <- apply(newWTE,1,function(x) length(which(is.na(x))))
new <- newWTE[which(index < 4),]
G <- apply(new,1,function(x) length(which(x=="G")))
A <- apply(new,1,function(x) length(which(x=="A")))
T <- apply(new,1,function(x) length(which(x=="T")))
C <- apply(new,1,function(x) length(which(x=="C")))

new$G <- G
new$A <- A
new$T <- T
new$C <- C

delete <- new[apply(new,1,function(x) length(which(x==0))!=2),]
select <- new[apply(new,1,function(x) length(which(x==0))==2),]

for (i in 1:14089) {
  if (is.na(select[i,5])){
    select[i,5] <- select[i,7]
    select[i,6] <- select[i,8]
  } 
}

rm(A)
rm(C)
rm(E)
rm(G)
rm(i)
rm(index)
rm(T)
rm(W)
rm(new)
rm(newE)
rm(newT)
rm(newW)
rm(newWT)
rm(newWTE)
rm(World)

newData <- select[,c(1,4,5)]

WorldwideResult <- read.table("WorldwideResult",header=T,stringsAsFactors = F)

deleteID <- as.character(delete[,2])
all <- WorldwideResult[,1]
B.ID <- B.pos[,1]

num1 <- which(all %in% deleteID)
WorldFilterd <- WorldwideResult[-num1,]

num2 <- which(B.ID %in% deleteID)
Bpos.filted <- B.pos[-num2,]


pos.merge <- merge(Bpos.filted, select, by="psid",all.x=TRUE)
pos.merge[,7] <- as.character(pos.merge[,7])
pos.merge[,8] <- as.character(pos.merge[,8])


for (i in 1:102690) {
  if (!is.na(pos.merge[i,7])){
    pos.merge[i,2] <- pos.merge[i,7]
    pos.merge[i,3] <- pos.merge[i,8]
  } 
}

pos <- pos.merge[,c(1,2,3)]

write.table(WorldFilterd, "WorldFilterdResult",quote = F,row.names = F)
write.table(pos, "WorldFilterd.pos",quote = F,row.names = F)

AA <- paste(pos[,2],pos[,2],sep="")
AB <- paste(pos[,2],pos[,3],sep="")
BB <- paste(pos[,3],pos[,3],sep="")

all <- cbind(AA,AB,BB,WorldFilterd)
write.table(all,file = "all.txt",quote = F,sep=" ",row.names = F)
