###XPCLR的结果，焦雨铃老师
####这是画所有的 不是比较的*********************************************************
setwd("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/permutation/XPCLR_recom_rate/norm20M_10k")
cent = read.table("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatEvolution/pi/centromerer.txt",head =T)
col = rgb(105/255,105/255,105/255,alpha=0.8) 
xpxlrA = read.table("Chr_A.xpclr.txt")
chr = xpxlrA[,1]
chr.num <- unique(chr)
par(mfrow=c(7,1), oma=c(4.5,3, 0, 4), mar=c(0.3,1,0.3,1), cex=1)
twoPoints.x = c(0,1)
twoPoints.y = c(4.126591745,4.126591745)
for(c in chr.num){
  score = xpxlrA[chr == c,4]
  pos = xpxlrA[chr == c,2]
  index = cent$chr==c
  plot(pos/1000000,score,ylim = c(-2,50),axes = F,col = "blue",type="l",cex = 0.15,xlim=c(0,800))
  points(cent[index,2],0,pch=20,cex=2,col=col)
  twoPoints.x[2] = max(pos/1000000)
  points(twoPoints.x,twoPoints.y,type="l",col="black",lty=3,lwd=1.5)
  axis(2,cex.axis =0.6,at = c(0, 10, 20,30))
  #abline(h=4.126591745, lty=2, lwd=1, col="black")
  mtext(c,side = 3 ,las = 1, line  = -2.3 ,cex = 1,at = 380)
  if(c == "chr1A"){
    abline(v=505680794/1000000,col = "red",cex= 2)
  }
}
axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
mtext("XPCLR",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 210)
mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
##
mtext("Genome position",side = 1 ,las = 1, line  = 3 ,cex = 1,at = 450)
mtext("Wild einkorn - Domesticated einkorn",side = 3 ,las = 1, line  = 3 ,cex = 1,at = 450)

###########
xpxlrA = read.table("Chr_AB_15.xpclr.txt")
chr = xpxlrA[,1]
chr.num <- unique(chr)
par(mfrow=c(7,2), oma=c(4.5,4, 0, 0), mar=c(0.3,1,0.3,2), cex=1)
twoPoints.x = c(0,1)
twoPoints.y = c(2.776164597,2.776164597)
for(c in chr.num){
  score = xpxlrA[chr == c,4]
  pos = xpxlrA[chr == c,2]
  index = cent$chr==c
  plot(pos/1000000,score,ylim = c(-2,50),axes = F,col = "blue",type="l",cex = 0.15,xlim=c(0,800))
  points(cent[index,2],0,pch=20,cex=2,col=col)
  axis(2,cex.axis =0.6,at = c(0, 10, 20,30))
  twoPoints.x[2] = max(pos/1000000)
  points(twoPoints.x,twoPoints.y,type="l",col="black",lty=3,lwd=1.5)
  mtext(c,side = 3 ,las = 1, line  = -2.3 ,cex = 1,at = 380)
  if(c == "chr1A"){
    abline(v=505680794/1000000,col = "red",cex= 2)
    text(505680794/1000000, 36,expression(TraesCS1A02G314200),cex = 0.75)
  }
  if(c == "chr1B"){
    abline(v=551869725/1000000,col = "red",cex= 2)
    text((40+551869725/1000000), 20,expression(TraesCS1B02G326500),cex = 0.75)
  }
  if(c == "chr7A"){
    axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
    mtext("XPCLR",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 210)
  }
}
axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
#mtext("XPCLR",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 210)
mtext("Genome position",side = 1 ,las = 1, line  = 3.2 ,cex = 1,at = -150)
mtext("Wild emmer - Domesticated emmer",side = 3 ,las = 1, line  = 3 ,cex = 1,at = 450)

##################
xpxlrA = read.table("Chr_AB_45.xpclr.txt")
chr = xpxlrA[,1]
chr.num <- unique(chr)
par(mfrow=c(7,2), oma=c(4.5,4, 0, 0), mar=c(0.3,1,0.3,2), cex=1)
twoPoints.x = c(0,1)
twoPoints.y = c(3.203784603,3.203784603)
for(c in chr.num){
  score = xpxlrA[chr == c,4]
  pos = xpxlrA[chr == c,2]
  index = cent$chr==c
  plot(pos/1000000,score,ylim = c(-2,50),axes = F,col = "blue",type="l",cex = 0.15,xlim=c(0,800))
  points(cent[index,2],0,pch=20,cex=2,col=col)
  axis(2,cex.axis =0.6,at = c(0, 10, 20,30))
  twoPoints.x[2] = max(pos/1000000)
  points(twoPoints.x,twoPoints.y,type="l",col="black",lty=3,lwd=1.5)
  mtext(c,side = 3 ,las = 1, line  = -2.3 ,cex = 1,at = 380)
  if(c == "chr1A"){
    abline(v=505680794/1000000,col = "red",cex= 2)
    #text(505680794/1000000, 36,expression(TraesCS1A02G314200),cex = 0.75)
  }
  if(c == "chr1B"){
    abline(v=551869725/1000000,col = "red",cex= 2)
    #text((40+551869725/1000000), 20,expression(TraesCS1B02G326500),cex = 0.75)
  }
  if(c == "chr7A"){
    axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
    mtext("XPCLR",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 210)
  }
}
axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
mtext("Genome position",side = 1 ,las = 1, line  = 3.2 ,cex = 1,at = -150)
mtext("Domesticated emmer VS Durum",side = 1 ,las = 1, line  = 3.2 ,cex = 1,at = -150)

######
xpxlrA = read.table("Chr_ABD.xpclr.txt")
chr = xpxlrA[,1]
chr.num <- unique(chr)
par(mfrow=c(7,3), oma=c(4.5,4, 0, 0), mar=c(0.3,1,0.3,2), cex=1)
twoPoints.x = c(0,1)
twoPoints.y = c(4.126591745,4.126591745)
for(c in chr.num){
  score = xpxlrA[chr == c,4]
  pos = xpxlrA[chr == c,2]
  index = cent$chr==c
  plot(pos/1000000,score,ylim = c(-2,50),axes = F,col = "blue",type="l",cex = 0.15,xlim=c(0,800))
  points(cent[index,2],0,pch=20,cex=2,col=col)
  axis(2,cex.axis =0.6,at = c(0, 10, 20,30))
  twoPoints.x[2] = max(pos/1000000)
  points(twoPoints.x,twoPoints.y,type="l",col="black",lty=3,lwd=1.5)
  mtext(c,side = 3 ,las = 1, line  = -2.3 ,cex = 1,at = 340)
  if(c == "chr1A"){
    abline(v=505680794/1000000,col = "red",cex= 2)
    #text(505680794/1000000, 36,expression(TraesCS1A02G314200),cex = 0.75)
  }
  if(c == "chr1B"){
    abline(v=551869725/1000000,col = "red",cex= 2)
    #text((40+551869725/1000000), 20,expression(TraesCS1B02G326500),cex = 0.75)
  }
  if(c == "chr1D"){
    abline(v=409798580/1000000,col = "red",cex= 2)
    #text((40+551869725/1000000), 20,expression(TraesCS1B02G326500),cex = 0.75)
  }
  if(c == "chr7A"){
    axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
    mtext("XPCLR",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 210)
  }
  if(c == "chr7B"){
    axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
    mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
    mtext("XPCLR",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 210)
  }
  #num = num + sum(score > 5.89983,na.rm=T)
}
axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
mtext("Genome position",side = 1 ,las = 1, line  = 3.2 ,cex = 1,at = -150)
mtext("Bread wheat landrace - Bread wheat cultivar",side = 1 ,las = 1, line  = 3.2 ,cex = 1,at = -550)

#####
xpxlrA = read.table("Chr_D.xpclr.txt")
chr = xpxlrA[,1]
chr.num <- unique(chr)
par(mfrow=c(7,1), oma=c(4.5,3, 0, 4), mar=c(0.3,1,0.3,1), cex=1)
twoPoints.x = c(0,1)
twoPoints.y = c(4.126591745,4.126591745)
for(c in chr.num){
  score = xpxlrA[chr == c,4]
  pos = xpxlrA[chr == c,2]
  index = cent$chr==c
  plot(pos/1000000,score,ylim = c(-2,50),axes = F,col = "blue",type="l",cex = 0.15,xlim=c(0,800))
  points(cent[index,2],0,pch=20,cex=2,col=col)
  twoPoints.x[2] = max(pos/1000000)
  points(twoPoints.x,twoPoints.y,type="l",col="black",lty=3,lwd=1.5)
  axis(2,cex.axis =0.6,at = c(0, 10, 20,30))
  #abline(h=4.126591745, lty=2, lwd=1, col="black")
  mtext(c,side = 3 ,las = 1, line  = -2.3 ,cex = 1,at = 380)
  if(c == "chr1D"){
    abline(v=409798580/1000000,col = "red",cex= 2)
  }
}
axis(1, pos = -10,cex.axis =0.8,at = c(0, 200, 400, 600,800))
mtext("XPCLR",side = 2 ,las = 0, line  = 2.5 ,cex = 1,at = 210)
mtext("(Mb)",side = 1 ,las = 1, line  = 1.4 ,cex = 0.8,at = 850)
##
mtext("Genome position",side = 1 ,las = 1, line  = 3 ,cex = 1,at = 450)
mtext("Strangulata - Bread wheat landrace",side = 3 ,las = 1, line  = 3 ,cex = 1,at = 450)

