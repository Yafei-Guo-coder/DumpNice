#循环：for
v <- LETTERS[1:6]
for ( i in v) {
  if (i == "D") {  # D 不会输出，跳过这次循环，进入下一次
    next
  }
  print(i)
}
#连接两个字符串：paste
paste (a, sep = " ", collapse = NULL)
paste0(a, collapse = NULL)
#定义一个空向量
x=vector()
x<-numeirc(0) #长度可变的存储数字的向量
x=character() #创建出来的为字符串向量
x[1]=1
x[2]=3
x<-NULL; x[1]<-2
vector(mode="numeric",length=0) #定义一个空向量，往里面添加元素即可
x=matrix(nrow = 2,ncol=3) #创建空矩阵
#添加元素
c1 <- c(1,2,3,4,5) #创建一个向量
c1 <- c(c1,5) #追加一个元素
c1 <- c(c1,c(5,6)) #追加一个向量
c1 <- c(c1[1:2],c(5,6),c[3:5]) #指定位置来添加的元素
c1 <- append(c1,8) #在向量最后追加一个元素8
c1 <- append(c1,c(11,22)) #在向量后追加向量
c1 <- append(c1,35,3) #在第3个元素后插入新元素，也可以插入向量

#删除元素
c1 <- c1[-1] #从向量中指定位置为1的元素
c1 <- c1[-c(2:3)] #可以给定一个位置向量来删除多个元素
c1 <- c1[c(3:5)] #与上面的方式相反，保留想要的元素
#展示调色板
library(RColorBrewer)
display.brewer.all()
brewer.pal(9, "Set3")[c(9,7,2,3)]
scale_fill_brewer(palette = "Accent")
scale_fill_manual(values = c("#FDC086","#BEAED4")) 


