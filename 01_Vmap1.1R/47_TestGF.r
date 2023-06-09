library(gradientForest)
library(RColorBrewer)
library(rasterVis)
library(LEA)
library(raster)
library(adegenet)
library(maps)
library(dismo)
library(gplots)
library(gdistance)

gfData <- read.csv("/data1/home/yafei/003_Project3/Structure/gradientForest/Climate-Change-Genomics-main/datasets/input/gradient_forest_input.csv",sep="\t",row.names=1) 
dim(gfData) #这个文件一共使用了13个群体，X指的是经度，Y指的是纬度，bio_1-bio19指的是环境因子，之后是555-19-2=534个SNP
#也就是说，这个文件一共使用了500个SNP，在小麦里面的话，可以每个亚基因挑200或者300个SNP进行计算
#可以对SNP进行分类，但是在这里是不需要的，因为这边没有需要分类的必要，有需要
candidate <- gfData
#把文件中的经度，纬度和环境变量提取出来
present <- gfData[,c(1,2,grep("bio",names(gfData)))]
#把19个环境变量重新命名，叫做bio_1-bio19
bioclimatic <- paste("bio_",1:19,sep = "")
#设置每个变量的排列分布的重要性。Type help(gradientForest for more details)
maxLevel <- log2(0.368*nrow(candidate)/2)
############################################现在开始run梯度森林
#要是第一次做需要Run，但是要是gf_runs.R这个文件产生之后(10M),就不需要再来这么一通了，因为有gf_runs.R这个中间文件出来
if(FALSE){ # FALSE if there is no need to run the analyses 
  gf_candidate <- gradientForest(cbind(present[,bioclimatic], candidate),
                                 predictor.vars=colnames(present[,bioclimatic]),
                                 response.vars=colnames(candidate), ntree=500,
                                 maxLevel=maxLevel, trace=T, corr.threshold=0.50)
  #将GF模型合并成一个列表并保存它
  gf_runs <- (gf_reference=gf_candidate)
  if(!dir.exists("/data1/home/yafei/003_Project3/Structure/gradientForest/Climate-Change-Genomics-main/datasets/output")){
    dir.create("/data1/home/yafei/003_Project3/Structure/gradientForest/Climate-Change-Genomics-main/datasets/output")
  }
  save(gf_runs,file = "/data1/home/yafei/003_Project3/Structure/gradientForest/Climate-Change-Genomics-main/datasets/output/gf_runs.R")
}
#####要是之前run过，可以在这一步直接load
load("/data1/home/yafei/003_Project3/Structure/gradientForest/Climate-Change-Genomics-main/datasets/output/gf_runs.R")
gf_candidate <- gf_runs$gf_candidate

##################################计算genetic offset
#建了当前和未来气候数据的矩阵，这些矩阵将用于外推梯度森林分析在整个地理景观中构建的功能 
#创建一个raster layer，将值设为零，但我们使用与生物气候变量相同的范围和分辨率。
mask <- raster("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/present/bio_1.asc") %>% replace(.,.,0)
#使用了convert_env_trns这个函数(需要提前source或者是跑一下这个climate_change_functions.R里面的这个function)
#该函数将任何给定数量的enviromental raster layers转换为data frame
turn_score <- data.frame(gfData[,c("X","Y",most_cand)],temp) #这里是得到经纬度，SNP解释的最多环境变量，还有impact的值提取出来
#导入环境变化的信息，这样的话就可以实现现在，2050，2070年的环境预测
#2050-RCP4.5，不是太剧烈；2070-RCP8.5，这样就剧烈很多
present_mat <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/present/")
future_mat_50 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/year_50/") 
future_mat_70 <- convert_env_trns(path = "/Users/xuebozhao/Documents/LuLab/wheatSpeciation/SDM/year_70/") 
#预测整个地区的等位基因数据
pred_paSNPs <- predict(gf_candidate,present_mat[grep("bio",names(present_mat))])
pred_paSNPs_future_50 <- predict(gf_candidate,future_mat_50[grep("bio",names(future_mat_50))])
pred_paSNPs_future_70 <- predict(gf_candidate,future_mat_70[grep("bio",names(future_mat_70))])
#估计两个矩阵之间的欧氏距离，这是遗传偏移-genetic offset; 
#欧几里得距离函数也是一个climate_change_functions.R这个里面的附属函数
euclidian_50 <-euclidian_distance(proj_fut=pred_paSNPs_future_50,pred_pres=pred_paSNPs)
euclidian_70 <-euclidian_distance(proj_fut=pred_paSNPs_future_70,pred_pres=pred_paSNPs) 
#创建一个raster layer，这里面包含每一个像素的genetic offset的值，之后我们往这个raster layer添加信息
offset_ras_50 <- mask
offset_ras_50[present_mat$cell]<- euclidian_50 #the present_mat$cell里面的每个cell包含每个种群的像素
offset_ras_70 <- mask
offset_ras_70[present_mat$cell]<- euclidian_70
#获得已知群体的genetic offset
genetic_off_50 <- raster::extract(offset_ras_50,present[,1:2])
genetic_off_70 <- raster::extract(offset_ras_70,present[,1:2])
#创建包含种群坐标和genetic offset的表
pop_vul_50 <- data.frame(present[,1:2], genetic_off_50)
pop_vul_70 <- data.frame(present[,1:2], genetic_off_70)
#将genetic offset的类别，就是之前分出来的冷和热，添加到表格中，以便将来进行分析。 
pop_vul_50$temp <- NA
pop_vul_50$temp[warm]<-"warm"
pop_vul_50$temp[cold]<-"cold"
pop_vul_50$temp <- factor(pop_vul_50$temp,levels = c("cold","warm"))
pop_vul_70$temp <- NA
pop_vul_70$temp[warm]<-"warm"
pop_vul_70$temp[cold]<-"cold"
pop_vul_70$temp <- factor(pop_vul_70$temp,levels = c("cold","warm"))
#得到最大的max offset
max_val <- max(c(pop_vul_50$genetic_off_50,pop_vul_70$genetic_off_70))
#创建一个表格，显示地图中每个cell的遗传偏移量
#2050
vulnerable_areas <- as.data.frame(offset_ras_50)
vulnerable_areas <- data.frame(cell=1:nrow(vulnerable_areas),vul=vulnerable_areas[,1])
vulnerable_areas_50 <- vulnerable_areas[which(vulnerable_areas$vul>0),] #把可以预测的vulnerability的cell提取出来
#2070
vulnerable_areas <- as.data.frame(offset_ras_70)
vulnerable_areas <- data.frame(cell=1:nrow(vulnerable_areas),vul=vulnerable_areas[,1])
vulnerable_areas_70 <- vulnerable_areas[which(vulnerable_areas$vul>0),]
gen_off_stack <- stack(offset_ras_50,offset_ras_70)
names(gen_off_stack) <- paste(c("Year_2050","Year_2070"))
#plot genetic offset and the known populations according to the environmental cluster in which they grow
pal=colorRampPalette(c("#4575B4","#74ADD1","#E0F3F8","white","#FEE090","#F46D43","#D73027"))
rasterVis::levelplot(gen_off_stack,margin=FALSE, colorkey=list(space="bottom"),xlab=NULL, ylab=NULL, 
                     scales=list(draw=FALSE), main = "Genomic offset",
                     col.regions=pal)