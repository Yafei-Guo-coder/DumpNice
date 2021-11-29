library(raster)
library(rgdal)
library(rasterVis)
library(RColorBrewer)
library(ggplot2)
setwd("/Volumes/HP x750w/wc2.1_30s_bio")
setwd("/Volumes/HP x750w/wc2.1_2.5m_tmax_2010-2018")
temp1 <- raster("wc2.1_30s_bio_1.tif")
temp2 <- raster("wc2.1_30s_bio_2.tif")
temp3 <- raster("wc2.1_30s_bio_3.tif")
temp4 <- raster("wc2.1_30s_bio_4.tif")
temp5 <- raster("wc2.1_30s_bio_5.tif")
temp6 <- raster("wc2.1_30s_bio_6.tif")
temp7 <- raster("wc2.1_30s_bio_7.tif")
temp8 <- raster("wc2.1_30s_bio_8.tif")
temp9 <- raster("wc2.1_30s_bio_9.tif")
temp10 <- raster("wc2.1_30s_bio_10.tif")
temp11 <- raster("wc2.1_30s_bio_11.tif")
temp12 <- raster("wc2.1_30s_bio_12.tif")
temp13 <- raster("wc2.1_30s_bio_13.tif")
temp14 <- raster("wc2.1_30s_bio_14.tif")
temp15 <- raster("wc2.1_30s_bio_15.tif")
temp16 <- raster("wc2.1_30s_bio_16.tif")
temp17 <- raster("wc2.1_30s_bio_17.tif")
temp18 <- raster("wc2.1_30s_bio_18.tif")
temp19 <- raster("wc2.1_30s_bio_19.tif")
setwd("/Volumes/HP x750w/wc2.1_30s_elev")
elev <- raster("/Volumes/HP x750w/wc2.1_30s_elev/wc2.1_30s_elev.tif")
setwd("/Volumes/HP x750w/wc2.1_30s_prec")
prec20 <- raster("wc2.1_30s_prec_01.tif")
prec21 <- raster("wc2.1_30s_prec_02.tif")
prec22 <- raster("wc2.1_30s_prec_03.tif")
prec23 <- raster("wc2.1_30s_prec_04.tif")
prec24 <- raster("wc2.1_30s_prec_05.tif")
prec25 <- raster("wc2.1_30s_prec_06.tif")
prec26 <- raster("wc2.1_30s_prec_07.tif")
prec27 <- raster("wc2.1_30s_prec_08.tif")
prec28 <- raster("wc2.1_30s_prec_09.tif")
prec29 <- raster("wc2.1_30s_prec_10.tif")
prec30 <- raster("wc2.1_30s_prec_11.tif")
prec31 <- raster("wc2.1_30s_prec_12.tif")
setwd("/Volumes/HP x750w/wc2.1_30s_tavg")
tavg32 <- raster("wc2.1_30s_tavg_01.tif")
tavg33 <- raster("wc2.1_30s_tavg_02.tif")
tavg34 <- raster("wc2.1_30s_tavg_03.tif")
tavg35 <- raster("wc2.1_30s_tavg_04.tif")
tavg36 <- raster("wc2.1_30s_tavg_05.tif")
tavg37 <- raster("wc2.1_30s_tavg_06.tif")
tavg38 <- raster("wc2.1_30s_tavg_07.tif")
tavg39 <- raster("wc2.1_30s_tavg_08.tif")
tavg40 <- raster("wc2.1_30s_tavg_09.tif")
tavg41 <- raster("wc2.1_30s_tavg_10.tif")
tavg42 <- raster("wc2.1_30s_tavg_11.tif")
tavg43 <- raster("wc2.1_30s_tavg_12.tif")

taxa <- read.table("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/select_taxa4.txt", sep="\t",header=T, stringsAsFactors = F)
lon <- as.numeric(taxa$Logititude)
lat <- as.numeric(taxa$Latitude)
samples <- data.frame(lon,lat)
temp.data <- samples
temp.data$elev <- extract(elev, samples)
temp.data$temp1 <- extract(temp1, samples)
temp.data$temp2 <- extract(temp2, samples)
temp.data$temp3 <- extract(temp3, samples)
temp.data$temp4 <- extract(temp4, samples)
temp.data$temp5 <- extract(temp5, samples)
temp.data$temp6 <- extract(temp6, samples)
temp.data$temp7 <- extract(temp7, samples)
temp.data$temp8 <- extract(temp8, samples)
temp.data$temp9 <- extract(temp9, samples)
temp.data$temp10 <- extract(temp10, samples)
temp.data$temp11 <- extract(temp11, samples) 
temp.data$prec1 <- extract(temp12, samples)
temp.data$prec2 <- extract(temp13, samples)
temp.data$prec3 <- extract(temp14, samples)
temp.data$prec4 <- extract(temp15, samples)
temp.data$prec5 <- extract(temp16, samples)
temp.data$prec6 <- extract(temp17, samples)
temp.data$prec7 <- extract(temp18, samples)
temp.data$prec8 <- extract(temp19, samples)
rownames(temp.data) <- taxa[,1]
all <- temp.data[,-c(1,2)]
setwd("/Users/guoyafei/Documents/01_个人项目/02_Migration/02_数据表格/01_Vmap1-1/01_Add_ZNdata/05_Environment/")
write.table(all, "select_bio2.txt", sep="\t", row.names = T,quote=F)




temp6_2018 <- raster("wc2.1_2.5m_tmax_2018-06.tif")
temp6_2017 <- raster("wc2.1_2.5m_tmax_2017-06.tif")
temp6_2016 <- raster("wc2.1_2.5m_tmax_2016-06.tif")
temp6_2015 <- raster("wc2.1_2.5m_tmax_2015-06.tif")
temp6_2014 <- raster("wc2.1_2.5m_tmax_2014-06.tif")
temp6_2013 <- raster("wc2.1_2.5m_tmax_2013-06.tif")

colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))
oregon <- readOGR('.', 'Oregon_10N')
levelplot(temp6_2015, 
          margin=FALSE,# suppress marginal graphics
          ylim=c(-50,60),
          colorkey=list(
            space='bottom',                   # plot legend at bottom
            labels=list(at=-5:5, font=4)      # legend ticks and labels 
          ),    
          par.settings=list(
            axis.line=list(col='transparent') # suppress axes and legend outline
          ),
          scales=list(draw=FALSE),            # suppress axis labels
          col.regions=colr,                   # colour ramp
          at=seq(-10, 50)) 
  #layer(sp.polygons(oregon, lwd=3))           # add oregon SPDF with latticeExtra::layer

