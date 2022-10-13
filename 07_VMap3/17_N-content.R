setwd("/Users/guoyafei/Documents/02_VmapIII/02_DataSource/05_WEGA/")
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2)
library(ggmap)
library(RColorBrewer)

library(rasterVis)

data <- read.table("369sample.txt",header=T, stringsAsFactors = F,sep="\t")

N <- nc_open("/Users/guoyafei/公共资源/01_数据资源/TN5min.nc", write=FALSE, readunlim=TRUE, verbose=FALSE,
              auto_GMT=TRUE, suppress_dimvals=FALSE, return_on_error=FALSE )
lon <- ncvar_get(N, "lon")
lat <- ncvar_get(N,"lat",verbose=F)
t <- ncvar_get(N,"time")
ndvi.array <- ncvar_get(N,"TN")
dim(ndvi.array)
filevalue <- ncatt_get(N,"TN")
ndvi.array[ndvi.array==filevalue$value] <- NA
ndvi.slice <- ndvi.array[,,1]
dim(ndvi.slice)
r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
r <- raster(t(ndvi.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))

plot(r)

points(landrace$Longitude,landrace$Latitude,pch = 0, bg = "grey",cex=0.5)
points(cultivar$Longitude,cultivar$Latitude,pch = 8, bg = "grey",cex=0.5)

sample <- data[,c(9,12,11)]
sample <- sample[which(sample$Longitude != "NA"),]
landrace <- sample[which(sample$type.USDA_8. == "Landrace"),]
cultivar <- sample[which(sample$type.USDA_8. == "Cultivar"),]




sample <- data[,c(12,11)]
sample$NC <- extract(r, sample)
sample$source <- data$DataSource
sample <- sample[which(sample$NC != "NA"),]
#绘图
mp<-NULL
mapworld<-borders("world",colour = "gray70",fill="gray70") 
mp<-ggplot()+mapworld+ylim(-60,90) +theme_classic()
color <- brewer.pal(8, "Dark2")[c(1,2,3,4,6)]
#AABBDD----
mp+geom_point(aes(x=sample$Longitude, y=sample$Latitude, color=sample$NC,shape = sample$source),size=1,alpha=0.7)+
  scale_size(range=c(1,1))+ 
  #scale_color_manual(values = color) +
  scale_colour_gradient2(low = "dark red", mid = "white", high = "dark blue", midpoint = 15)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20))+
  #scale_fill_manual(values=c("#97FFFF", "#FFD700", "#FF6347", "#8470FF","#D8BFD8"), 
  #                   breaks=c("AA", "SS", "DD", "AABB","AABBDD"),
  #                   labels=c("AA", "SS", "DD", "AABB","AABBDD"))+
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank()) 
#theme(panel.grid =element_blank()) + 
#theme(axis.ticks.y = element_blank())




colr <- colorRampPalette(brewer.pal(11, 'RdYlBu')[3:9])
levelplot(r, 
          margin=FALSE,# suppress marginal graphics
          ylim=c(-50,80),
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





