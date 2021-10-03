data <- read.table("/Users/guoyafei/Documents/Lulab/Project-4-VmapIII/test",header=T,stringsAsFactors = F)
library(spatstat)
library(maps)
library(maptools)
replyinsert_emoticon
more_vert
usdat <- data.frame(x=runif(50, -115, -85), y=runif(50, 33, 41), z=runif(50, 0, 100))
map('usa')
points(usdat$x, usdat$y)
Create an owin object from the USA outline(s)
usmap <- map('usa', fill=TRUE, col="transparent", plot=FALSE)
uspoly <- map2SpatialPolygons(usmap, IDs=usmap$names)
spatstat.options(checkpolygons=FALSE)
usowin <- as.owin.SpatialPolygons(uspoly)
spatstat.options(checkpolygons=TRUE)
Create a spatstat ppp object
pts <- as.ppp(usdat, W=usowin)
plot(pts)
Plot a a density surface
plot(density(pts))

