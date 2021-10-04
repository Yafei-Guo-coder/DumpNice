#######现在做的是migration的图
setwd("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems")
library(rEEMSplots)
extdata_path <- system.file("extdata", package = "rEEMSplots")
eems_results <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/test_out1")
name_figures <- file.path("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/test_figout","EEMS-example")
if (!file.exists(eems_results)) {
  stop("Check that the rEEMSplots package is installed without errors.")
}
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/test_figout/out1",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-default"),longlat = TRUE,out.png = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-axes-flipped"),longlat = FALSE,out.png = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-output-PNGs"),longlat = TRUE,plot.height = 8,plot.width = 7,res = 600,out.png = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-demes-and-edges"),longlat = TRUE,add.grid = TRUE,col.grid = "gray90",lwd.grid = 2,
           add.outline = TRUE,col.outline = "blue",lwd.outline = 5,add.demes = TRUE,col.demes = "red",pch.demes = 5,min.cex.demes = 0.5,max.cex.demes = 1.5,out.png = FALSE)
##
library(rgdal)
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"
pdf("/Users/xuebozhao/Documents/LuLab/wheatSpeciation/migration/eems/test_figout/out1",height = 6, width = 6, onefile = FALSE)
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-cartographic-projections"),
           longlat = TRUE,projection.in = projection_none,projection.out = projection_mercator,out.png = FALSE)
library("rworldmap")
library("rworldxtra")
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-geographic-map"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,add.map = TRUE,col.map = "black",lwd.map = 5,out.png = FALSE)
library("RColorBrewer")
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-new-eems-colors"),longlat = TRUE,
           eems.colors = brewer.pal(11, "RdBu"),out.png = FALSE)
##
library("rgdal") ## Defines functions to transform spatial elements
library("rworldmap") ## Defines world map
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"
map_world <- getMap() ## Add the map of Africa explicitly by passing the shape file
map_africa <- map_world[which(map_world@data$continent == "Africa"), ] 
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile"),longlat = TRUE, 
           m.plot.xy = { plot(map_africa, col = NA, add = TRUE) },q.plot.xy = { plot(map_africa, col = NA, add = TRUE) },out.png = FALSE)
## Apply the Mercator projection and add the map of Africa ## Don't forget to apply the same projection to the map as well
map_africa <- spTransform(map_africa, CRSobj = CRS(projection_mercator))
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-shapefile-projected"),longlat = TRUE,
           projection.in = projection_none,projection.out = projection_mercator,m.plot.xy = { plot(map_africa, col = NA, add = TRUE) },
           q.plot.xy = { plot(map_africa, col = NA, add = TRUE) },out.png = FALSE)
## Similarly we can add points, lines, labels, etc. ## Here is how to add a set of colored "labels" on top of the ## migration/diversity rates and the Africa map
coords <- matrix(c(-10, 10, 10, 10, 30, 0, 40, -10, 30, -20), ncol = 2, byrow = TRUE) 
colors <- c("red", "green", "blue", "purple", "orange") 
labels <- LETTERS[1:5]
coords_merc <- sp::spTransform(SpatialPoints(coords, CRS(projection_none)), CRS(projection_mercator))
## `coords_merc` is a SpatialPoints structure ## but we only need the coordinates themselves
coords_merc <- coords_merc@coords
eems.plots(mcmcpath = eems_results,plotpath = paste0(name_figures, "-labels-projected"),longlat = TRUE,projection.in = projection_none,projection.out = projection_mercator,
           m.plot.xy = { plot(map_africa, col = NA, add = TRUE);text(coords_merc, col = colors, pch = labels, font = 2); },
           q.plot.xy = { plot(map_africa, col = NA, add = TRUE);text(coords_merc, col = colors, pch = labels, font = 2); },out.png = FALSE)
dev.off()