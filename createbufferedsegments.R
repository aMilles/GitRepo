library(rgdal)
library(ggplot2)
library(ggmap)
library(MBA)
library(geosphere)
library(sp)
library(raster)

rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- readOGR("Z:/GEC/segments_GEE.shp")

centroids <- rgeos::gCentroid(segments, byid = T)
centroids_102024 <- spTransform(centroids, "+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
buffered_centroids_102024 <- rgeos::gBuffer(centroids_102024, width = 5000)

buffered_centroids <- spTransform(buffered_centroids_102024, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
bc <- as(buffered_centroids, "SpatialPolygonsDataFrame")

writeOGR(bc, "Z:/GEC/buffered_segments.shp", "buffered_segments.shp", "ESRI Shapefile")