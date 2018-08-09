library(RColorBrewer)
library(raster)
library(rgeos)
library(rgdal)
library(foreign)

setwd("Z:/")
shapes <- list.files(pattern = "WDPA_Apr2017-shapefile-polygons.shp", recursive = T)
shapes <- shapes[grep("WDPA", shapes)]

cc<-ccodes()
cc<-cc[which(cc$CONTINENT == "Africa" & !cc$SOVEREIGN %in% c("France", "United Kingdom")),]

a <- read.dbf("PA/WDPA_Apr2018-shapefile-polygons.dbf")

cc<-ccodes()
cc<-cc[which(cc$CONTINENT == "Africa" & !cc$SOVEREIGN %in% c("France", "United Kingdom")),]

WDPAID <- a$WDPAID[which(a$ISO3 %in% cc$ISO3)]

PA <- base::merge(cc, a)
