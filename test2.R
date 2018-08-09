xy <- read.csv("Z:/modelling/yxtable.csv")[1:100, ]

write.csv(xy, "xy.csv")

segments <- rgdal::readOGR("Z:/GEC/segments_GEE.shp")

rgdal::writeOGR(segments[match(xy$ID, segments$ID),], "segments.shp", "segments.shp", "ESRI Shapefile")
      
                