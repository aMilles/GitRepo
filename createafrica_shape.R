library(raster)
library(rgdal)
library(rmapshaper)

#CREATE A SIMPLE SHAPEFILE FOR AFRICA AND REMOVE SMALL ISLANDS

cc<-ccodes()
cc<-cc[which(cc$CONTINENT == "Africa" & !cc$SOVEREIGN %in% c("France", "United Kingdom")),]
for(country in cc$ISO3) assign(paste0("shp_", country), getData(name = "GADM", country = country, level = 0))

all.countries <- mget(x = ls(pattern = "shp_"))
africa <- do.call("rbind", all.countries)
africa2 <- rmapshaper::ms_dissolve(africa)
africa3 <- rmapshaper::ms_filter_islands(africa2, min_area = 2e+13)

#only outline of the continent
writeOGR(africa3, dsn = "africa", layer = "africa_continent", driver = "ESRI Shapefile")

#borders of countries included
writeOGR(africa, dsn = "africa", layer = "africa_countries", driver = "ESRI Shapefile")