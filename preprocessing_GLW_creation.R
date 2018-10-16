library(RColorBrewer)
library(raster)
library(rgeos)
library(rgdal)
library(doParallel)
library(fasterize)
setwd("Z:/")

africa <- readOGR("Z:/africa/africa_continent.shp")

#read African cattle, goat and sheep data
for(tif in list.files(recursive = T, pattern = "1km_AD_2010_v2_1.tif$")) assign(basename(tif), raster(tif))

#stack the 3 datasets and multiply each density by the TLU (Tropical livestock unit) conversion factor
tlu_total <- stack(AF_Cattle1km_AD_2010_v2_1.tif * 0.7,
                   AF_Goats1km_AD_2010_v2_1.tif * 0.1,
                   AF_Sheep1km_AD_2010_v2_1.tif * 0.1)

#calculate the sum
tlu_total_single <- stackApply(tlu_total, indices = c(1,1,1), fun = sum)

#limit unreasonably high values (livestock densities > 400 units per / km² to 400)
tlu_total_single[tlu_total_single > 400] <- 400
tlu <- crop(tlu_total_single, africa)
africa <- fasterize(sf::st_as_sf(africa), tlu)
tlu <- africa * tlu

#save raster file (will be uploaded to GEE)
writeRaster(x = tlu, filename = "Z:/GLW/TLU.tif", overwrite = T)