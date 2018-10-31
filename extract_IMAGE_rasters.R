library(raster)
library(velox)
library(rgdal)
setwd("Z:/IMAGE/")

tif.files <- list.files(pattern = ".grd$")

names(stack(tif.files[1]))
# read segments.shp and create different buffer sizes for the extraction.
segments <- readOGR("Z:/GEC/segments_GEE.shp")
segments_sf <- sf::st_as_sf(segments)
segments_sf_05 <- sf::st_buffer(segments_sf, 0.5)
segments_sf_01 <- sf::st_buffer(segments_sf, 0.1)


#prepare and name the output list
image.ext <- vector("list", 3)
names(image.ext) <- c("y2015", "y2050", "y2100")
image.ext <- mget(rep("image.ext", length(tif.files)))
names(image.ext) <- tif.files
tif = 1
#scroll trough the singe raster stacks
for(tif in seq(length(tif.files))){
  print(tif)
  file <- list.files(pattern = ".grd$")[tif]
  r <- stack(file)
  
  #adjust buffer size of the segements according to the resolution of the tif
  if(res(r)[1] >= 0.5){
    seg <- segments_sf_05
  }else{
    seg <- segments_sf_01
  }   
  
  #scroll through the single years and extract all values (land cover) or calculate the mean
  for(year in 1:3){
    
    year_bands = 1:(length(names(r))/3)+ (length(names(r))/3 * (year - 1))
    print(year_bands)
    rb <- r[[year_bands]]
    vtif <- velox(rb)
    
    if(stringi::stri_count_regex(names(image.ext)[tif], "land_cover") > 0){
      image.ext[[tif]][[year]] <- (vtif$extract(seg))
      print(paste0("Image:", names(image.ext)[tif], " and band ", year, " use all values"))

    }else{
      image.ext[[tif]][[year]] <- data.frame(matrix(vtif$extract(seg, function(x) mean(x, na.rm = T)), nc = length(year_bands)))
      print(paste0("Image:", names(image.ext)[tif], " and band ", year, " use mean"))
      names(image.ext[[tif]][[year]]) <- names(r)[year_bands]
    }
  }
} 


rm(list = ls()[-which(ls() %in% c("image.ext"))])
save.image("image_extraction.RData")