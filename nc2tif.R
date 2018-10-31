library(RNetCDF)
library(raster)
setwd("Z:/IMAGE/")

nc.file = file
#Function to read a nc-file, extract the year 2015, 2050, 2100 and convert it to 
nc2tif <- function(nc.file){
  #open nc file and read lon, lat and time variables
  nc <- ncdf4::nc_open(nc.file)
  lat <- ncdf4::ncvar_get(nc, "latitude")
  lon <- ncdf4::ncvar_get(nc, "longitude")
  time <- ncdf4::ncvar_get(nc, "time")
  time_2015 <- which.min(abs(35 - time/365.25))
  time_2050 <- which.min(abs(80 - time/365.25))
  time_2100 <- length(time)
  
  #create vector for each year
  for(time in ls(pattern = "time_")){
    year <- as.numeric(strsplit(time, "_")[[1]][2])
    start_ini <- c(1,1,get(time))
    count_ini <- c(length(lon),length(lat),1)
    
    #create multiple vectors if variable has more than 1 value, check if variable is in first or last position
    
    if(length(nc$var[[1]]$varsize) == 4){
      varlength <- which.min(nc$var[[1]]$varsize)
      for(i in seq(nc$var[[1]]$varsize[varlength])){
        
        if(varlength == 1){
          start = c(i, start_ini); count = c(1, count_ini)
          }else{
          start = c(start_ini, i); count = c(count_ini, 1)
          }
        assign(paste0("map_", year, "_cat_", i), as.vector(ncdf4::ncvar_get(nc, start = start, count = count)))
      }
    }else{
      assign(paste0("map_", year), as.vector(ncdf4::ncvar_get(nc, start = start_ini, count = count_ini)))
    }
  }
  
  #combine vectors and lon, lat to a xyz files
  vars <- do.call(cbind, mget(ls(pattern = "map_"))) 
  names(vars) <- ls(pattern = "map_")
  
  xyz <- data.frame(
    x = rep(lon, length(lat)),
    y = rep(lat, each = length(lon)),
    vars)
  
  #create a raster and assign epsg 4326
  r <- raster::rasterFromXYZ(xyz)

  if(ncol(vars) > 3){
    class_names <- NULL
    for(i in c(2015, 2050, 2100))  class_names <- append(class_names, paste0(nc$var[[1]]$dim[[varlength]]$vals, i))
    names(r) <- class_names
  }else{
    names(r) <- paste0("y_", c(2015, 2050, 2100))
  }
  
  raster::crs(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
  raster::writeRaster(r, paste0(tools::file_path_sans_ext(basename(nc.file)), ".grd"), format = "raster", overwrite = T)
  print(paste0(paste0(tools::file_path_sans_ext(basename(nc.file)), ".grd"), " written!"))
}

#apply the function on all nc.files
nc.files <- list.files(pattern = ".nc$")

for(file in nc.files) nc2tif(file)

for(tif in list.files(pattern = ".tif$")) print(raster::stack(tif))