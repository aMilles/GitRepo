library(maptools)
library(rgdal)
library(sp)
library(raster)
library(lubridate)
library(plyr)
library(XML)
library(geosphere)
library(rgeos)
library(matlib)
library(sf)
library(spatialEco)
library(fasterize)
library(velox)
library(parallel)
library(doParallel)
library(foreach)
library(psych)
setwd("Z:/GEC")

#read previously created GEC data
GEC <- readOGR(dsn = "GEC_points.shp")
GEC@data[,"pht_cr_"] <- as.numeric(gsub("[^0-9]", "", as.data.frame(GEC@data[,"pht_cr_"])[,1]))
summary(GEC$pht_cr_)

#get all the shape and xlsx. files of the GEC

shp.files <- list.files(pattern = "_transects.shp$", recursive = T)
shp_names <- unlist(lapply(strsplit(basename(shp.files), "_transect"), function(x) x[1]))
xlsx.files <- list.files(pattern = "elephants.xlsx$", recursive = T)
xlsx.files <- xlsx.files[which(seq(length(xlsx.files)) %in% grep("normalized", xlsx.files) &
      !seq(length(xlsx.files)) %in% grep("send_to", xlsx.files))]
GEC$trnr <- NA
GEC$flight_degree <- NA

#some summary statistics
GEC_split <- as.data.frame(GEC)
GEC_split <- base::split(GEC_split, GEC_split$srvy_cd)
unlist(lapply(GEC_split, NROW))
unlist(lapply(GEC_split, function(x) sum(x[which(x$obsrvt_ %in% c("bh", "mh", "ele_unkown")),"pht_cr_"], na.rm = T)))

##### COMBINE ALL TRANSECTS TO A SINGLE OBJECT ##### 
out <- vector("list", length = length(shp.files))

for(i in shp.files) assign(paste0("shape_", basename(i)), readOGR(i))
std_crs<- crs(shape_AGO_Lui_transects.shp)
#STEP 1: UNIFY THE SHAPEFILES TO ONE FORMAT WITH ID AND SURVEY_CODE GIVEN, TRANSFORM CRS IF NECCESSARY TO LONLAT
for(i in ls(pattern = "^shape_")){
  shp <- get(i)
  shp$ID <- seq(NROW(shp))
  shp$survey_code <- paste(strsplit(i, "_")[[1]][2:3], collapse = "_")
  shp <- shp[,which(names(shp) %in% c("ID", "survey_code"))]

  #if crs is not lonlat, transform!
  if(length(grep("longlat", as.character(crs(get(i))))) == 0){ assign("shp", spTransform(shp, std_crs))}
  assign(paste0("unified_", i), shp)
} 

#combine all transformed shapefiles to a list and stack via rbind
transects <- mget(ls(pattern = "unified_shape_"))
transects <- do.call(rbind, transects)


#### CALCULATE FLIGHT ANGLE FOR EVERY TRANSECT ####

## STEP 1: CALCULATE THE CHANGE BETWEEN TWO POINTS OF A TRANSECT
#extract the first two coordinates of transects (assuming all transects are a straight line)
cts <- lapply(coordinates(transects), function(x) x[[1]][1:2, 1:2])
#calculate the change between the two points
diffs <- lapply(cts, function(x) apply(x, 2, diff))
#stack those changes in lon and lat via rbind
a <- do.call(rbind, diffs)


### STEP 2: CONVERT TO UNIT VECTOR
#function to convert coordinats to unit vector
scalar1 <- function(x) {
  uv <- x / sqrt(sum(x^2))
  return(uv)
}

uvs <- t(apply(a, 1, scalar1))
#STEP 3: CALCULATE DEGREE
degree <- round(asin(uvs[,1])*180/pi, 0)
degree <- ifelse(uvs[,2] < 0, 90  + 90 - degree, 
                   ifelse(uvs[,1] < 0 & uvs[,2] > 0, 90 + degree + 270, degree)) 
degree <- ifelse(degree == 360, 0, degree)
plot(uvs)
text(uvs, labels = degree)
degree[which(degree < 0)]

#calculate degree of each transect, with 0 being northwards, clockwise increase
transects$DEGREE <- abs(degree)
transects$dx <- uvs[,1]
transects$dy <- uvs[,2]

#### GET TRANSECT CLOSEST TO AN OBSERVATION OF THE GEC ####
# THESE STEPS ARE MOSTLY OVERWRITTEN BY ASSIGNING OBSERVATIONS TO TRANSECTS LATER ONWARDS #
# TO AVOID ANY POTENTIAL ERRORS AFTERWARDS THESE LINES SHOULD STILL BE RUN#

# ------ start of obsolete code ------ #

  #FIRST STEP: TRANSFORM LONLAT TO PLANAR COORDINATES 
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  GEC_XYZ <- data.frame(t(apply(coordinates(GEC), 1, function(x){
    X = 6371000 * cos(deg2rad(x[2])) * cos(deg2rad(x[1]))
    Y = 6371000 * cos(deg2rad(x[2])) * sin(deg2rad(x[1]))
    Z = 6371000 * sin(deg2rad(x[2]))
    return(c(X,Y,Z))
  })))
  
  ts_coords <- lapply(coordinates(transects), function(x) x[[1]])
  ts_coords <- lapply(ts_coords, function(x) cbind(x = approx(x[,1], n = 200)$y, y = approx(x[,2], n = 200)$y))
  #ts_coords <- data.frame(do.call(rbind, lapply(coordinates(transects), function(x) x[[1]])))
  ts_length <- rep(NROW(ts_coords[[1]]), length(ts_coords))
  ts_coords <- data.frame(do.call(rbind, ts_coords))
  ts_coords$ID <- rep(seq(length(transects)), each = NROW(ts_coords[[1]])/length(transects))
  
  TS_XYZ <- data.frame(t(apply(ts_coords, 1, function(x){
    X = 6371000 * cos(deg2rad(x[2])) * cos(deg2rad(x[1]))
    Y = 6371000 * cos(deg2rad(x[2])) * sin(deg2rad(x[1]))
    Z = 6371000 * sin(deg2rad(x[2]))
    return(c(X,Y,Z))
  })))
  
  
  TS_XYZ$ID <- rep(seq(NROW(transects)), each =NROW(ts_coords[[1]])/length(transects))
  TS_XYZ$srvy <- rep(transects$survey_code, ts_length)
  GEC_XYZ$srvy <- GEC$srvy_cd
  
  ### STEP 2: ITERATE TROUGH EACH SURVEY REGION WITH A TRANSECT.shp AND FIND THE TRANSECT CLOSEST TO EACH SINGLE OBSERVATION
  ID <- rep(NA, length(GEC))
  for(srvy in unique(transects$survey_code)){
    print(srvy)
    for(i in which(GEC$srvy_cd == srvy)){
     ID[i] <- TS_XYZ$ID[TS_XYZ$srvy == srvy][which.min(sqrt((GEC_XYZ[i,1] - TS_XYZ[TS_XYZ$srvy == srvy,1])^2 + (GEC_XYZ[i,2] - TS_XYZ[TS_XYZ$srvy == srvy,2])^2))]
     if(i%%100 == 0) print(paste(i, "and", ID[i]))
    }
  }
  
  ### STEP 3: ADD A SUPER ID FOR EACH GIVEN TRANSECT TO THE GEC DATASET, THAT IS UNIQUE IN THE TRANSECT SHAPEFILE
  GEC$trnr <- ID
  GEC$flight_degree <- transects$DEGREE[GEC$trnr]
  GEC$dx <- transects$dx[GEC$trnr]
  GEC$dy <- transects$dy[GEC$trnr]
  
  transects$superID <- seq(NROW(transects))  
  
  rm(list = ls(pattern = "shape_"))

# ------ end of obsolete code ------ #


#### CREATE SEGMENTED POLYGONS PARALLEL TO TRANSECT LINES ####
## THESE STEPS ARE VALID ### 
#convert transects from geographic to planar
transects_102024 <- spTransform(transects, "+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")

#extract coordinates of transect lines
se_coords <- vector("list", length = length(transects))
for(i in seq(length(transects))){
  coords <- coordinates(transects_102024[i,])[[1]]
  se_coords[[i]] <- coords[[1]][c(1,NROW(coords[[1]])),]
}


for(i in 1:length(se_coords)) se_coords[[i]] <- cbind(se_coords[[i]], rep(i, NROW(se_coords[[i]])))


#### CREATE EQUALLY SPACED POLYGONS PARALLEL TO TRANSECT LINE ####
create_polys <- function(x){
  X_S <- x[1,1]
  X_E <- x[2,1]
  Y_S <- x[1,2]
  Y_E <- x[2,2]
  
  #0.75 aus Mittelwert der Ruanda-Daten geschätzt!! 
  # w_strp = Flight Area * 1e+06 / Strip Length / 2
  w_strp <- 0.75 * 1e+06 / 2500 / 2 

 #Get the distance rquired to cover the whole length of the transect by 2500 segments
 ext = ceiling(dist(x)/2500) * 2500
 
 #calculate the changes in x and y direction for two points on the transects with a distance of 1m
 dx = ((X_E-X_S) / dist(x))[1]
 dy = ((Y_E-Y_S) / dist(x))[1] 

 #estimate the extended end point where distance to the end point is equal to ext and therefore %% 2500 == 0,
 #optim minimizes the absolute difference between the extended point and ext
 d <- optim(par = c(ext), 
                 fn = function(pars, dx_ = dx, dy_ = dy, x1 = X_S, y1 = Y_S, ext_ = ext){
                   return(abs(sqrt( (x1 - (x1 + dx_ * pars))^2 + (y1 - (y1 + dy_ * pars))^2 ) - ext_))
                   }
   , method = "Brent", lower = 0, upper = ext + 1)$par
 
 #new end points
 X_E <- X_S + d * dx
 Y_E <- Y_S + d * dy
 
 X <- approx(c(X_S,X_E), n = 1 + ext/2500)$y
 Y <- approx(c(Y_S,Y_E), n = 1 + ext/2500)$y
 i = 2
 ts_pol <- vector("list", length = length(X)-1)
 
 for(i in seq(length(X))[-1]){
   x_s <- X[i-1]
   x_e <- X[i]
   y_s <- Y[i-1]
   y_e <- Y[i]
   
   alpha <- atan((y_e - y_s) / (x_e - x_s))
   beta_f <- cos(pi/2 - alpha)
   alpha_f <- cos(alpha)
   
   x_1 <- c(x_s + 100 * beta_f, x_s - 100 * beta_f)
   y_1 <- c(y_s - 100 * alpha_f, y_s + 100 * alpha_f)
   
   x_2 <- c(x_e + 100 * beta_f, x_e - 100 * beta_f)
   y_2 <- c(y_e - 100 * alpha_f, y_e + 100 * alpha_f)
   
   x_3 <- c(x_e + (100 + w_strp) * beta_f, x_e - (100 + w_strp) * beta_f)
   y_3 <- c(y_e - (100 + w_strp) * alpha_f, y_e + (100 + w_strp) * alpha_f)
   
   x_4 <- c(x_s + (100 + w_strp) * beta_f, x_s - (100 + w_strp) * beta_f)
   y_4 <- c(y_s - (100 + w_strp) * alpha_f, y_s + (100 + w_strp) * alpha_f)
   
   pol1 <- Polygon(matrix(c(x_1[1], x_2[1], x_3[1], x_4[1], y_1[1], y_2[1], y_3[1], y_4[1]), nc = 2))
   pol2 <- Polygon(matrix(c(x_1[2], x_2[2], x_3[2], x_4[2], y_1[2], y_2[2], y_3[2], y_4[2]), nc = 2))
 x  
   pol <- Polygons(list(pol1, pol2), ID = paste0("ts", x[1,3], "p", i-1))
   ts_pol[[i-1]] <- SpatialPolygonsDataFrame(SpatialPolygons(list(pol), proj4string=CRS("+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")), data = data.frame(ID = paste0("ts", x[1,3], "p", i-1)), match.ID = FALSE)
   
 }
 
 return(do.call(rbind, ts_pol))}
 
segments <- lapply(se_coords, create_polys)

transects_4326 <- spTransform(transects_102024, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
segments_4326 <- lapply(segments, function(x) spTransform(x, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
for(i in seq(length(segments_4326))) segments_4326[[i]]$SC <- transects[i,]$survey_code
for(i in seq(length(segments_4326))) segments_4326[[i]]$CC <- strsplit(segments_4326[[i]]$SC, "_")[[1]][1]
segments_4326 <- do.call(rbind, segments_4326)

rgdal::writeOGR(segments_4326, "Z:/GEC/segments.shp", "segments.shp", "ESRI Shapefile")