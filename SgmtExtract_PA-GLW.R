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

segments_4326 <- readOGR(dsn = "GEC_points.shp")
transects_4326 <- readOGR(dsn = "transects.shp")
GEC <- readOGR(dsn = "GEC_points.shp")

what <- sample(1:4000, size = 1)
plot(segments_4326[[what]])
segments_4326 <- do.call(bind, segments_4326)
rm(ts_coords)

plot(transects_4326[what,], col = "red", add = T)
GEC[1,]
intersect(GEC[1,], transects_4326)
centroids <- centroid(segments_4326)
sample <- sample(1:4000, size = 1)


what <- which.min(geosphere::distGeo(GEC[sample,], centroids))
plot(segments_4326[what, ])
plot(GEC[sample,], add = T)
which.min(gDistance(GEC[1,], centroi))

GEC$sg_nr <- NA
GEC$d2sg <- NA

for(i in seq(NROW(GEC))){
  dists <- geosphere::distGeo(GEC[i,], centroids)
  GEC[i,"sg_nr"] <- segments_4326$ID[which.min(dists)]
  GEC[i,"d2sg"] <- dists[which.min(dists)]
  if(i %% 100 == 0) print(i)
}


#PREPARE GEC DATA

#add segment IDs, remove carcasses, tracks etc, create an ID for GEC, replace photo-corrected numbers if they are NA with observed numbers 

GEC$sg_nr <- segments_4326$ID[GEC$sg_nr]
GEC <- GEC[GEC$obsrvt_ %in% c("mh", "bh", "ele_unkown"), ]
GEC$ID <- seq(NROW(GEC))
GEC$pht_cr_ <- ifelse(is.na(GEC$pht_cr_), GEC$obsrvd_, GEC$pht_cr_) 

#check which GEC points have a nearby segments, use those data points that are reported without transects, have a block definition or are more than 1500m distant to the next segment as points. less than 10 percent

GEC_aspg <- GEC[!GEC$srvy_cd %in% c("KEN_Mar", "MLI_Gou", "MWI_KAS", "MWI_KAS", "MWI_LIW", "ZAF_Kru", "ZAF_SHA", "ZAF_TUL"),]
GEC_aspg <- GEC_aspg[is.na(GEC_aspg$blck_df), ]
GEC_aspg <- GEC_aspg[GEC_aspg$d2sg < 1500, ]

GEC_aspo <- GEC[-GEC_aspg$ID, ]

GEC_aspo <- as.data.frame(GEC_aspo)
GEC_aspg <- as.data.frame(GEC_aspg)

# segments_4326$obs_type <- "nothing"
# segments_4326$obs_count <- 0
# segments_4326$utc_time <- NA

aspg_sums <- data.frame(aggregate(GEC_aspg$pht_cr_, by = list(GEC_aspg$sg_nr), FUN = function(x) mean(x, na.rm = T)))
aspg_types <- aggregate(GEC_aspg$obsrvt_, by = list(GEC_aspg$sg_nr), FUN = function(x) x[1])
aspg_times <- aggregate(GEC_aspg$utc_dt_, by = list(GEC_aspg$sg_nr), FUN = function(x) x[1])

aspg_sums$ID <- aspg_sums$Group.1
aspg_types$ID <- aspg_types$Group.1 
aspg_times$ID <- aspg_times$Group.1

df.segments_4326 <- data.frame(segments_4326)
temp.df <- plyr::join(df.segments_4326, aspg_sums)
segments_4326$obs_count <- temp.df$x
segments_4326$obs_count[is.na(segments_4326$obs_count)] <- 0


df.segments_4326 <- data.frame(segments_4326)
temp.df <- plyr::join(df.segments_4326, aspg_types)
segments_4326$obs_types <- temp.df$x

df.segments_4326 <- data.frame(segments_4326)
temp.df <- plyr::join(df.segments_4326, aspg_times)
segments_4326$utc_time <- temp.df$x


### Exract Information about wether the transect (partly) covers a protected area or not
# use fasterize to speed up the process
#read the portected area shapefile created in ...
PA <- rgdal::readOGR("Z:/PA/PA_filtered.shp")
#unify crs
PA <- spTransform(PA, crs(segments_4326))
#create a dummy raster, that cover whole africa for rasterization, increase resolution to 500x500m
dummy <- raster("Z:/GLW/AF_Cattle1km_AD_2010_v2_1.tif")
res(dummy) <- res(dummy) / 2
#rasterize
PA_raster <- fasterize(sf::st_as_sf(PA), dummy)
PA_raster[is.na(PA_raster)] = 0
#extract the values of the segments from PA_raster
intersects <- extract(x = PA_raster, y = segments_4326)
check <- as.numeric(unlist(lapply(intersects, function(x) !all(x == 0))))
segments_4326$PA <- check

raster::writeRaster(PA_raster, "Z:/PA/PA_raster.tif", format = "GTiff", overwrite = T)


### Extract TLU
TLU <- raster("Z:/GLW/TLU.tif")
TLU_extract <- extract(TLU, segments_4326)

# cl <- makeCluster(detectCores() - 12)
# registerDoParallel(cl)
# out <- foreach(i=1:NROW(segments_4326), .combine=rbind, .packages = c("raster")) %dopar% {
#   return(mean(unlist(extract(TLU, segments_4326[i,])), na.rm = T))
#   if(i %% 100 == 0) print(i)
# }
# stopCluster(cl)
# 

#This takes quite some time...
TLUs <- vector(length = NROW(transects_4326))
for(i in 1:NROW(segments_4326)){
  TLUs[i] <- mean(unlist(extract(TLU, segments_4326[i,])), na.rm = T)
  print(i)
}
segments_4326$TLU <- TLUs




v.TLU <- velox(TLU)
TLU_extract <- v.TLU$extract(segments_4326, fun = function(x) mean(x, na.rm = T))
summary(TLU_extract)
hist(TLU_extract)

rgdal::writeOGR(segments_4326, "Z:/GEC/segments.shp", "segments.shp", "ESRI Shapefile")


# 
# 
# 
# dataset <- as.data.frame(segments_4326)
# dataset$obs_types
# out <- glm(obs_count ~ PA + TLU, data = dataset, family = "poisson")
# summary(out)
# psych::pairs.panels(dataset[,c(2,3,6)])
# plot(log10(dataset$TLU), dataset$obs_count)
# library(ggplot2)
# ggplot(dataset, aes(x = log10(TLU), fill = PA))+
#   geom_histogram()+
#   facet_wrap(~PA)
# hist(log10(TLUs), breaks = 100, col = ifelse)
# par(mfrow = c(2,2))
# plot(out)
# ?glm

