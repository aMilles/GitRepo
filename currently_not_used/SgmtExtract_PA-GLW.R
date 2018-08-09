### THIS SHOULD CURRENTLY NOT BE USED!!!

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

segments_4326 <- readOGR(dsn = "segments.shp")
GEC <- readOGR(dsn = "GEC_points_eleonly.shp")
GEC <- GEC[GEC$obsrvt_ %in% c("bh", "mh", "ele_unknown"), ]
summary(GEC$obsrvt_)
# what <- sample(1:4000, size = 1)
# plot(segments_4326[what,])
# 
# rm(ts_coords)
# 
# plot(transects_4326[what,], col = "red", add = T)
# GEC[1,]
# intersect(GEC[1,], transects_4326)
centroids <- centroid(segments_4326)
# sample <- sample(1:4000, size = 1)
# 
# 
# what <- which.min(geosphere::distGeo(GEC[sample,], centroids))
# plot(segments_4326[what, ])
# plot(GEC[sample,], add = T)
# which.min(gDistance(GEC[1,], centroi))

GEC$sg_nr <- NA
GEC$d2sg <- NA
GEC$td2sg <- NA


for(SC in casefold(as.character(unique(segments_4326$SC)))){
  seg_subset <- segments_4326[casefold(as.character(segments_4326$SC)) == SC,]
  SC_cents <- centroid(seg_subset)
  GEC_points <- GEC[casefold(as.character(GEC$srvy_cd)) == SC, ]
  print(head(GEC_points))
  for(p in seq(nrow(GEC_points))){
    dists <- geosphere::distGeo(GEC_points[p, ], SC_cents)
    # 
    # if(abs(as.numeric(difftime(as.POSIXct(GEC_points$utc_dt_[p], origin = "1970-01-01"),
    #             as.POSIXct(segments_4326$gpx_time[which.min(dists)], origin = "1970-01-01"),units = "hour"))) > 6){
    #   print('date is strange')
    #   print(as.POSIXct(segments_4326$gpx_time[which.min(dists)], origin = "1970-01-01"))
    #   print(as.POSIXct(GEC_points$utc_dt_[p], origin = "1970-01-01"))
    # }
    
    GEC_points[p,"sg_nr"] <- as.character(seg_subset$ID[which.min(dists)])
    GEC_points[p,"d2sg"] <- dists[which.min(dists)]
    GEC_points[p,"td2sg"]<- as.numeric(difftime(as.POSIXct(GEC_points$utc_dt_[p], origin = "1970-01-01"),
                                        as.POSIXct(seg_subset$gpx_time[which.min(dists)], origin = "1970-01-01"),
                                        units = "hour"))
     
  }
  assign(paste0("GEC_", SC), GEC_points)
  print(SC)
}

#remove points that are mor than 14 days or 2000m apart from next gpx.points
GEC_test <- do.call(rbind, mget(paste0("GEC_", casefold(as.character(unique(segments_4326$SC))))))
GEC_test <- GEC_test[abs(GEC_test$td2sg) < 14*24,]
GEC_test <- GEC_test[GEC_test$d2sg < 2000,]

for(SC in casefold(as.character(unique(segments_4326$SC)))) plot(ecdf(GEC_test$td2sg[casefold(as.character(GEC_test$srvy_cd)) == SC]), main = SC)
for(SC in casefold(as.character(unique(segments_4326$SC)))) plot(ecdf(GEC_test$d2sg[casefold(as.character(GEC_test$srvy_cd)) == SC]), main = SC, lty = 1, pch = "+", log = "x")

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

aspg_sums <- data.frame(aggregate(GEC_aspg$pht_cr_, by = list(GEC_aspg$sg_nr), FUN = function(x) sum(x, na.rm = T)))
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


times <- as.POSIXct(segments_4326$utc_time, tz = "UTC")
segments_4326$time_backup <- times

#SC = "AGO_Lui"
segments_4326$utc_time_int <- NA
segments_4326$inttime <- as.POSIXct(segments_4326$utc_tme, tz = "UTC")
#TS = "ts1"
for(TS in unique(unlist(lapply(strsplit(as.character(segments_4326$ID), "p"), function(x) x[1])))){
  
  which.ts <- which(stri_count_regex(as.character(segments_4326$ID), paste0(TS, "p")) == 1)
  if(sum(!is.na(times[which.ts])) > 1){
    if(abs(as.numeric(difftime(max(times[which.ts], na.rm = T), min(times[which.ts], na.rm = T), units = ("hours"))) > 1)){
      which(TS)
    }
    ts.times <- na.fill(times[which.ts], "extend")
    segments_4326$utc_time_int[which.ts] <- ts.times
  }
  if(sum(!is.na(times[which.ts])) == 1){
    ts.times <- times[which.ts]
    segments_4326$utc_time_int[which.ts] <- rep(ts.times[!is.na(ts.times)], length(ts.times)) 
  }
}


for(SC in unique(segments_4326$SC)){
  test <- data.frame(segments_4326[segments_4326$SC == SC,])
  plot(x = as.POSIXct(test$utc_time_int, origin = "1970-01-01"), y = 1:nrow(test), main = SC)
  #print(sort(na.omit(as.POSIXct(test$utc_time))))
}

sum(is.na(segments_4326$utc_time_int))

for(SC in unique(segments_4326$SC)){
  which.sc <- which(segments_4326$SC == SC)
  missing.Dates <- as.Date(na.fill(times[which.sc], "extend"))[is.na(segments_4326$utc_time_int[which.sc])]
  hist(as.numeric((times[which.sc]))%%86400/3600, main = SC, freq = F, breaks = 0:24,xlim= c(0,24))
  density <- hist(as.numeric((times[which.sc]))%%86400/3600, main = SC, freq = F, breaks = 0:24,xlim= c(0,24))
  missing.Times <- sample(density$mids, size = length(missing.Dates), prob = density$density, replace = T) * 3600
  missing.Datetimes <- as.POSIXct(as.numeric(missing.Dates) * 86400 + missing.Times, origin = "1970-01-01", tz = "UTC")
  segments_4326$utc_time_int[which.sc[is.na(segments_4326$utc_time_int[which.sc])]] <- missing.Datetimes
}
max(as.POSIXct(segments_4326$utc_tm_))
#
segments_4326$inttime <- as.numeric(as.POSIXct(segments_4326$utc_tm_, origin = "1970-01-01"))
segments_4326$inttime <- segments_4326$inttime/1000000000
#writeOGR(segments_4326, "Z:/GEC/segments.shp", layer = "segments.shp", driver = "ESRI Shapefile", layer_options = "RESIZE=YES")

for(SC in unique(segments$SC)){
  test <- data.frame(segments[segments$SC == SC,])
  plot(x = as.POSIXct(test$utc_tm_), y = 1:nrow(test), main = SC)
}


rgdal::writeOGR(segments_4326, "Z:/GEC/segments.shp", "segments.shp", "ESRI Shapefile")

