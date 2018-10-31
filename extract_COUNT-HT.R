library(rgdal)
library(sp)
library(raster)
library(geosphere)
library(rgeos)
library(sf)
library(ggplot2)
setwd("Z:/GEC")
segments_4326 <- readOGR(dsn = "segments_GEE.shp")
segments_SC <- readOGR(dsn = "segments.shp")
summary(segments_4326$time)
#read GEC data.. select the observations codes you're interested in. 
GEC <- readOGR(dsn = "GEC_points_eleonly.shp")
GEC <- GEC[GEC$obsrvt_ %in% c("bh", "mh", "ele_unknown"), ]

centroids <- centroid(segments_4326)
segments_4326$ID <- as.character(segments_4326$ID)
segments_SC$ID <- as.character(segments_SC$ID)
segments_4326 <- merge(segments_4326, segments_SC)
#segments_4326$time <- segments_4326$time * 1e+9
GEC$sg_nr <- NA
GEC$d2sg <- NA
GEC$td2sg <- NA
table(GEC$srvy_cd)

# merge each observation with the spatially closest spatial subunit
for(SC in casefold(as.character(unique(segments_SC$SC)))){
  seg_subset <- segments_4326[casefold(as.character(segments_4326$SC)) == SC,]
  SC_cents <- centroid(seg_subset)
  GEC_points <- GEC[casefold(as.character(GEC$srvy_cd)) == SC, ]
  print(head(GEC_points))
  for(p in seq(nrow(GEC_points))){
    dists <- geosphere::distGeo(GEC_points[p, ], SC_cents)
    
    GEC_points[p,"sg_nr"] <- as.character(seg_subset$ID[which.min(dists)])
    GEC_points[p,"d2sg"] <- dists[which.min(dists)]
    GEC_points[p,"td2sg"]<- as.numeric(difftime(as.POSIXct(GEC_points$utc_dt_[p], origin = "1970-01-01"),
                                                as.POSIXct(seg_subset$time[which.min(dists)] * 1e+9, origin = "1970-01-01"),
                                                units = "hour"))
    
  }
  assign(paste0("GEC_", SC), GEC_points)
  print(SC)
}

#remove points that are mor than 14 days or 2000m apart from next gpx.points
GEC_test <- do.call(rbind, mget(paste0("GEC_", casefold(as.character(unique(segments_4326$SC))))))

temp_dist <- ggplot(data.frame(GEC_test), aes(x = abs(td2sg)))+
  geom_vline(xintercept = 14*24, linetype = "dashed", color = "gray", size = 1)+
  geom_freqpoly(bins = 20, fill = "gray")+
  scale_x_log10(breaks = c(0.1, 10, 1000))+
  theme_bw()+
  xlab("log10(temporal distance to spatial subunits) [h]")
?scale_x_continuous
spatial_dist <- ggplot(data.frame(GEC_test), aes(x = abs(d2sg)))+
  geom_vline(xintercept = 2000, linetype = "dashed", color = "gray", size = 1)+
  geom_freqpoly(bins = 50, fill = "gray")+
  scale_x_log10()+
  theme_bw()+
  xlab("log10(spatial distance to spatial subunits) [m]")

pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/temporal_spatial_distance.pdf")
gridExtra::grid.arrange(spatial_dist, temp_dist)
dev.off()

GEC_test <- GEC_test[abs(GEC_test$td2sg) < 14 * 24,]
GEC_test <- GEC_test[GEC_test$d2sg < 2000,]

#check how the spatiotemporal distances are distributed on a ecdf
for(SC in casefold(as.character(unique(segments_4326$SC)))) plot(ecdf(GEC_test$td2sg[casefold(as.character(GEC_test$srvy_cd)) == SC]), main = SC)
for(SC in casefold(as.character(unique(segments_4326$SC)))) plot(ecdf(GEC_test$d2sg[casefold(as.character(GEC_test$srvy_cd)) == SC]), main = SC, lty = 1, pch = "+")

GEC_test <- as.data.frame(GEC_test)

IDS <- data.frame(ID = segments_4326$ID)

GEC_test$pht_cr_ <- as.integer(as.character(GEC_test$pht_cr_))
GEC_test$pht_cr_[is.na(GEC_test$pht_cr_)] <- as.integer(as.character(GEC_test$obsrvd_[is.na(GEC_test$pht_cr_)]))
HT <- GEC_test[,c("sg_nr", "obsrvt_")]

agg_HT <- data.frame(ID = unique(HT$sg_nr), HT = NA)

both <- 0 
for(ID in unique(HT$sg_nr)){
  types <- as.character(HT$obsrvt_[HT$sg_nr == ID])
  
  if(all.equal.character("bh", types) == T){
    agg_HT[agg_HT$ID == ID, 2] <- "bh"
  }else{
    if(all(c("bh","mh") %in% types)){
     agg_HT[agg_HT$ID == ID, 2] <- "both"
    }else{
     agg_HT[agg_HT$ID == ID, 2] <- "mh"
    } 

  }     
}
HT <- agg_HT
summary(as.factor(HT$HT))
OBS <- COUNT <- GEC_test[,c("sg_nr", "pht_cr_")]
COUNT <- aggregate(COUNT$pht_cr_, by = list(COUNT$sg_nr), FUN = function(x) sum(x, na.rm = T))
OBS <- aggregate(as.numeric(OBS$pht_cr_> 0), by = list(OBS$sg_nr), FUN = function(x) sum(x, na.rm = T))
names(OBS)[1] <- names(COUNT)[1] <- "ID"
names(OBS)[2] <- "weight"
COUNT <- base::merge(IDS, COUNT, by = "ID", all.x = T)
OBS <- base::merge(IDS, OBS, by = "ID", all.x = T)
COUNT$x[is.na(COUNT$x)] <- 0

HT <- base::merge(IDS, HT, by = "ID", all.x = T)

summary(HT)
summary(OBS)
summary(COUNT)
write.csv(OBS, file = "Z:/predictors/REPS.csv")
write.csv(HT, file = "Z:/predictors/HT.csv")
write.csv(COUNT, file = "Z:/predictors/COUNT.csv")