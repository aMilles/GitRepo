library(googledrive)
library(rgdal)
library(zoo)
library(stringi)
library(pscl)
library(MASS)
library(ggplot2)
library(reshape2)

seg = readOGR("Z:/GEC/segments_GEE_scenarios.shp")

preds <- drive_ls("Climate_Change")
pred.choice <- vector(length = NROW(unique(preds$name)))

for(i in seq(length(pred.choice))){
  pred  <- data.frame(unique(preds[,1]))[i,]
  subset <- which(preds[,1] == pred)
  times <- unlist(lapply(preds$drive_resource, function(x) x$createdTime))[subset]
  times <- as.POSIXct(times, format = "%FT%T")
  pred.choice[i] <- subset[which.max(times)]
}

#download and read predictors
for(i in pred.choice) drive_download(path = paste0("Z:/climate_change/", preds[i,1]), file = preds[i,], overwrite = T)
for(pred in list.files("Z:/climate_change")) assign(gsub(".csv", "", pred), read.csv(paste0("Z:/climate_change/", pred)))

dry_season_2015_withgeom <- dry_season_2015_withgeom[match(as.character(seg$ID), as.character(dry_season_2015_withgeom$ID)), ]
dry_season_2050_withgeom <- dry_season_2050_withgeom[match(as.character(seg$ID), as.character(dry_season_2050_withgeom$ID)), ]

preds <- gsub(".csv", "", list.files("Z:/climate_change"))
ds = preds[1]
for(ds in preds){
  ds_df <- get(ds)[,-c(1,2)]
  names(ds_df)[1] <- "median_0"
  apply(ds_df[1:100,]*86400, 1, sum)
}

for(ds in preds){
  ds_df <- get(ds)
  names(ds_df)[3] <- "median_0"
  ds_df <- ds_df[,names(ds_df)[3:72][order(as.numeric(gsub("median_", "", names(ds_df)[3:72])))]]
  assign(paste0("sum_", ds), t(apply(ds_df, 1, function(x) rollapply(x, sum, align= "right", width = 17))))
}
#sum_dry_season_2015_withgeom = cbind(dry_season_2015_withgeom$ID, t(sum_dry_season_2015_withgeom))
#sum_dry_season_2050_withgeom = cbind(dry_season_2050_withgeom$ID, t(sum_dry_season_2050_withgeom))

for(ds in paste0("sum_", preds)){
  ds_df <- get(ds)
  assign(paste0("min_", ds), apply(ds_df, 1, which.min))
}

for(ds in paste0("sum_", preds)){
  ds_df <- get(ds)
  assign(paste0("min_", ds), apply(ds_df, 1, which.min))
}
 

 
#gg_2015 <- melt(sum_dry_season_2015_withgeom)
#names(gg) <- c("segment", "week", "precipitation")
# gg_cc <- data.frame("scn_2015" = min_sum_dry_season_2015_withgeom, 
#                     "scn_2050" = min_sum_dry_season_2050_withgeom)

# samples <- sample(nrow(sum_dry_season_2015_withgeom), n = 1000)
# samples <- sample(1:6000)
# 
# sample_sum_dry_season_2015_withgeom <- melt(sum_dry_season_2015_withgeom[samples,])
# sample_sum_dry_season_2050_withgeom <- melt(sum_dry_season_2050_withgeom[samples,])
# 
# ggplot(sample_sum_dry_season_2015_withgeom, aes(x = Var2, y = value*3600, group = Var1, col = Var1))+
#   geom_line()
# 
# sum_dry_season_2015_withgeom
# hist(sample_sum_dry_season_2050_withgeom$Var1)
# 
# ggplot(gg_cc, aes(x = scn_2015 - scn_2050))+
#          geom_histogram(bins = 100)


mins <- data.frame(
  y_2015 = apply(sum_dry_season_2015_withgeom[,-1], 1, min),
  y_2050 = apply(sum_dry_season_2050_withgeom[,-1], 1, min))

write.csv(mins, "Z:/Scenario_Synthesis/SC.csv")

times_2015 <- as.POSIXct(floor(zoo::rollapply(as.numeric(seq.Date(as.Date("2014-09-01"), as.Date("2016-01-06"), by = 7)), FUN = mean, width = 2))*86400 + 43200, origin = "1970-01-01", tz = "UTC")[18:70]

times_2050 <- as.POSIXct(floor(zoo::rollapply(as.numeric(seq.Date(as.Date("2049-09-01"), as.Date("2051-01-06"), by = 7)), FUN = mean, width = 2))*86400 + 43200, origin = "1970-01-01", tz = "UTC")[18:70]


min.times <- data.frame(
  t_2015 = apply(sum_dry_season_2015_withgeom[,-1], 1, function(x) as.numeric(times_2015[which.min(x)])/1000000000),
  t_2050 = apply(sum_dry_season_2050_withgeom[,-1], 1, function(x) as.numeric(times_2050[which.min(x)])/1000000000))


seg$time_2015 <- min.times$t_2015
seg$time_2050 <- min.times$t_2050

writeOGR(seg, "Z:/GEC/segments_alltimes.shp", "segments_alltimes.shp", "ESRI Shapefile")
seg$time_2015
hist(((seg$time - seg$time_2015)*1000000000/86400) %% 365)

