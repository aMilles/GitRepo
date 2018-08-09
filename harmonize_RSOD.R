library(maptools)
library(rgdal)
library(sp)
library(raster)
library(lubridate)
library(plyr)
library(XML)

setwd("Z:/GEC")

#use all standardized GEC xlsx.files
files <- list.files(pattern = ".xlsx", recursive = T)
files <- files[grep(files, pattern = "normalized")]
files <- files[-grep(files, pattern = "send_to")]
#ZAF_Kru is empty
files <- files[-grep(files, pattern = "ZAF_Kru")]

#read excel files
datasets <- vector("list", length(files)) 
for(i in seq(length(files))) datasets[[i]] <- as.data.frame(readxl::read_xlsx(files[i], skip = 1, sheet = "RSO"))

#############################
#### WEIRD CASE HANDLING ####
#############################

#### IF POSSIBLE, HANDLE CASES WHERE DATASETS DO NOT AGREE WITH STANDARD FORMAT OR SOMETHING IS JUST WRONG!

#skip 3 rows if column names are double but keep the 2nd row as column names and fill it with names in the third row if not given otherwise, this is done by checking wether the second and third element in one column are equal in any case
for(i in which(unlist(lapply(datasets, function(x) any(casefold(as.character(as.data.frame(x[1,]))) == names(x), na.rm = T))))){
  df.names <- names(datasets[[i]])
  df.names[grep("X__", df.names)] <- datasets[[i]][1,][grep("X__", df.names)]
  df <- readxl::read_xlsx(files[i], skip = 3, sheet = "RSO", col_names = F)
  names(df) <- unlist(df.names)
  datasets[[i]] <- df
  }

#check if 2 lat and lon cols exist and replace them so the number of na is reduced, case in COD_VIR

for(i in seq(length(files))){
  df <- datasets[[i]]
  lon.cols <- grep("lon", casefold(names(datasets[[i]])))
  if(length(lon.cols) == 2) df$lon[is.na(df$lon)] <- df[is.na(df$lon), lon.cols[2]]
  
  lat.cols <- grep("lat", casefold(names(datasets[[i]])))
  if(length(lon.cols) == 2) df$lat[is.na(df$lat)] <- df[is.na(df$lat), lat.cols[2]]
  df -> datasets[[i]]
  }

ds.names <- do.call(rbind, strsplit(basename(files), "_eleph"))[,1]

#Fill missing survey code names

for(i in seq(length(files))){
  df <- datasets[[i]]
    
  df$survey_code <- ds.names[i]
  
  df -> datasets[[i]]
}


#check wether there are NAS in count
for(i in seq(length(files))){
  df <- datasets[[i]]
  
  if(sum(is.na(df$photo_corrected_count) & is.na(df$observed_count)) > 0){
    print(sum(is.na(df$photo_corrected_count) & is.na(df$observed_count)))
    if("number_out" %in% names(df)) df$photo_corrected_count[is.na(df$photo_corrected_count)] <- df$number_out[is.na(df$photo_corrected_count)]
    print(sum(is.na(df$photo_corrected_count) & is.na(df$observed_count)))
  } 
  
  df -> datasets[[i]]
}

#####################################
#### HARMONIZE OBSERVATION CODES ####
#####################################

#get all the original observation codes of all datasets
oc <- lapply(datasets, function(x) x$observation_code)
oc <- data.frame(original = unique(do.call("c", oc)), name = NA)
oc$name <- casefold(oc$original)

#replace individual observation codes with consistent codes

#carcass categories
oc$name[grep(1, oc$name)] <- "car1"
oc$name[grep(2, oc$name)] <- "car2"
oc$name[grep(3, oc$name)] <- "car3"
oc$name[grep(4, oc$name)] <- "car4"

#from KEN_Mar, F= Fresh ...
oc$name[grep("^f$", oc$name)] <- "car1"
oc$name[grep("^r$", oc$name)] <- "car2"
oc$name[grep("^o$", oc$name)] <- "car3"
oc$name[grep("^vo$", oc$name)] <- "car4"

oc$name[grep("fresh", oc$name)] <- "car1"
oc$name[grep("recent", oc$name)] <- "car2"
oc$name[grep("old", oc$name)] <- "car3"
oc$name[grep("very old", oc$name)] <- "car4"
oc$name[grep("carcass", oc$name)] <- "car_unknown"
oc$name[grep("bones", oc$name)] <- "car_unknown"


#mixed herds (mh)
oc$name[grep("breeding herd", oc$name)] <- "mh"
oc$name[grep("bh", oc$name)] <- "mh"  #breeding herd in KEN_TSV
oc$name[grep("family", oc$name)] <- "mh"
oc$name[grep("cow", oc$name)] <- "mh"
oc$name[grep("female", oc$name)] <- "mh"
oc$name[grep("elef$", oc$name)] <- "mh"
oc$name[grep("ele-f$", oc$name)] <- "mh"


#unknown elephant type
oc$name[grep("^e$", oc$name)] <- "ele_unknown"
oc$name[grep("^el$", oc$name)] <- "ele_unknown"
oc$name[grep("^ele$", oc$name)] <- "ele_unknown"
oc$name[grep("^elephant$", oc$name)] <- "ele_unknown"

#bull herd (bh)
oc$name[grep("^b$", oc$name)] <- "bh"
oc$name[grep("bull", oc$name)] <- "bh"
oc$name[grep("male", oc$name)] <- "bh"
oc$name[grep("^elem$", oc$name)] <- "bh"
oc$name[grep("ele-m$", oc$name)] <- "bh"

#tracks
oc$name[grep("sand", oc$name)] <- "track"
oc$name[grep("track", oc$name)] <- "track"
oc$name[grep("corridor", oc$name)] <- "track"
oc$name[grep("dung", oc$name)] <- "track"
oc$name[grep("habitat", oc$name)] <- "track"
oc$name[grep("eletra$", oc$name)] <- "track"
oc$name[grep("eletrf$", oc$name)] <- "track"
oc$name[grep("trace", oc$name)] <- "track"


#misc, will rather not be considered
oc$name[grep("buf", oc$name)] <- "misc"
oc$name[grep("buffalo", oc$name)] <- "misc"
oc$name[grep("tusk", oc$name)] <- "misc"

#remaining unknowns, atm there are none!
oc$name[!oc$name %in% c("bh", "mh", "track", "misc", "car1", "car2", "car3", "car4", "car_unknown", "ele_unknown")] <- "ele_uncertain"

#replace the original observation codes in each dataset with the new classes
for(i in seq(length(datasets))){
  if(length(datasets[[i]]$observation_code) > 0) datasets[[i]]$observation_code <- oc$name[match(datasets[[i]]$observation_code, oc$original)]
}

#############################
#### HARMONIZE POSITIONS ####
#############################

#not.empty, with.pos,
not.empty <- which(unlist(lapply(datasets, NROW)) > 0 &
                   unlist(lapply(datasets, function(x) all(!is.na(x$lon)))))

ds.names <- unlist(lapply(datasets, function(x) x$survey_code[1]))

#5 files with trouble, COD_GAR hat utm only, COD_VIR has 2 missing lonlats, ETH_NW has one missing lonlat, KEN_TSV has ~10 missing longlat in 680 obs, ZWE_ZV has 2 missing observations in 700 obs. 

#if datasets has utm, such as COD_GAR convert to lonlat wgs84, remove points without any spatial information

for(i in seq(length(datasets))[-not.empty]){
  df <- datasets[[i]]
  
  #remove points with no spatial information
  if(length(which(is.na(df$lon) & is.na(df$utm_x)) > 0)) df <- df[-which(is.na(df$lon) & is.na(df$utm_x)),]
  
  #which points have utm information, THIS IS CURRENTLY NOT GENERALIZED FOR ALL CASES, UTM WILL CURRENTLY BE TRANSFORMED FROM ONE UTM ZONE TO WGS84 - IF MORE UTM COORDINATES APPEAR, THIS HAS TO BE ADJUSTED!!!!!
  missing.lonlat <- which(is.na(df$lon) & !is.na(df$utm_x))
  
  if(length(missing.lonlat) > 0){
    df <- df[missing.lonlat, ]
  
    xy <- data.frame(ID = 1:nrow(df), X = df$utm_x, Y = df$utm_y)
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=utm +zone=35 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))

    df[, c("lon", "lat")] <- coordinates(res) 
  }
  
  long2UTM <- function(long) {
    (floor((long + 180)/6) %% 60) + 1
  }
  
  xy <- data.frame(ID = 1:nrow(df), X = df$lon, Y = df$lat)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
  res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  
  df[, c("utm_x", "utm_y")] <- coordinates(res) 
  
  datasets[[i]] <- df
    
  }

#no more missing positions!
lapply(datasets, function(x) sum(is.na(x$lon)))


##########################
#### HARMONIZE TIMES  ####
##########################


backup -> datasets

#not empty and contains lonlat, should be everything except for kurgernationalpark

ds.names <- do.call(rbind, strsplit(basename(files), "_eleph"))[,1]

not.empty <- which(unlist(lapply(datasets, NROW)) > 0 &
                     unlist(lapply(datasets, function(x) sum(is.na(x$lon)) == 0)))    

#not.empty[-which(ds.names[not.empty] == "KEN_Mar")]
df <- datasets[[10]]
ds.names

#harmonize local date

#many special cases are handled here

for(dataset in not.empty){

  df <- datasets[[dataset]]
  
  #some datasets contain local times and utc times with am and pm - those are converted to numerics by string split to seconds of the day. pm is added with 12h or 43200seconds
  if(is.character(df$utc_date_time)){
    if(length(grep("am", casefold(df$utc_date_time))) > 0){
      df$utc_date_time <- NA
      am <- grep("am", casefold(df$utc_date_time))
      pm <- grep("pm", casefold(df$utc_date_time))
      df$local_time[am] <- as.numeric(unlist(strsplit(casefold(df$local_time[am]), " am")))
      df$local_time[pm] <- as.numeric(unlist(strsplit(casefold(df$local_time[pm]), " pm"))) + 43200
    }
  }  

    
  #delete RSO where no information about date and time is given, this is rarely the case
  no.date <- which(is.na(df$local_date) & is.na(df$utc_date_time))
  if(length(no.date) > 0) df <- df[-no.date, ]
  
  #no time but date information, you have to estimate time in this cases
  estm.time <- which(is.na(df$local_time) &
                       !is.na(df$local_date))
  
  #delete where time zone is a posixct, time zone should be a character, therefore those columns are resetted
  if(is.POSIXct(df$tz)){
    df$tz <- NA
    df$tz <- as.character(df$tz)
  } 
  
  #convert time to seconds of the day
  
  #where time is available convert it to numeric as seconds of the day...
  
    #work-around for negative indexing (if not times need to be estimated)
  if(length(estm.time) == 0) estm.time <- nrow(df) + 1 
  
    #convert time to seconds of the day, this is done dependant on the current format
  if(is.POSIXct(df$local_time)){
    df$local_time <- as.numeric(df$local_time) %% 86400
  }
  else{
    if(!is.numeric(df$local_time)){
      df$local_time <- as.numeric(hms(df$local_time))
    } 
  }
  
  
  # ... if not assume its 10:30 AM (flights were mostly in the morning)
  if(estm.time[1] != (nrow(df) + 1)) df$local_time[estm.time] <- 43200 - 3600 - 1800
  
  
  #no tz given but utc, retrieve from coordinates if available - do this for 100 positions..
  if(length(which(is.na(df$tz) & !is.na(df$lon))) > 0){
    
    for(i in sample(which(is.na(df$tz) & !is.na(df$lon) & !is.na(df$local_date)), size = 100, replace = T)){
      apiurl <- sprintf("https://maps.googleapis.com/maps/api/timezone/%s?location=%s,%s&timestamp=%d&sensor=%s",
                        "xml", 
                        df$lon[i], 
                        df$lat[i], 
                        as.POSIXct(df$local_date[i]), 
                        "false")
      df$tz[i] <- xmlParse(readLines(apiurl))[["string(//time_zone_id)"]]
      
      if(i %% 100 == 0) print(i)
    }
  }  
  
  
  #convert local time and local date to utc datetime, add local time to local date, respect timezone. Multiple time zones not considered!!
  local_tz <- df$tz[which(!is.na(df$tz) & df$tz != "")][1]
  

  ## if no time information is given approx if possible... NOT IMPLEMENTED AT THE MOMENT!
  # approx.try <- na.approx(df$local_time, maxgap = 2)
  # if(length(approx.try) == length(df$local_time))  df$local_time <- approx.try 
  ## ... or fill with 12:00 AM in seconds
  # 
  # df$local_time[is.na(df$local_time)] <- 43200
  
  #if local times and timezone given but no utc time, convert to numeric time format with local time...
  no.utc.time <- which(is.na(df$utc_date_time) & !is.na(df$local_date))
  
  if(length(no.utc.time) == nrow(df)){
    
  df$utc_date_time[no.utc.time] <- as.POSIXct(
    x = as.numeric(df$local_time[no.utc.time]) + 
      as.numeric(as.POSIXct(df$local_date[no.utc.time], 
                            tz = local_tz)),
    origin = as.POSIXct("1970-01-01 00:00:00", tz = local_tz),
    tz = local_tz)
  
  
  #... approximate missing points if difftime is smaller than 3 hours... NOT IMPLEMENTED
  #df$utc_date_time <- na.approx(df$utc_date_time, maxgap = 1)

  #... and convert it to utc time format
  df$utc_date_time <- as.POSIXct(
    x = format.POSIXct(as.POSIXct(df$utc_date_time, tz = local_tz, origin = "1970-01-01 00:00:00"), 
                       tz = "utc", 
                       origin = "1970-01-01 00:00:00"), tz = "utc")
  }
  
  
  datasets[[dataset]] <- df
  not.empty <- not.empty[-match(dataset, not.empty)]
  
}

#now there should be no more missing utc times!
lapply(datasets, function(x) sum(is.na(x$utc_date_time)))

###### CREATE A COMPLETE CSV FILE ####


# get the column names that are inside all rsods
min.cols <- Reduce(intersect, lapply(datasets, names))
red.df <- lapply(datasets, function(x, cols) x[,cols], cols = min.cols)

red.df <- do.call(rbind, lapply(red.df, function(x) x[,-c(5,6,7)]))

#apparently in ZAF_TUL lat and lon was confused
red.df[red.df$survey_code == "ZAF_TUL", c("lat", "lon")] <- red.df[red.df$survey_code == "ZAF_TUL", c("lon", "lat")]

write.csv(red.df, "agg_rso.csv")


###### CREATE A SHAPEFILE #####

#lonlat extraction
xy <- red.df[, c("lon","lat")]

spdf <- SpatialPointsDataFrame(coords = xy, data = red.df,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
red.df$photo_corrected_count[is.na(red.df$photo_corrected_count)] <- red.df$observed_count[is.na(red.df$photo_corrected_count)] 
summary(as.numeric(red.df$photo_corrected_count))
writeOGR(spdf, "GEC_points.shp", "GEC_points.shp", "ESRI Shapefile")

red.df_eleonly <- red.df[red.df$observation_code %in% c("mh", "bh", "ele_unkown"), ]
writeOGR(spdf, "GEC_points_eleonly.shp", "GEC_points_eleonly.shp", "ESRI Shapefile")


