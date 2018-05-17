library(maptools)
library(rgdal)
library(sp)
library(raster)
library(XML)
library(zoo)

setwd("Z:/GEC")
files <- list.files(pattern = ".xlsx", recursive = T)
files <- files[grep(files, pattern = "normalized")]
files <- files[-grep(files, pattern = "send_to")]
#ETH_BAB fso is empty / Strata Count
files <- files[-grep(files, pattern = "ETH_BAB")]
#ETH_NW fso is empty / Total Count (?)
files <- files[-grep(files, pattern = "ETH_NW")]
#ZAF_Kru is empty
files <- files[-grep(files, pattern = "ZAF_Kru")]
#COD_GAR is did a total count, no FSO given
files <- files[-grep(files, pattern = "COD_GAR")]
#KEN_LAM did transect count but did not assess fso, no altitude
files <- files[-grep(files, pattern = "KEN_Lam")]
#KEN_LAI fso is mostly empty, altitude inside RSO as GSPAlt in ft.
files <- files[-grep(files, pattern = "KEN_LAI")]
#XWA_TBC has no fso information
files <- files[-grep(files, pattern = "XWA_TBC")]

a <- FSOD[[8]]

#read excel files
FSOD <- vector("list", length(files)) 
for(i in seq(length(files))) FSOD[[i]] <- readxl::read_xlsx(files[i], skip = 1, sheet = "FSO")
unlist(lapply(FSOD, function(x) x$survey_code[1]))

#
list <- lapply(FSOD, function(x) as.data.frame(apply(x[sample(nrow(x), size = ifelse(nrow(x) > 5, 5, nrow(x))),], 2, as.character)))
FSOD.df <- data.frame(do.call(rbind.fill, list))
#FSOD.df$survey_code <- rep(ds.names, each = 5)


#only 10 files contain a position
unique(FSOD.df$survey_code[which(!is.na(FSOD.df$lon))]) 

#only 6 files contain position and transect number
unique(FSOD.df$survey_code[which(!is.na(FSOD.df$lon) & !is.na(FSOD.df$transect_number))]) 

#only 2 files contain position, transect number and flight direction
unique(FSOD.df$survey_code[which(!is.na(FSOD.df$lon) & !is.na(FSOD.df$transect_number) & !is.na(FSOD.df$flight_direction))]) 


#only select those FSO where longitude is given and FSO is not empty.. currently only 9 out of 24 files!
not.empty <- which(unlist(lapply(FSOD, NROW)) > 0 & unlist(lapply(FSOD, function(x) all(!is.na(x$lon)))))
not.empty
#harmonize altitude unit and create utcdatetime
for(dataset in not.empty){
  
  df <- FSOD[[dataset]]
  
  
  #delete FSO where no information about date and time is given, 3 cases in AGO Lui
  no.date <- which(is.na(df$local_date) & is.na(df$utc_date_time))
  if(length(no.date) > 0) df <- df[-no.date, ]
  
  #no time information, ~ 2000 cases
  estm.time <- which(is.na(df$local_time) &
                     is.na(df$utc_date_time) &
                     !is.na(df$local_date))
  
  
  #no tz given but utc, retrieve from coordinates if available - do this for 100 positions
  if(length(which(is.na(df$tz) & !is.na(df$lon))) > 0){
    
    for(i in sample(which(is.na(df$tz) & !is.na(df$lon) & !is.na(df$utc_date_time)), size = 100, replace = T)){
      apiurl <- sprintf("https://maps.googleapis.com/maps/api/timezone/%s?location=%s,%s&timestamp=%d&sensor=%s",
                        "xml", 
                        df$lon[i], 
                        df$lat[i], 
                        as.POSIXct(df$utc_date_time[i]), 
                        "false")
      df$tz[i] <- xmlParse(readLines(apiurl))[["string(//time_zone_id)"]]
      
      if(i %% 100 == 0) print(i)
      
    }
  }  
  
  #no tz given and no utc, retrieve from coordinates if available - do this for 100 positions
  # if(length(which(is.na(df$tz) & !is.na(df$lon)) &is.na(df$utc_date_time)) > 0){
  #   
  #   for(i in sample(which(is.na(df$tz) & !is.na(df$lon)), size = 100, replace = T)){
  #     apiurl <- sprintf("https://maps.googleapis.com/maps/api/timezone/%s?location=%s,%s&timestamp=%d&sensor=%s",
  #                       "xml", 
  #                       df$lon[i], 
  #                       df$lat[i], 
  #                       as.POSIXct(df$utc_date_time[i]), 
  #                       "false")
  #     df$tz[i] <- xmlParse(readLines(apiurl))[["string(//time_zone_id)"]]
  #     
  #     if(i %% 100 == 0) print(i)
  #     
  #   }
  # } 
  #  
  #convert local time and local date to utc datetime, add local time to local date, respect timezone. Multiple time zones not considered!!
  local_tz <- df$tz[which(!is.na(df$tz) & df$tz != "")][1]
  
  #handle Angola utc, it contains milliseconds, two time zones (i choose only one)..
  if(df$survey_code[1] == "AGO_Lui"){
    local_tz = "utc"
    df[, c("local_date", "local_time")] <- do.call(rbind, strsplit(df$utc_date_time, "T"))
    df$local_time <- do.call(rbind, strsplit(df$local_time, "\\."))[,1]
    
    df$utc_date_time <- 
      as.numeric(as.POSIXct(df$local_time, tz = local_tz, format = "%H:%M:%S"))%%86400 + 
      as.numeric(as.POSIXct(df$local_date, tz = local_tz))
  }
  
  
  #if local times and timezone given but no utc time yet convert to numeric time format with local time...
  no.utc.time <- which(is.na(df$utc_date_time) & !is.na(df$local_time) & !is.na(df$local_date) & !is.na(df$tz))

  df$utc_date_time[no.utc.time] <- as.POSIXct(
    x = (as.numeric(as.POSIXct(df$local_time[no.utc.time], tz = local_tz)))%%86400 + 
         as.numeric(as.POSIXct(df$local_date[no.utc.time], tz = local_tz)),
    origin = as.POSIXct("1970-01-01 00:00:00", tz = local_tz),
    tz = local_tz)
  
  
  #... approximate missing points ...
  df$utc_date_time <- na.approx(df$utc_date_time)
  
  #... and convert it to utc time format
  df$utc_date_time <- as.POSIXct(
    x = format.POSIXct(as.POSIXct(df$utc_date_time, tz = local_tz, origin = "1970-01-01 00:00:00"), 
    tz = "utc", 
    origin = "1970-01-01 00:00:00"), tz = "utc")
  
  
  FSOD[[dataset]] <- df
  not.empty <- not.empty[-match(dataset, not.empty)]

}