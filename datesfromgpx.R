library('plotKML')
library(rgdal)
library(stringi)
library(XML)
library(xml2)
segments <- readOGR("Z:/GEC/segments.shp")


#extract dates of the gpx files, times are set to 8am - this will be used if gpx does not contain times
ext <- regmatches(basename(gpx.files), gregexpr("\\d+\\d", basename(gpx.files)))
dates <- lapply(ext, function(x){
  x <- x[which(stringi::stri_length(x) == 8)]
  year = stringi::stri_sub(x, 1, 4)
  month = stringi::stri_sub(x, 5, 6)
  day = stringi::stri_sub(x, 7, 8)
  stringi::stri_datetime_create(year, month, day, hour = 8,tz = "UTC")
})

#read all .gpx files - if there is no date replace with dates given above
gpx.files <- list.files("Z:/GEC", pattern = ".gpx$", recursive = T)
SC <- "AGO_Lui"
for(SC in as.character(unique(segments$SC))){
  SC.gpx.files <- gpx.files[stri_count_regex(gpx.files,SC) >= 1]
  
  for(gpx.file in SC.gpx.files){
    assign(basename(gpx.file), do.call(rbind, lapply(track.readGPX.element(paste0("Z:/GEC/", gpx.file)), 
                                      function(x){
      if(length(x[[1]]$time) == 0) x[[1]]$time <- dates[[which(gpx.file == gpx.files)]]
        return(x[[1]][, c("lon", "lat", "time")])
    } )))
    
  }
  
  assign(paste0(SC, "_gpx.tot"), do.call(rbind, mget(ls(pattern = ".gpx$"))))
  rm(list = ls(pattern = ".gpx$"))
  print(SC)
}       


#THERE ARE .GPX-Files - handle these case by using the date(times) given in the .xlsx files

##GET BWA_NOR as it has no .gpx but complete FSO
BWA_NOR_gpx.tot <- readxl::read_xlsx("Z:/GEC/raw_data/BWA_NOR/normalized/BWA_NOR_elephants.xlsx", col_names = T, skip = 1, sheet = "FSO")
attr(BWA_NOR_gpx.tot$local_date, "tzone") <- "Africa/Gaborone"
attr(BWA_NOR_gpx.tot$local_time, "tzone") <- "Africa/Gaborone"
BWA_NOR_gpx.tot$time <- as.POSIXct(as.numeric(BWA_NOR_gpx.tot$local_date)+zoo::na.fill(as.numeric(BWA_NOR_gpx.tot$local_time) %% 86400, "extend"), origin = "1970-01-01", tz = "UTC")
BWA_NOR_gpx.tot <- BWA_NOR_gpx.tot[, c("lon", "lat", "time")]

#GET KEN_LAM as it has no.gpx but abundant RSO data (~5000)
KEN_Lam_gpx.tot <- readxl::read_xlsx("Z:/GEC/raw_data/KEN_Lam/normalized/KEN_Lam_elephants.xlsx", col_names = T, skip = 1, sheet = "RSO_flight")
KEN_Lam_gpx.tot$time <- as.POSIXct(as.numeric(KEN_Lam_gpx.tot$local_date)+zoo::na.fill(as.numeric(KEN_Lam_gpx.tot$local_time) %% 86400, "extend"), origin = "1970-01-01", tz = "UTC")
KEN_Lam_gpx.tot <- KEN_Lam_gpx.tot[, c("lon", "lat", "time")]

#GET XWA_TBC as it has no.gpx but some RSO data (~400), time is missing so it is filled with seconds of day equal to 08:00 am
XWA_TBC_gpx.tot <- readxl::read_xlsx("Z:/GEC/raw_data/XWA_TBC/normalized/XWA_TBC_elephants.xlsx", col_names = T, skip = 1, sheet = "RSO")
XWA_TBC_gpx.tot$time <- as.POSIXct(as.numeric(XWA_TBC_gpx.tot$local_date)+86400/3, origin = "1970-01-01", tz = "UTC")
XWA_TBC_gpx.tot <- XWA_TBC_gpx.tot[, c("lon", "lat", "time")]



#HARMONIZE THE DATES AND COMBINE THEM INTO A SINGLE .GPX-FILE
rm(all.gpx.tot)
all.gpx.tot <- mget(ls(pattern = "gpx.tot"))
SC <- rep(names(all.gpx.tot), lapply(all.gpx.tot, function(x) NROW(x)))
all.gpx.tot <- do.call(rbind, all.gpx.tot)
all.gpx.tot$SC <- gsub("_gpx.tot", "", SC)
all.gpx.tot$time <- gsub("T", " ", all.gpx.tot$time)
all.gpx.tot$time <- gsub("Z", "", all.gpx.tot$time)
all.gpx.tot$time <- gsub(" UTC", "", all.gpx.tot$time)
all.gpx.tot$time <- gsub(".000", "", all.gpx.tot$time)
all.gpx.tot$time[which(stringi::stri_length(all.gpx.tot$time) != 19)] <- NA
all.gpx.tot$time <- as.POSIXct(all.gpx.tot$time, tz = "UTC")

#SOME DATETIMES HAVE A DIFFERENT THAT HAS TO BE DEALT WITH DIFFERENTLY
SC.with.missing.dates <- unique(all.gpx.tot$SC[is.na(all.gpx.tot$time)])
AGO_Lui.times <- as.POSIXct(as.numeric(AGO_Lui_gpx.tot$time[match(rownames(all.gpx.tot[all.gpx.tot$SC == "AGO_Lui" & is.na(all.gpx.tot$time), ]), paste0("AGO_Lui_gpx.tot.", rownames(AGO_Lui_gpx.tot)))]), origin = "1970-01-01", tz = "UTC")
BWA_NOR.times <- as.POSIXct(BWA_NOR_gpx.tot$time[match(rownames(all.gpx.tot[all.gpx.tot$SC == "BWA_NOR" & is.na(all.gpx.tot$time), ]), paste0("BWA_NOR_gpx.tot.", rownames(BWA_NOR_gpx.tot)))])
COD_GAR.times <- as.POSIXct(COD_GAR_gpx.tot$time[match(rownames(all.gpx.tot[all.gpx.tot$SC == "COD_GAR" & is.na(all.gpx.tot$time), ]), paste0("COD_GAR_gpx.tot.", rownames(COD_GAR_gpx.tot)))])
COD_VIR.times <- as.POSIXct(as.numeric(COD_VIR_gpx.tot$time[match(rownames(all.gpx.tot[all.gpx.tot$SC == "COD_VIR" & is.na(all.gpx.tot$time), ]), paste0("COD_VIR_gpx.tot.", rownames(COD_VIR_gpx.tot)))]), origin = "1970-01-01", tz = "UTC")
KEN_TSV.times <- as.POSIXct(as.numeric(KEN_TSV_gpx.tot$time[match(rownames(all.gpx.tot[all.gpx.tot$SC == "KEN_TSV" & is.na(all.gpx.tot$time), ]), paste0("KEN_TSV_gpx.tot.", rownames(KEN_TSV_gpx.tot)))]), origin = "1970-01-01", tz = "UTC")
KEN_Lam.times <- KEN_Lam_gpx.tot$time
XWA_TBC.times <- XWA_TBC_gpx.tot$time
gpx.list <- split.data.frame(all.gpx.tot, all.gpx.tot$SC)


for(SC in SC.with.missing.dates){
  df <- gpx.list[[SC]]
  times <- df$time
  times[is.na(times)] <- as.POSIXct(get(paste0(SC, ".times")), tz = "UTC", origin = "1970-01-01")
  df$time <- times
  gpx.list[[SC]] <- df
  }

all.gpx.tot <- do.call(rbind, gpx.list)
rm(gpx.list)

#NOW ALL POINT DATA THAS HAS BEEN COLLECTED DURING THE GEC HAS BEEN COMBINED INTO ONE DATAFRAME (all.gpx.tot) WITH DATETIMES - THOSE CAN NOW BE USED TO ASSIGN DATES TO THE TRANSECT SEGMENTS.

segments$gpx_time <- NA
segments$TS <-unlist(lapply(strsplit(as.character(segments$ID), "p"), function(x) x[1]))

#this loop iterates through the transects of the single sites and assigns the datetimes of the points that are closest to the transect, by calculating the median of the datetimes of the closest points
for(SC in casefold(as.character(unique(segments$SC)))){
  gpx <- all.gpx.tot[casefold(all.gpx.tot$SC) == SC, ]
  for(TS in unique(segments$TS[which(casefold(as.character(segments$SC)) == SC)])){
    SC <- as.character(segments$SC[segments$TS == TS][1])
    if(NROW(gpx) != 0){
    ts_cents <- geosphere::centroid((segments[segments$TS == TS,]))
    times <- as.POSIXct(rep(0, NROW(ts_cents)), origin = "1970-01-01", tz = "UTC")
    for(i in 1:length(times)){
      dist <- geosphere::distGeo(ts_cents[i,], gpx[,c("lon", "lat")])
      times[i] <- gpx$time[which.min(dist)]
    }
    segments$gpx_time[casefold(as.character(segments$TS)) == TS] <- median(times)
    
    }
    print(TS)
  }
}

#check how the assigned dates look like and if they produce consistent results
for(i in unique(segments$SC)) plot(as.POSIXct(segments$gpx_time[segments$SC == i], tz = "UTC", origin = "1970-01-01"), main = i)

#are there any NAS (should not be the case!!)
unique(segments$SC[which(is.na(segments$gpx_time))])


segments$time <- segments$gpx_time/1e+9
segments_GEE <- segments[, c("ID", "time")]

writeOGR(segments_GEE, "segments_GEE.shp", "segments_GEE.shp", "ESRI Shapefile")
