library(googledrive)
library(rgdal)
library(zoo)
library(stringi)
library(pscl)
library(MASS)
library(corrplot)
library(reshape2)
library(ggplot2)
library(reshape2)
library(raster)
library(velox)
setwd("Z:/")

dropbox.file = "C:/Users/amilles/Dropbox/"
xy <- read.csv(paste0(dropbox.file, "modelling/yxtable.csv"))[,-c(1)]
segments <- readOGR(paste0(dropbox.file, "modelling/segments_GEE.shp"))
files <- paste0(dropbox.file, "/modelling/CC/", list.files(paste0(dropbox.file, "/modelling/CC")))

lc <- stack("Z:/IMAGE/PBL_IMAGE_SSP2_SPA2_RCP_6.0_categorial_land_cover_type_5min.tif")


v2015 <- velox(lc$PBL_IMAGE_SSP2_SPA2_RCP_6.0_categorial_land_cover_type_5min.1)
v2100 <- velox(lc$PBL_IMAGE_SSP2_SPA2_RCP_6.0_categorial_land_cover_type_5min.3)

ext_2015 <- v2015$extract(segments, small = T)
ext_2100 <- v2100$extract(segments, small = T)


class.change <- function(ext.2015 = ext_2015, ext.2100 = ext_2100, class, .xy = xy){
  out_2100 <- out_2015 <- vector("list", length = length(unique(.xy$Site)))
  names(out_2015) <- names(out_2100) <- unique(.xy$Site)
  for(i in unique(.xy$Site)){
    exts <- match(as.character(.xy$ID[.xy$Site == i]), segments$ID)
    out_2015[[i]] <- sum(do.call(c, ext.2015[exts]) == class)
    out_2100[[i]] <- sum(do.call(c, ext.2100[exts]) == class)
  }
  ref <- do.call(c, out_2015)
  ref[is.na(ref)] <- 0
  return((do.call(c, out_2100) - ref)/  ref * 100)
}

out <- vector("list", length = 20)
for(i in seq(20)) out[[i]] <- class.change(class = i)
out[[1]]

install.packages('xlsx')
wpp <- xlsx::read.xlsx("C:/Users/amilles/Dropbox/modelling/WPP2017_POP_F06_POPULATION_DENSITY.xlsx", sheetName = "MEDIUM VARIANT", startRow =  18, header = T)


countries <- c("Angola", 
  "Botswana", 
  "Zimbabwe", 
  "Kenya", 
  "Ethiopia",
  "Democratic Republic of the Congo",
  "Chad",
  "Benin", 
  "Burkina Faso",
  "Niger")

pop_2100 <- wpp[match(countries, wpp$Region..subregion..country.or.area..),]$X2100
pop_2015 <- wpp[match(countries, wpp$Region..subregion..country.or.area..),]$X2015


pop_growth <- data.frame(growth = 100 * pop_2100 / pop_2015)
pop_growth$country <- c("AGO", 
  "BWA", 
  "ZWE", 
  "KEN", 
  "ETH",
  "COD",
  "TCD",
  "XWA", 
  "XWA",
  "XWA")
pop_growth

pop_growth <- aggregate(growth ~ country, pop_growth, mean)

agriculture <- data.frame(agriculture = out[[1]])
agriculture$site <- row.names(agriculture)
agriculture$agriculture[is.na(agriculture$agriculture)] = 0
agriculture$country <- do.call(c, lapply(strsplit(agriculture$site, "_"), function(x) x[1]))

lc_change <- merge(agriculture, pop_growth)

write.csv(lc_change, "C:/Users/amilles/Dropbox/modelling/lc_change.R")

