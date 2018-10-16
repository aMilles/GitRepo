library(googledrive)
library(rgdal)
library(zoo)
library(stringi)
library(pscl)
library(MASS)
library(corrplot)

rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- rgdal::readOGR("Z:/GEC/segments.shp")


#Google Drive does not overwrite old predictors automatically (versioning), so make sure to download latest update of the predictors
preds <- drive_ls("GEC")
pred.choice <- vector(length = NROW(unique(preds$name)))

for(i in seq(length(pred.choice))){
  pred  <- data.frame(unique(preds[,1]))[i,]
  subset <- which(preds[,1] == pred)
  times <- unlist(lapply(preds$drive_resource, function(x) x$createdTime))[subset]
  times <- as.POSIXct(times, format = "%FT%T")
  pred.choice[i] <- subset[which.max(times)]
}

#download and add predictors to the global environment (including predictor and response from extract_COUNT-HT.R)
for(i in pred.choice) drive_download(path = paste0("Z:/predictors/", preds[i,1]), file = preds[i,], overwrite = T)
for(pred in list.files("Z:/predictors")) assign(gsub(".csv", "", pred), read.csv(paste0("Z:/predictors/", pred)))
HT <- HT[,2:3]
COUNT <- COUNT[,2:3]

#tif of protected areas does only contain 1 and NA. replace NAs with 0 in protected areas predictor
PA$max[is.na(PA$max)] <- 0

#combine all predictors + count data, do some unit conversions if necessary
preds <- gsub(".csv", "", list.files("Z:/predictors"))
preds.list <- mget(preds)
all.preds <- Reduce(function(d1, d2) merge(d1, d2, by = "ID", all.x = T, all.y = FALSE), 
                    preds.list) 
names(all.preds) <- c("ID", preds)
all.preds$PA <- as.factor(all.preds$PA)
summary(all.preds)

all.preds <- all.preds[match(as.character(segments$ID), as.character(all.preds$ID)), ]
all.preds$Site <- segments$SC
all.preds$Country <- segments$CC
all.preds$Transect <- unlist(lapply(strsplit(as.character(all.preds$ID), "p"), function(x) x[1]))
all.preds$REPS[is.na(all.preds$REPS)] <- 0

all.preds$HT <- as.character(all.preds$HT)
all.preds$HT <- ifelse(is.na(all.preds$HT), "none", all.preds$HT)
all.preds <- na.omit(all.preds)

all.preds$WA <- all.preds$WA / 1000 # scale to km
all.preds$AR <- all.preds$PI; all.preds$PI <- NULL
all.preds$SS <- all.preds$`NA`; all.preds$`NA` <- NULL
all.preds$AD <- all.preds$AI; all.preds$AI <- NULL
all.preds$TC <- all.preds$TC - 273.15 # scale to °C
all.preds$VD <- all.preds$VD / 10000 # scale to 0 to 1

#save all predictors and the response
write.csv(all.preds, "C:/Users/amilles/Dropbox/modelling/yxtable.csv")