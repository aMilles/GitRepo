library(googledrive)
library(rgdal)
library(zoo)
library(stringi)
library(pscl)
library(MASS)
library(corrplot)
library(reshape2)
library(ggplot2)
setwd("Z:/")

dropbox.file = "C:/Users/amilles/Dropbox/"
xy <- read.csv(paste0(dropbox.file, "modelling/yxtable.csv"))[,-c(1)]
segments <- readOGR("Z:/GEC/segments_GEE.shp")

#make sure to download latest update of the predictor
preds <- drive_ls("5years")
pred.choice <- vector(length = NROW(unique(preds$name)))

for(i in seq(length(pred.choice))){
  pred  <- data.frame(unique(preds[,1]))[i,]
  subset <- which(preds[,1] == pred)
  times <- unlist(lapply(preds$drive_resource, function(x) x$createdTime))[subset]
  times <- as.POSIXct(times, format = "%FT%T")
  pred.choice[i] <- subset[which.max(times)]
}
#download and read predictors
#for(i in pred.choice) drive_download(path = paste0("cc_model/", preds[i,1]), file = preds[i,], overwrite = T)
for(pred in list.files("cc_model/")) assign(gsub(".csv", "", pred), read.csv(paste0("cc_model/", pred)))
preds <- stringi::stri_replace_all_regex(unique(preds$name), ".csv", "")

names(IMAGE_2015) <- names(IMAGE_2100)
IMAGE <- rbind(IMAGE_2015, IMAGE_2100)
rows <- as.character(IMAGE$ID) %in% as.character(xy$ID)
IMAGE <- IMAGE[rows,]
IMAGE$Site <- xy$Site[match(as.character(xy$ID), as.character(IMAGE$ID))]
IMAGE$year <- rep(as.factor(c(2015, 2100)), each = nrow(IMAGE)/2)
mIMAGE <- melt(IMAGE[,c("aridity", "leaves2100", "npp", "Site", "year")], id.vars = c("year", "Site"))


gg <- ggplot(mIMAGE, aes(x = year, y = value))+
  geom_boxplot()+
  facet_grid(variable ~ Site, scales = "free")
mIMAGE$variable[]

pdf("test.pdf", width = 30, height = 6)
  gg
dev.off()

for(pred in preds) assign(paste0(pred,"_med", apply(get(pred)[,-1], 1, function(x) median(x, na.rm = T))))

months <- seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "month")
WA$month <- VD$month <-  rep(months, each = nrow(WA)/12)
VDgg <- melt(VD, id.vars = c("ID", "month"))
VDmonth <- aggregate(VDgg$value, by = list(VDgg$ID, VDgg$month), FUN = median)
VDgg <- aggregate(VDgg$value, by = list(VDgg$ID), FUN = sum)
VDgg$x <- VDgg$x / 5
VDgg <- VDgg[match(as.character(segments$ID),as.character(VDgg$Group.1)),]

WAgg <- melt(VD, id.vars = "ID")
WAgg <- aggregate(WAgg$value, by = list(WAgg$ID), FUN = mean)
WAgg <- WAgg[match(as.character(segments$ID),as.character(WAgg$Group.1)),]

WAdf <- data.frame(WA = WAgg$x, IMAGE_2015[,-1])
VDdf <- data.frame(VD = VDgg$x, IMAGE_2015[,-1])

rF <- ranger::ranger(VD ~., data = VDgg)




months <- levels(unique(lubridate::month(as.POSIXct(seq(1, 86400*365, by = 86400), origin = "1970-01-01"), label = T, locale = "English")))
summary(VD)

unify <- function(df, select.Site = "ZWE_MAT", select.Month = 1:12, select.Year = 1:5){
  names(df)[-1] <- paste0("y", 2010:2014)
  months <- seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "month")
  df$month <- rep(months, each = nrow(WA)/12)
  df$ID <- as.character(df$ID)
  xy$ID <- as.character(xy$ID)
  df$Site <- xy$Site[match(df$ID, xy$ID)]
  df$Country <- xy$Country[match(df$ID, xy$ID)]
  df <- na.omit(df)
  df <- df[df$Site == select.Site & df$month %in% months[select.Month], paste0("y", 2010:2014)[select.Year]]
  return(df)
}
VD$month <- NULL
SC$month <- NULL

TC_list <- WA_list <- SC_list <- vector("list", 5)
for(i in 1:5){
  TC_list[[i]] <- apply(data.frame(matrix(unify(TC, select.Year =  i, select.Month = 1:5, select.Site = "KEN_LAI"), nc = 5)), 2, median)
  SC_list[[i]] <- apply(data.frame(matrix(unify(SC, select.Year =  i, select.Month = 1:5, select.Site = "KEN_LAI"), nc = 5)), 2, median)
  WA_list[[i]] <- median((unify(VD, select.Month = 5, select.Year =  i, select.Site = "KEN_LAI")))
} 

for(i in 1:5){
  TC_list[[i]] <- data.frame(matrix(unify(TC, select.Year =  i, select.Month = 1:5, select.Site = "KEN_LAI"), nc = 5))
  SC_list[[i]] <- data.frame(matrix(unify(SC, select.Year =  i, select.Month = 1:5, select.Site = "KEN_LAI"), nc = 5))
  WA_list[[i]] <- (unify(VD, select.Month = 5, select.Year =  i, select.Site = "KEN_LAI"))
} 
unique(xy$Site)
sum_SC <- apply(do.call(rbind, SC_list), 1, summary)
sum_TC <- apply(do.call(rbind, TC_list), 1, summary)
median_WA <- apply(do.call(cbind, WA_list), 1, median)
model.df <- data.frame(WA = do.call(c, WA_list), do.call(rbind, TC_list), do.call(rbind, SC_list), t(sum_SC), t(sum_TC))


lines(WA_list[[1]])
 f.gam <- as.formula(paste0("WA ~", paste0("s(",names(model.df)[-1], ")", collapse = "+")))
 rF <- mgcv::gam(f.gam, data = model.df)


rF <- ranger::ranger(WA ~ . , data = model.df)

pr <- NULL
actual <- NULL
for(i in seq(5) - 1){
  rows <- seq(nrow(model.df)/5) + i * nrow(model.df)/5
  rF <- ranger::ranger(WA ~ . , data = model.df[- rows,])
  p <- predict(rF, data  = model.df[rows,])
  print(plot(p$predictions, model.df$WA[rows]))
  print(cor(p$predictions, model.df$WA[rows]))
  print(sum(p$predictions))
  pr <- append(pr, sum(p$predictions))
  actual <- append(actual, sum(model.df$WA[rows]))
  print(sum(model.df$WA[rows]))
}

plot(pr, actual)
plot(rF$predictions, actual)
SC_ss <- unify(SC, select.Year =  2, select.Month = 1:5)
SC_ss <- matrix(SC_ss, nc = 9)
TC_ss <- unify(TC, select.Year =  2, select.Month = 1:9)
TC_ss <- matrix(TC_ss, nc = 9)
model.df <- data.frame(SC_ss, TC_ss)

p <- predict(rF, data  = model.df)
plot(p$predictions, unify(WA, select.Month = 10, select.Year =  2))

glm.()

plot(SC_ss$y2010 ~ WA_ss$y2010)

df.agg <- aggregate(y2010 ~ Site + month, FUN = median, data = df)

names(VD)[-1] <- paste0(2010:2014)
VD$month <- rep(seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "month"), each = nrow(WA)/12)
WAgg <- melt(VD, id.vars = c("ID", "month"))
WAgg$ID <- as.character(WAgg$ID)
xy$ID <- as.character(xy$ID)

WAgg$Site <- xy$Site[match(WAgg$ID, xy$ID)]
WAgg$Country <- xy$Country[match(WAgg$ID, xy$ID)]

WAgg.agg <- aggregate(value ~ Site + month + variable, FUN = median, data = WAgg)


ggplot(WAgg.agg, aes(x = month, y = value, group = variable, col = variable))+
  geom_line()+
  facet_wrap(~Site)




names(SC)[-1] <- paste0(2010:2014)
SC$month <- rep(seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "month"), each = nrow(TC)/12)
SCgg <- melt(SC, id.vars = c("ID", "month"))
SCgg$ID <- as.character(SCgg$ID)
xy$ID <- as.character(xy$ID)

SCgg$Site <- xy$Site[match(SCgg$ID, xy$ID)]
SCgg$Country <- xy$Country[match(SCgg$ID, xy$ID)]

SCgg.agg <- aggregate(value ~ Site + month + variable, FUN = median, data = SCgg)


ggplot(SCgg.agg, aes(x = month, y = value, group = variable, col = variable))+
  geom_line()+
  facet_wrap(~Site)

