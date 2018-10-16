library(googledrive)
library(rgdal)
library(zoo)
library(stringi)
library(pscl)
library(MASS)
library(corrplot)
setwd("Z:/")

#make sure to download latest update of the predictor
preds <- drive_ls("Seasonal Shift")
pred.choice <- vector(length = NROW(unique(preds$name)))

for(i in seq(length(pred.choice))){
  pred  <- data.frame(unique(preds[,1]))[i,]
  subset <- which(preds[,1] == pred)
  times <- unlist(lapply(preds$drive_resource, function(x) x$createdTime))[subset]
  times <- as.POSIXct(times, format = "%FT%T")
  pred.choice[i] <- subset[which.max(times)]
}
#download and read predictors
for(i in pred.choice) drive_download(path = paste0("2014_Season/", preds[i,1]), file = preds[i,], overwrite = T)
for(pred in list.files("2014_Season/")) assign(gsub(".csv", "", pred), read.csv(paste0("2014_Season/", pred)))

preds <- stringi::stri_replace_all_regex(unique(preds$name), ".csv", "")
months <- levels(unique(lubridate::month(as.POSIXct(seq(1, 86400*365, by = 86400), origin = "1970-01-01"), label = T, locale = "English")))

xy <- read.csv("C:/Users/amilles/Dropbox/modelling/yxtable.csv")
xy$ID <- as.character(xy$ID)
xy_future <- do.call(rbind, mget(rep("xy", 12)))
xy_future$month <- rep(months, each = nrow(xy))
xy_future <- xy_future[, - match(preds, names(xy))]

which(xy_future$ID == "ts1p1")

for(pred in preds){
  var <- get(pred)
  var$ID <- as.character(var$ID)
  var <- var[var$ID %in% xy$ID,]
  var$month <- rep(months, each = nrow(xy))
  names(var)[2] <- pred
  assign(pred, var)
}

all.preds <- do.call(cbind, lapply(mget(preds), function(x) data.frame(x[,2])))
all.preds$ID_adv <- paste0(var$ID, var$month)
names(all.preds)[1:4] <- preds
xy_future$ID_adv <- paste0(xy_future$ID, xy_future$month)

xy_future[,preds] <- all.preds[match(xy_future$ID_adv, all.preds$ID_adv),preds]
xy_future$VD <- 10000 * xy_future$VD
xy_future$SC <- 3600 * 3 * xy_future$SC


table(xy_future$Site[is.na(xy_future$VD)])

write.csv(xy_future, "C:/Users/amilles/Dropbox/modelling/yxtable_season.csv")

gg_season <- data.frame(xy_future[,c("ID", "Site", "Country", "month", preds)])
gg_season <- reshape2::melt(gg_season, id.vars = c("ID", "Site", "Country", "month"))
gg_season$month <- factor(gg_season$month, levels = months)

gg_season.agg <- aggregate(gg_season$value, by = list(gg_season$Site, gg_season$month, gg_season$variable, gg_season$Country), FUN = function(x) quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = T))
gg_season.agg <- cbind(gg_season.agg[,-5], data.frame(gg_season.agg$x))
names(gg_season.agg) <- c("site", "month", "predictor", "country", "q2.5", "q50", "q97.5")

colors <- viridis::cividis(length(unique(gg_season$Site)))
colors <- sample(colors, length(colors))
library(ggplot2)
season.plot <- ggplot(gg_season.agg, aes(x = month, y = q50, ymin = q2.5, ymax = q97.5, col = site, fill = site, group = site))+
  geom_ribbon(alpha = .5)+
  geom_line()+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  facet_grid(predictor ~ country, scales = "free")+
  xlab("month")+
  ylab("predictor range")+
  scale_x_discrete(breaks=levels(gg_season.agg$month)[c(2,7,11)])+
  theme_bw()

tf_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/transform_sheet.csv", stringsAsFactors = F)
scale_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/scale_sheet.csv", stringsAsFactors = F)

for(pred in seq(length(tf_sheet$name))){
  
  name = tf_sheet$name[pred]
  shift = tf_sheet$shift[pred]
  tf = tf_sheet$value[pred]
  
  scale = scale_sheet$scale[pred]
  center = scale_sheet$center[pred]
  
  if(is.na(name)) name <- "NA."
  
  cov <- xy_future[,name]
  
  cov <- cov - shift
  cov <- cov^tf
  cov <- (cov - center) / scale
  
  print(summary(cov))
  
}

xy_season <- xy[,c("ID", "Site", "Country")]



old_Country <- "Start"
site_number <- data.frame(Site = NA, count = NA)
for(Site in as.character(unique(xy$Site))){
  
  Country = strsplit(Site, "_")[[1]][1]
  
  if(Country != old_Country){ 
    count = 1
  }else{
    count = count + 1
  } 
  
  site_number <- rbind(site_number, data.frame(Site = Site, count = count))
  old_Country = Country
  
}

site_number <- site_number[-1,]

all.preds <- do.call(rbind, mget(preds))
xy_season$ID <- as.character(xy_season$ID)
all.xy <- base::merge(xy_season, all.preds, by = "ID")
gg_season <- reshape2::melt(all.xy, id.vars = c("Site", "pred", "Country", "ID"))
gg_season$value[gg_season$pred == "SC"] <- gg_season$value[gg_season$pred == "SC"] * 3 *3600
gg_season$value[gg_season$pred == "WA"] <- gg_season$value[gg_season$pred == "WA"] / 1000
gg_season$variable <- rep(seq.Date(as.Date("2014-01-01"), as.Date("2014-12-31"), by= "month"), each = nrow(gg_season)/12)
gg_season <- merge(gg_season, site_number, by = "Site")

pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/seasonality_predictors.pdf")
ggplot(gg_season, aes(group = Site, x = variable, y = value, color = Site))+
  geom_smooth(alpha = .1, se = F)+
  facet_grid(pred~Country, scales = "free")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b")+
  xlab("month")+
  ylab("predictor value")+
  scale_color_manual(values = sample(viridis::viridis(16), 16))+
  guides(linetype = guide_legend(title = "Site Number"))+
  theme_light()+
  theme(legend.position = "top")
dev.off()






#make sure to download latest update of the predictor
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


TC_cc$year <- SC_cc$year <- rep(as.factor(c(1950, 2015, 2100)), each = nrow(SC_cc)/3)
TC_cc$month <- SC_cc$month <- rep(rep(months, each = nrow(SC_cc)/36))
SC_cc <- SC_cc[as.character(SC_cc$ID) %in% as.character(xy$ID),]
TC_cc <- TC_cc[as.character(TC_cc$ID) %in% as.character(xy$ID),]
all(as.character(SC_cc$ID) == as.character(TC_cc$ID))
plot(as.numeric(factor(SC_cc$month[which(SC_cc$ID == "ts1p1")], levels = months)))
TC_cc$Site <- SC_cc$Site <- xy$Country[match(as.character(TC_cc$ID), as.character(xy$ID))]

gg_SC <- reshape2::melt(SC_cc[,-1], id.vars = c("year", "month", "Site"))

agg.SC <- aggregate(gg_SC$value, by = list(as.factor(gg_SC$month), gg_SC$Site, gg_SC$year), FUN = quantile, prob = c(0.025, 0.5, 0.975))
agg.SC <- cbind(agg.SC[, -4], data.frame(agg.SC$x))
names(agg.SC) <- c("Month", "Site", "Year", "lower", "mid", "upper")
agg.SC$Month <- factor(agg.SC$Month, levels = months)
#agg.SC$Median <- agg.SC$Median - agg.SC$Median[agg.SC$Year == "2015"]


gg_TC <- reshape2::melt(TC_cc[,-1], id.vars = c("year", "month", "Site"))

agg.TC <- aggregate(gg_TC$value - 273.15, by = list(as.factor(gg_TC$month), gg_TC$Site, gg_TC$year), FUN = quantile, prob = c(0.025, 0.5, 0.975))
agg.TC <- cbind(agg.TC[, -4], data.frame(agg.TC$x))
names(agg.TC) <- c("Month", "Site", "Year", "lower", "mid", "upper")
agg.TC$Month <- factor(agg.TC$Month, levels = months)
#agg.SC$Median <- agg.SC$Median - agg.SC$Median[agg.SC$Year == "2015"]

pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/climate_change_temperature.pdf")
Temp <- ggplot()+
  geom_ribbon(data = agg.TC, aes(x = Month, ymin = lower, ymax = upper, group = Year, fill = Year), alpha = .5)+
  geom_line(data = agg.TC, aes(x = Month, y = mid, group = Year, col = Year), size = 1.3, alpha = .8)+
  facet_wrap(~Site)+
  ylab("daily maximum temperature [°C]")+
  scale_x_discrete(breaks=levels(agg.TC$Month)[c(2,7,11)])+
  theme_bw()+
  theme(text = element_text(size = 14), legend.position = "none")
dev.off()

pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/climate_change_precipitation.pdf")
Precip <- ggplot()+
  geom_ribbon(data = agg.SC, aes(x = Month, ymin = lower, ymax = upper, group = Year, fill = Year), alpha = .5)+
  geom_line(data = agg.SC, aes(x = Month, y = mid, group = Year, col = Year), size = 1.3)+
  facet_wrap(~Site)+
  ylab("4-month precipitation [mm]")+
  scale_x_discrete(breaks=levels(agg.SC$Month)[c(2,7,11)])+
  theme_bw()+
  theme(text = element_text(size = 14), legend.position = "top")
dev.off()



comb <- rbind(agg.TC, agg.SC)
comb$mode <- rep(c("daily Tmax [°C]", "4-month precipitation [mm]"), each = nrow(agg.SC))

pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/climate_change.pdf", width = 12, height = 7)
ggplot()+
  geom_ribbon(data = comb, aes(x = Month, ymin = lower, ymax = upper, group = Year, fill = Year), alpha = .5)+
  geom_line(data = comb, aes(x = Month, y = mid, group = Year, col = Year), size = 1.3, alpha = .8)+
  facet_grid(mode ~ Site, scales = "free_y")+
  ylab("")+
  scale_x_discrete(breaks=levels(agg.TC$Month)[c(3,10)])+
  theme_bw()+
  theme(text = element_text(size = 20), strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14), legend.position = "top")+
  guides(fill = guide_legend(title = "year"), color = guide_legend(title = "year"))+
  xlab("month")
dev.off()



gridExtra::grid.arrange(Precip, Temp)
