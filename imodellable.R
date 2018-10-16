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
setwd("Z:/")

dropbox.file = "C:/Users/amilles/Dropbox/"
xy <- read.csv(paste0(dropbox.file, "modelling/yxtable.csv"))[,-c(1)]
segments <- readOGR(paste0(dropbox.file, "modelling/segments_GEE.shp"))
files <- paste0(dropbox.file, "/modelling/CC/", list.files(paste0(dropbox.file, "/modelling/CC")))

for(file in files) assign(tools::file_path_sans_ext(basename(file)), read.csv(file))



names(IMAGE_2015) <- names(IMAGE_2100)
IMAGE <- rbind(IMAGE_2015, IMAGE_2100)
rows <- as.character(IMAGE$ID) %in% as.character(xy$ID)
IMAGE <- IMAGE[rows,]
IMAGE$Site <- xy$Site[match(as.character(xy$ID), as.character(IMAGE$ID))]
IMAGE$year <- rep(as.factor(c(2015, 2100)), each = nrow(IMAGE)/2)
mIMAGE <- melt(IMAGE[,c("aridity", "leaves2100", "npp", "Site", "year")], id.vars = c("year", "Site"))

IMAGE[1,]
IMAGE[(nrow(IMAGE)/2) + 1,]


gg <- ggplot(mIMAGE, aes(x = year, y = value))+
  geom_boxplot()+
  facet_grid(variable ~ Site, scales = "free")
mIMAGE$variable[]

pdf("test.pdf", width = 30, height = 6)
gg
dev.off()

for(pred in preds) assign(paste0(pred,"_med", apply(get(pred)[,-1], 1, function(x) median(x, na.rm = T))))


months <- seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "month")
TC$month <- SC$month <- WA$month <- VD$month <-  rep(months, each = nrow(WA)/12)
VDgg <- melt(VD, id.vars = c("ID", "month"))
VDmonth <- aggregate(VDgg$value, by = list(VDgg$ID, VDgg$month), FUN = median)
VDgg <- aggregate(VDgg$value, by = list(VDgg$ID), FUN = sum)
VDgg$x <- VDgg$x / 60
VDgg <- VDgg[as.character(VDgg$Group.1) %in% xy$ID, ]
names(VDgg) <- c("ID", "VD")
VDgg <- merge(VDgg, xy[,c("ID", "Site")], by = "ID")

WAgg <- melt(WA[,1:5], id.vars = "ID")
WAgg <- aggregate(WAgg$value, by = list(WAgg$ID), FUN = mean)
WAgg <- WAgg[as.character(WAgg$Group.1) %in% xy$ID, ]
names(WAgg) <- c("ID", "VD")
WAgg <- merge(WAgg, xy[,c("ID", "Site")], by = "ID")


TC_cc <- read.csv("Z:/climate_change/TC_cc.csv")
SC_cc <- read.csv("Z:/climate_change/SC_cc.csv")
r <- nrow(TC_cc)
TC_2100 <- TC_cc[((2 * r/3)+1) : r, ]
SC_2100 <- SC_cc[((2 * r/3)+1) : r, ]
TC_2015 <- TC_cc[((1 * r/3)+1) : (2 * r/3), ]
SC_2015 <- SC_cc[((1 * r/3)+1) : (2 * r/3), ]


SCmonth <- aggregate(SC_2015$mean, by = list(SC$ID), FUN = sum)
SCmonth <- SCmonth[SCmonth$Group.1 %in% VDgg$ID,]
SCmonth <- SCmonth[match(as.character(SCmonth$Group.1), as.character(VDgg$ID)), ]
SC_2015 <- SCmonth$x
SCmonth <- aggregate(SC_2100$mean, by = list(SC$ID), FUN = sum)
SCmonth <- SCmonth[SCmonth$Group.1 %in% VDgg$ID,]
SCmonth <- SCmonth[match(as.character(SCmonth$Group.1), as.character(VDgg$ID)), ]
SC_2100 <- SCmonth$x

TCmonth <- aggregate(TC_2015$mean, by = list(TC$ID), FUN = mean)
TCmonth <- TCmonth[as.character(TCmonth$Group.1) %in% as.character(VDgg$ID),]
TCmonth <- TCmonth[match(as.character(TCmonth$Group.1), as.character(VDgg$ID)), ]
TC_2015 <- TCmonth$x
TCmonth <- aggregate(TC_2100$mean, by = list(TC$ID), FUN = mean)
TCmonth <- TCmonth[as.character(TCmonth$Group.1) %in% as.character(VDgg$ID),]
TCmonth <- TCmonth[match(as.character(TCmonth$Group.1), as.character(VDgg$ID)), ]
TC_2100 <- TCmonth$x



WAdf_2015 <- data.frame(merge(WAgg, IMAGE_2015[,-1], by = "ID"), TC = TC_2015, SC = SC_2015)
WAdf_2100 <- data.frame(merge(WAgg, IMAGE_2100[,-1], by = "ID"), TC = TC_2100, SC = SC_2100)
VDdf_2015 <- data.frame(merge(VDgg, IMAGE_2015[,-1]), TC = TC_2015, SC = SC_2015)
VDdf_2100 <- data.frame(merge(VDgg, IMAGE_2100[,-1], by = "ID"), TC = TC_2100, SC = SC_2100)

boxplot(TC ~ Site, data = WAdf_2015)
boxplot(Branches_2100 ~ Site, data = WAdf_2100)

change <- data.frame(WAdf_2015$Site)
names(change) <- "site"
change$aridity <- 100 * ((WAdf_2100$aridity - WAdf_2015$aridity)/WAdf_2015$aridity)
change$npp <- 100 * ((WAdf_2100$npp - WAdf_2015$npp)/WAdf_2015$npp)
change$SC <- 100 * ((WAdf_2100$SC - WAdf_2015$SC)/WAdf_2015$SC)
change$TC <- 100 * ((WAdf_2100$TC - WAdf_2015$TC)/(WAdf_2015$TC - 273.15))
ggchange <- melt(change, id.vars = "site")
ggchange <- aggregate(value ~ variable + site, ggchange, median)

lc_change <-read.csv("C:/Users/amilles/Dropbox/modelling/lc_change.R", row.names = 1)[,-1]
ggchange <- melt(merge(change, lc_change), id.vars = "site")
ggchange <- aggregate(value ~ variable + site, ggchange, median)
levels(ggchange$variable)

weird <- scales::trans_new("signed_log",
                           transform=function(x) sign(x)*log10(abs(x)),
                           inverse=function(x) sign(x)*(10^abs(x)))
10^log10(20)
pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/environmental_change.pdf", width = 10, height = 5)
ggplot(ggchange, aes(x = site, y = value, fill = variable))+
  geom_bar(stat="identity",position="dodge", col = "black")+
  ylab("median change [%]")+
  guides(fill = guide_legend(title = "", nrow = 1))+
  scale_fill_manual(labels=c( "humdity", "NPP", "precipitation", "temperature", "agriculture", "human density"), values = c("turquoise", "springgreen4", "royalblue2", "tomato", "sienna4", "cornsilk2"))+
theme_bw()+
  theme(text = element_text(size = 20), legend.position = "top", panel.grid.major.x = element_line(color = c("gray90", NA), size = 14), panel.grid.major.y = element_line(color = c("gray90"), linetype = "solid"), panel.grid.minor.y =  element_line(color = NA) ,axis.text.x = element_text(angle = 30, hjust = .9))+
  scale_y_continuous(trans = weird, breaks = c(-100, -50, -20, -10, - 5, 0, 5, 20, 10, 50, 100, 200, 500, 1000), limits = c(-100, 1000))
dev.off()



VDdf_2015$
VDdf_2100$pred_2100 <- NA
VDdf_2100$pred_2015 <- NA

for(Site in unique(VDdf_2015$Site)){
  nr <- which(VDdf_2015$Site == Site)
  rF <- ranger::ranger(VD ~ ., data = VDdf_2015[-nr, -1])
  VDdf_2100$pred_2100[nr] <- predict(rF, VDdf_2100[nr, -1])$predictions
  VDdf_2100$pred_2015[nr] <- predict(rF, VDdf_2015[nr, -1])$predictions
}

plot(`pred_2015` ~ VD, VDdf_2100)

VDchange <- cbind(aggregate(VD ~ Site, data = VDdf_2100, FUN = mean),
aggregate(pred_2015 ~ Site, data = VDdf_2100, FUN = mean)$pred_2015,
aggregate(pred_2100 ~ Site, data = VDdf_2100, FUN = mean)$pred_2100)
names(VDchange) <- c("Site", "originalVD", "2015_prediction", "2100_prediction")
VDchange
plot(`2100_prediction` ~ originalVD, VDchange)
head(WAdf_2100)
rF <- ranger::ranger(VD ~., data = WAdf_2015)
WA_pred <- predict(rF, WAdf_2100)
WAdf_2100$pred_2100 <-  WA_pred$predictions
WAdf_2100$pred_2015 <- rF$predictions
WAchange <- cbind(aggregate(VD ~ Site, data = WAdf_2100, FUN = mean),
                  aggregate(pred_2015 ~ Site, data = WAdf_2100, FUN = mean)$pred_2015,
                  aggregate(pred_2100 ~ Site, data = WAdf_2100, FUN = mean)$pred_2100)
names(WAchange) <- c("Site", "originalWA", "2015_prediction", "2100_prediction")

WAchange


names(VDchange) <- c("Site", "originalVD", "2015_prediction", "2100_prediction")
VDchange


plot(rF$predictions, VDdf$VD)
names(VDdf)[1] <- "VDresp"
apply(VDdf, 1, mean)


for(Site in unique(VDdf$Site)){
  rows <- which(VDdf$Site == Site)
  
  rF <- glm(VD ~., data = VDdf[-rows, c(-1, -3)])
  #p <- predict(rF, VDdf[rows,])$predictions
  p <- predict(rF, newdata = VDdf[rows,])
  plot(p, VDdf$VD[rows])
  print(Site)
  #print(cor(p, VDdf$VD[-rows]))
  print(median(p))
  print(median(VDdf$VD[rows]))
  print(mean(p))
  print(mean(VDdf$VD[rows]))
}



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
SC_ss <- unify(SC, select.Year =  2, select.Month = 1:%)
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

plot(VD$X2010_02_iterator_0)
