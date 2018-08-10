library(googledrive)
library(rgdal)
library(zoo)
library(stringi)
library(pscl)
library(MASS)
library(corrplot)

rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- rgdal::readOGR("Z:/GEC/segments.shp")


#make sure to download latest update of the predictor
preds <- drive_ls("GEC")
pred.choice <- vector(length = NROW(unique(preds$name)))

for(i in seq(length(pred.choice))){
  pred  <- data.frame(unique(preds[,1]))[i,]
  subset <- which(preds[,1] == pred)
  times <- unlist(lapply(preds$drive_resource, function(x) x$createdTime))[subset]
  times <- as.POSIXct(times, format = "%FT%T")
  pred.choice[i] <- subset[which.max(times)]
}
preds
#download and read predictors
for(i in pred.choice) drive_download(path = paste0("Z:/predictors/", preds[i,1]), file = preds[i,], overwrite = T)
for(pred in list.files("Z:/predictors")) assign(gsub(".csv", "", pred), read.csv(paste0("Z:/predictors/", pred)))
HT <- HT[,2:3]
COUNT <- COUNT[,2:3]
summary(TC)
#replace NA with 0 in protected areas predictor
PA$max[is.na(PA$max)] <- 0
summary(PA)

segments$gpx_time <- as.POSIXct(segments$gpx_time, origin = "1970-01-01", tz = "UTC")
data.frame(segments[segments$ID %in% SC$ID[is.na(SC$mean)],])
data.frame(segments[segments$ID %in% VD$ID[is.na(VD$mean)],])
data.frame(segments[segments$ID %in% NB$ID[is.na(NB$mean)],])
data.frame(segments[segments$ID %in% NB$ID[is.na(NB$mean)],])

summary(NB)
summary(VD)
summary(WA)

#TC <- read.csv("Z:/predictors/TC.csv")
TC$approx <- NA
summary(unique(TC$ID))
for(seg in seq(nrow(segments))){
  SoDs = seq(0, 86400, length = 5)
  SoD = as.numeric(as.POSIXct(segments$gpx_time[seg], tz = "UTC", origin = "1970-01-01")) %% 86400
  diff2SoDs = abs(SoD - SoDs) / (86400 / 4)
  SoDs_next = head(order(diff2SoDs), 2)
  TC$approx[seg] <- (1 - diff2SoDs[SoDs_next[1]]) * TC[seg, SoDs_next[1] + 1] + 
    (1 - diff2SoDs[SoDs_next[2]]) * TC[seg, SoDs_next[2] + 1]
}
TC <- TC[,c("ID", "approx")]
TC$approx <- unlist(TC$approx)
#WA
WA <- aggregate(WA$mean, by = list(WA$ID), FUN = function(x) mean(x, na.rm = T))
names(WA) <- c("ID", "WA")

#HD
HD2 <- HD[, !names(HD) %in% c("system.index", "ID", "polyCent", "time", "population",".geo")]
HD <- data.frame("ID" = HD$ID, "sum" = apply(HD2[,-1], 1, sum))


#combine all predictors + count data
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



###

write.csv(all.preds, "Z:/modelling/yxtable.csv")

#scale non-factor predictors
#for(col in which(!names(all.preds) %in% c("ID", "CC", "HT", "PA", "COUNT"))) all.preds[,col] <- scale(all.preds[,col])

#write.csv(all.preds, "Z:/modelling/yxtable_scaled.csv")

summary(all.preds)
#transform 
all.preds$HT <- as.character(all.preds$HT)
all.preds$HT <- ifelse(is.na(all.preds$HT), "none", all.preds$HT)
all.preds <- na.omit(all.preds)
transformed <- all.preds

bc_values <- NULL
for(col in names(all.preds)[which(!names(all.preds) %in% c("ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
  cov <- all.preds[,col]
  cov <- cov + abs(min(cov)) + 1e-64
  par(mfrow = c(1,1))
  #hist(cov, main = col)
  fm <- lm(cov ~ 1)
  
  lambdas <- c(0, tan(seq(-1.2,1.2, length = 200))) # length can be set even higher, increase range if bc_values are at max/min values
  bc <- boxcox(fm, lambda = lambdas, plotit = F)
  bc_value <- bc$x[which.max(bc$y)]
  if(bc_value == 0){
    transformed[,col] <- log(cov)
  }else{
    transformed[,col] <- cov^bc_value
  }
  print(bc_value)
  bc_values <- append(bc_values, bc_value)
  print(bc_values)

  par(mfrow = c(2,2))
  qqnorm(transformed[,col], main = paste0(col, " transformed"))
  qqline(transformed[,col])
  hist(transformed[,col])
  
  qqnorm(all.preds[,col], main = paste0(col, " not transformed"))
  qqline(all.preds[,col])
  hist(all.preds[,col])
}

#save and check transforming parameters - plots should look the same

transform_sheet <- 
data.frame(name = names(all.preds)[which(!names(all.preds) %in% c("ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))],
           value = bc_values, shift = 0)

for(pred in names(all.preds)[which(!names(all.preds) %in% c("ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
  
  bt_value = transform_sheet$value[which(transform_sheet$name == pred)]
  if(bt_value == 0) bt <- exp(transformed[,pred])
  if(bt_value != 0) bt <- transformed[,pred]^(1/bt_value)
  
  shift <- median(orig - (bt + abs(min(all.preds[,pred])) + 1e-64))
  transform_sheet$shift[which(transform_sheet$name == pred)] <- shift
  hist(orig, main = paste0("original ", pred))
  hist((bt + abs(min(all.preds[,pred])) + 1e-64 + shift), main = paste0("backtransformed ", pred))
} 


write.csv(transform_sheet, "Z:/modelling/transform_sheet.csv")

#scale non-factor predictors
scale_transformed <- transformed

centers <- NULL
scales <- NULL
names <- NULL
for(col in which(!names(transformed) %in% c("ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))){
  scale_transformed[,col] <- scale(transformed[,col])
  centers <- append(centers, attr(scale_transformed[,col], "scaled:center"))
  scales <- append(scales, attr(scale_transformed[,col], "scaled:scale"))
  names <- append(names, col)
} 

#save and check scaling parameters - plots should look the same

scale_sheet <- data.frame("center" = centers, "scale" = scales, "name" = names(all.preds)[which(!names(all.preds) %in% c("ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))])

for(pred in names(all.preds)[which(!names(all.preds) %in% c("ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
  orig <- transformed[, pred]
  scaled <- scale_transformed[, pred]
  
  center_value = scale_sheet$center[which(scale_sheet$name == pred)]
  scale_value = scale_sheet$scale[which(scale_sheet$name == pred)]
  
  
  us <- scaled * scale_value + center_value
  
  hist(orig, main = paste0("original ", pred))
  hist(us, main = paste0("unscaled ", pred))
}

write.csv(scale_sheet, "Z:/modelling/scale_sheet.csv")

#check how the scaled, transformed values are distributed in comparison to the original values
par(mfrow = c(2,2))
for(i in which(!names(transformed) %in% c("ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))){
  hist(scale_transformed[,i], main = names(scale_transformed)[i], breaks = 1000)
  hist(all.preds[,i], main = names(scale_transformed)[i], breaks = 1000)
}


scale_transformed$CC <- as.factor(scale_transformed$CC)


source("Z:/GitRepo/modified_cormat.R")

#check for correlation
#correlation between values before and after transforamtion + scaling, HD old and HD new are highly NEGATIVELY correlated (just to keep it in mind)
rquery.cormat.spearman(cbind(all.preds[, which(!names(scale_transformed) %in% c("ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))],
                             scale_transformed[, which(!names(scale_transformed) %in% c("ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))]))$r


#correlation between transformed + scaled values
spearman <- rquery.cormat.spearman(scale_transformed[, which(!names(scale_transformed) %in% c("ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))], graphType = "heatmap")$r

spearman <- data.frame(apply(spearman, 1, function(x) as.numeric(as.character(x))))
cornames <- names(spearman) 
spearman$names <- cornames

cors <- na.omit(reshape2::melt(spearman))
cors[which(abs(cors$value) > 0.5 & cors$names != as.character(cors$variable)), ]


write.csv(scale_transformed, "Z:/modelling/yxtable_scaled_transformed.csv")

# QUICK AND DIRTY MODEL CHECKS
# scale_transformed_f <- scale_transformed[which(scale_transformed$HT != "bh"), ]
# 
# f.ranger <- COUNT ~ AI + CC + HD + LD + NB + PA + PI + SC + SL + SM + TC + TD + TV + VD + WA + Site + Country 
# f <- COUNT ~ AI + CC + HD + LD + NB + PA + PI + SC + SL + SM + TC + TD + TV + VD + WA + Site + SC:WA + SC:VD + TC:VD + TC:WA + WA^2 + PI^2 + Site:WA + Site:VD
# rF <- ranger::ranger(f.ranger, data = scale_transformed_f, importance = "impurity")
# plot(scale_transformed_f$COUNT ~ rF$predictions)
# plot(order(scale_transformed_f$COUNT) ~ order(rF$predictions))
# hist(scale_transformed_f$COUNT - rF$predictions)
# 
# fm <- glm(f, data = scale_transformed, family = "poisson")
# fm2 <- glm(f, data = scale_transformed_f, family = "poisson")
# 
# nb <- glm.nb(f, data = scale_transformed_f)
# nb1 <- glm.nb(f, data = scale_transformed)
# 
# nb.f <- summary(nb)
# nb.all <- summary(nb1)
# 
# nb.f$deviance/nb.f$null.deviance
# nb.all$deviance/nb.all$null.deviance
# 
# p.all <- summary(fm)
# p.f <- summary(fm2)
# 
# p.f$deviance/p.f$null.deviance
# p.all$deviance/p.all$null.deviance
# 
# plot(order(scale_transformed$COUNT) ~ order(nb1$fitted.values))
#      
# ds <- read.csv("Z:/climate_change/dry_season_2015_withgeom.csv")
# test <- merge(ds, all.preds[, c("ID", "SC")], by = "ID")
# test$all <- apply(test[,3:72], 1, function(x) sum(x, na.rm = T))
# plot(test$SC, test$all*1440)
# cor(test$SC, test$all)
# summary(test)
