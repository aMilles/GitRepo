library(INLA)
library(ranger)
library(mgcv)
library(devtools)
library(FReibier)
library(randomForest)
library(effects)
library(glmnet)
library(rgdal)
library(spdep)
xy <- read.csv("Z:/modelling/yxtable.csv")[,-c(1)]
xy$ID <- as.character(xy$ID)

xy$TS <- unlist(lapply(strsplit(xy$ID, "p"), function(x) x[1]))
xy <- xy[, - c(1, 5)]
xy <- na.omit(xy)
sds <- aggregate(xy, by = list(xy$TS), FUN = sd)
summary(sds)
apply(xy, 2, sd)
apply(na.omit(sds[,-c(1, 17)]), 2, sd)

plot(sds$AI)
plot(sds$SC)
plot(sds$WA)
plot(sds$AI)
plot(sds$AI)


sds <- aggregate(xy, by = list(xy$TS), FUN = function(x) return(max(x) - min(x)))
apply(xy, 2, function(x) return(max(x) - min(x)))

?aggregate
