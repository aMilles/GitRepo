library(FReibier)
library(randomForest)
library(effects)
library(glmnet)
library(rgdal)
library(spdep)
library(INLA)
library(psych)

xy <- read.csv("Z:/modelling/yxtable.csv")[,-c(1)]
xy$ID <- as.character(xy$ID)
summary(xy)
segments <- readOGR("Z:/GEC/segments.shp")

segments$ID <- as.character(segments$ID)

xy <- xy[match(segments$ID, xy$ID), ]

xy <- xy[as.character(xy$ID) %in% as.character(segments$ID[stringi::stri_count_regex(as.character(segments$SC), "ZWE_") > 0]),]
specific_seg <- segments[stringi::stri_count_regex(as.character(segments$SC), "ZWE_") > 0,]
for(i in (2:ncol(xy))[-c(2,4)]) xy[,i] <- scale(xy[, i])

metric_segs <- spTransform(specific_seg, "+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")

dinla <- dnearneigh(as.matrix(geosphere::centroid(metric_segs[,])), d1 = 0.01, d2 = 6000)
nb2INLA("dinla.txt", dinla)

# test <- poly2nb(specific_seg) #full transects leads to cannot allocate vecotr of size 17,3 GB
# test2 <- nb2mat(test, style = "C", zero.policy = T)
# test2 <- as(test2, "dgTMatrix")


xy$ID <- 1:nrow(xy)
#testdata$idx <- 1:length.run
xy_m <- xy[, ]
n = nrow(xy_m)
E = sample(c(1,5,10,15), size=n, replace=TRUE)


form <- COUNT ~ AI + HD + LD + NB + PA + PI + SC + SM + TC + TV + VD + WA + I(WA^2) + I(PI^2) + SC:WA + SC:WA + TC:VD + TC:WA + f(ID, model = "besag", graph = "dinla.txt", hyper = list(prec = list(prior = "loggamma", param = c(0.1, 1), initial = 0.01)))
m1 <- inla(form, data = xy_m, control.predictor = list(compute = T), family = "zeroinflatedpoisson0", E = E, num.threads = 20)


summary(m1)
m1$cpu.used
