library(spdep)
library(INLA)

setwd("/home/fr/fr_fr/fr_am595/data")
xy <- read.csv("xy.csv")[,-c(1)]
segments <- rgdal::readOGR("segments.shp")
xy$ID <- as.character(xy$ID)

metric_segs <- sp::spTransform(segments, "+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")

dinla <- dnearneigh(as.matrix(geosphere::centroid(metric_segs[,])), d1 = 0.01, d2 = 6000)
nb2INLA("dinla.txt", dinla)

xy$ID <- 1:nrow(xy)

n = 100
E = sample(c(1,5,10,15), size=n, replace=TRUE)

form <- COUNT ~ AI + HD + LD + NB + PA + PI + SC + SM + TC + TV + VD + WA + I(WA^2) + I(PI^2) + SC:WA + SC:WA + TC:VD + TC:WA + f(ID, model = "besag", graph = "dinla.txt", hyper = list(prec = list(prior = "loggamma", param = c(0.1, 1), initial = 0.01)))
m1 <- inla(form, data = xy[1:100,], control.predictor = list(compute = T), family = "zeroinflatedpoisson1", E = E)

save.image("ouput.RData")
