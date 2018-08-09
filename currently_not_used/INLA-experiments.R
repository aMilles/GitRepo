library(FReibier)
library(randomForest)
library(effects)
library(glmnet)
library(rgdal)
library(spdep)
library(INLA)
library(psych)

xy <- read.csv("Z:/modelling/yxtable.csv")[,-c(1)]
summary(xy)
#segments <- readOGR("Z:/GEC/segments.shp")
xy <- xy[as.character(xy$ID) %in% as.character(segments$ID[stringi::stri_count_regex(as.character(segments$SC), "ZWE_") > 0]),]
specific_seg <- segments[stringi::stri_count_regex(as.character(segments$SC), "ZWE_") > 0,]

pairs.panels(xy)


test <- poly2nb(specific_seg) #full transects leads to cannot allocate vecotr of size 17,3 GB
test2 <- nb2mat(test, style = "C", zero.policy = T)
test2 <- as(test2, "dgTMatrix")

xy$ID <- 1:nrow(xy)
testdata$idx <- 1:length.run

form <- COUNT ~ AI + HD + LD + NB + PA + PI + SC + SM +SC + TC + TV + VD + WA + I(WA^2) + I(PI^2) + SC:WA + SC:WA + TC:VD + TC:WA + f(ID, model = "besag", graph = test2)
m1 <- inla(form, data = xy, control.predictor = list(compute = T))
Sys.time()
INLA::

i = 1
df <- (m1$marginals.linear.predictor)
df <- lapply(df, inla.smarginal)
df[[1]]
m1$size.linear.predictor
plot(df[[7000]], type = "l")
summary(m1)
