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
summary(xy)
segments <- readOGR("Z:/GEC/segments.shp")

segments$ID <- as.character(segments$ID)

xy <- xy[match(segments$ID, xy$ID), ]
xy <- xy[as.character(xy$ID) %in% as.character(segments$ID[stringi::stri_count_regex(as.character(segments$SC), "ZWE_") > 0]),]
specific_seg <- segments[stringi::stri_count_regex(as.character(segments$SC), "ZWE_") > 0,]
#xy$COUNT <- ifelse(xy$COUNT == 0, 0, log10(xy$COUNT))
xy$CC <- as.factor(xy$CC)

RMSE <- function(preds, obs){
  return(mean(sqrt((preds-obs)^2)))
}

# xy$HT <- as.character(xy$HT)
# xy$HT[is.na(xy$HT)] <- sample(c("mh", "bh"), size = sum(is.na(xy$HT)), replace = T)
#xy$HT <- as.factor(xy$HT)

f.xy <- na.omit(xy[xy$HT != "bh" | is.na(xy$HT),c(-1,-6)])
m.xy <- na.omit(xy[xy$HT != "mh" | is.na(xy$HT),c(-1,-6)])
xy <- na.omit(xy[,c(-1, -6)])

f.max <- FReibier::formulaMaker(xy[,-c(4,7)], y.col = 2)
f.min <- FReibier::formulaMaker(xy[,-c(4,7)], y.col = 2, interactions = F, quadratic = F)
f.optim <- formula(COUNT ~ AI + HD + LD + NB + PA + PI + SC + SL + SM + TC + TD + TV + VD + WA + I(WA^2) + I(PI^2) + SC:WA + SC:WA + TC:VD + TC:WA + CC)

p.glm <- glm(f.optim, data = xy, family = "poisson")
summary(p.glm)
#plot(p.glm)

p.glm <- glm(f.optim, data = f.xy, family = "poisson")
summary(p.glm)
p.glm <- glm(f.optim, data = m.xy, family = "poisson")
summary(p.glm)

gam(f.optim, data = xy, family = "poisson")

nb.glm <- MASS::glm.nb(f.optim, data = xy)
?MASS::glm.nb
summary(nb.glm)
nb.glm.f <- MASS::glm.nb(f.optim, data = f.xy)
summary(nb.glm.f)
nb.glm.m <- MASS::glm.nb(f.optim, data = m.xy)
summary(nb.glm.m)

cbind(summary(nb.glm.m)$coefficients[,1], summary(nb.glm.f)$coefficients[,1])


RMSE(nb.glm$fitted.values, xy$COUNT)
plot(order(nb.glm$fitted.values), order(xy$COUNT))
cor(order(nb.glm$fitted.values), order(xy$COUNT))^2

RMSE(nb.glm.f$fitted.values, f.xy$COUNT)
plot(order(nb.glm.f$fitted.values), order(f.xy$COUNT))
cor(order(nb.glm.f$fitted.values), order(f.xy$COUNT))^2

RMSE(nb.glm.m$fitted.values, m.xy$COUNT)
plot(order(nb.glm.m$fitted.values), order(m.xy$COUNT), pch = ".")
cor(order(nb.glm.m$fitted.values), order(m.xy$COUNT))^2
hist(residuals(nb.glm))

plot(allEffects(nb.glm))
plot(allEffects(nb.glm.f))
plot(allEffects(nb.glm.m))

effects::effectsTheme()
?allEffects
effect(SC:WA, nb.glm.m)
effect(nb.glm.m)
effects(nb.glm.m)

rF <- ranger(f.min, data = f.xy, importance = "impurity")
RMSE(predictions(rF), f.xy$COUNT)
names(xy)
rF <- randomForest(COUNT ~ AI + HD + LD + NB + PA + PI + SC + SM +SC + TC + TV + VD + WA + I(WA^2) + I(PI^2) + SC:WA + SC:WA + TC:VD + TC:WA, data = f.xy)

ranger(f.optim, data = xy)
names(xy)
rF


xy.binom <- read.csv("Z:/modelling/yxtable.csv")[,-c(1,2)]
xy.binom <- xy
xy.binom$COUNT <- as.numeric(xy.binom$COUNT > 0)

p.glm <- glm(f.optim, data = xy.binom, family = "binomial")
summary(p.glm)

plot(allEffects(p.glm))



test.bayes <- inla(COUNT ~ AI + HD + LD + NB + PA + PI + SC + SM +SC + TC + TV + VD + WA + I(WA^2) + I(PI^2) + SC:WA + SC:WA + TC:VD + TC:WA, data = xy)

?inla
