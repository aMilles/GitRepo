library(MASS)
library(ranger)
library(INLA)
library(spdep)
library(rgdal)

xy <- read.csv("Z:/modelling/yxtable.csv")
xy$PA[is.na(xy$PA)] <- 0
xy$PA <- as.factor(xy$PA)

poly <- readOGR('Z:/GEC/segments.shp')
poly <-poly[poly$ID %in% xy$ID,]
w <- poly2nb(poly, row.names=poly$ID)
wm<- nb2mat(w, style='B', zero.policy = T)
rwm <- mat2listw(wm, style='W')
moran.plot(poly$obs_cnt, rwm)
ww <- nb2listw(w, style='B',zero.policy = T)
moran(poly$obs_cnt, ww,n=length(ww$neighbours), S0=Szero(ww), zero.policy = T)
moran.test(poly$obs_cnt, ww, randomisation = F, zero.policy = T)
moran.mc(poly$obs_cnt, ww, nsim=99, zero.policy = T)

n <- length(poly)
ms <- cbind(id=rep(1:n, each=n), y=rep(y, each=n), value=as.vector(wm * y))
ms <- ms[ms[,3] > 0, ]
ams <- aggregate(ms[,2:3], list(ms[,1]), FUN=mean)
ams <- ams[,-1]
colnames(ams) <- c('y', 'spatially lagged y')
plot(ams)
reg <- lm(ams[,2] ~ ams[,1])
abline(reg, lwd=2)
abline(h=mean(ams[,2]), lt=2)
abline(v=ybar, lt=2)

?nb2mat
xy <- na.omit(xy)
summary(xy)
xy_binom <- xy

xy_binom$y <- as.numeric(xy_binom$y > 0)


test <- matrix(1.23, 48174, 48174)
rm(test)
nb.xy <- glm.nb(y ~ AI + HD + LD + PA + PI + SC +TC + VD + WA, data = xy)
summary(nb.xy)

nb.xy <- glm.nb(y ~ AI + HD + LD + PA + PI + SC +TC + VD + WA + (PI^2) + (WA)^2 + (SC:WA) + (TC:WA) + (TC:VD) , data = xy)

summary(nb.xy)

ranger.xy <- ranger(y ~ AI + HD + LD + PA + PI + SC +TC + VD + WA + (PI^2) + (WA)^2 + (SC:WA) + (TC:WA) + (TC:VD) , data = xy)
ranger.xy
ranger::importance(ranger.xy)
ranger::importance_pvalues(ranger.xy)
rF.xy <- ranger::holdoutRF(y ~ AI + HD + LD + PA + PI + SC +TC + VD + WA + (PI^2) + (WA)^2 + (SC:WA) + (TC:WA) + (TC:VD) , data = xy)
importance_pvalues(rF.xy)
?importance_pvalues

sqrt(mean((ranger.xy$predictions - xy$y)^2))
sqrt(mean((ranger.xy$predictions - xy$y)^2))
plot(ranger.xy$predictions, xy$y)

ranger.xy_binom <- ranger(y ~ AI + HD + LD + PA + PI + SC +TC + VD + WA + (PI^2) + (WA)^2 + (SC:WA) + (TC:WA) + (TC:VD) , data = xy, probability = T)
sqrt(mean((ranger.xy_binom$predictions - xy_binom$y)^2))
plot(ranger.xy_binom$predictions[,1], rep(xy_binom$y, 1))

ranger.xy