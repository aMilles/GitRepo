library(spdep)
library(INLA)
data("Scotland")
Scotland$Region
setwd("Z:/NEMO_in")
load("input_ZWE.RData")
# f
# Y ~ mu.z + mu.y + HD + LD + NB + PA + PI + SC + SL + SM + 
#   TC + TD + VD + WA + Site
nb2INLA("dinla.txt", dinla)

ldat$ID <- 1:nrow(xy.adv)
f
m1 <- inla(f, data = ldat, control.predictor = list(compute = T), family = c("binomial", "zeroinflated.nbinomial.0"))

f
m1
summary(m1)
plot(m1)
plot(m1$summary.fitted.values$`0.5quant`[seq(length(n.z))], Y[seq(length(n.z)), 1], pch = "|", xlab = "predictions", ylab = "observations", main = "b")
plot(m1$summary.fitted.values$`0.5quant`[-seq(length(n.z))], Y[-seq(length(n.z)), 2], pch = "+", xlab = "predictions", ylab = "observations", main = "zinb")
cor(m1$summary.fitted.values$`0.5quant`[seq(length(n.z))], Y[seq(length(n.z)), 1], method = "spearman")^2
cor(m1$summary.fitted.values$`0.5quant`[-seq(length(n.z))], Y[-seq(length(n.z)), 2], method = "spearman")^2
plot(m1$summary.fitted.values$`0.5quant`[seq(length(n.z))], Y[seq(length(n.z)), 1], pch = "|")

m1$summary.fixed[,1]
m1$summary.fixed[sign(m1$summary.fixed$`0.025quant`) == sign(m1$summary.fixed$`0.975quant`), 1:2]
plot(m1)
save("output_ZWE.RData")

test <- glm(COUNT ~ HD + LD + NB + PA + PI + SC + SL + SM + 
              TC + TD + VD + WA + Site, data = xy, family = "poisson")
summary(test)
logLik(test)

hist(xy$SC)
WA <- read.csv("Z:/predictors/SC.csv")
TEST <- merge(WA, xy[, c("ID", "SC")], by = "ID")
plot(SC ~ mean, TEST)

WA <- read.csv("Z:/predictors/LD.csv")
TEST <- merge(WA, xy[, c("ID", "LD")], by = "ID")
plot(LD ~ mean, TEST)
