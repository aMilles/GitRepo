library(INLA)
library(inlabru)
size = 50
y = 1:size + rnorm(size)
x = (1:size)/10

y[sample(1:size, size/10)] <- NA 

yx <- data.frame(y = y, x = x)
res <- bru(y ~ x, data = yx, family = "gaussian")
plot(res$summary.fitted.values$mode)
yx 
res$formula = y~x
xpost = predict(res, NULL, ~x)
plot(res)
plot(y ~ xpost$median)
predict(bru)

