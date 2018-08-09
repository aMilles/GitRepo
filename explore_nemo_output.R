library(INLA)
library(inlabru)
load("Z:/NEMO_out/output_all_poisson.RData")

rm(m1)
load("Z:/NEMO_out/output_all_pardiso_parallel_nb.RData")
nb <- m1
plot(nb)
summary(nb)
rm(nb)
rm(zip0)
load("Z:/NEMO_out/output_all_pardiso_parallel_zip0.RData")
zip0 <- m1
plot(zip0)
summary(zip0)
load("Z:/NEMO_out/output_all_pardiso_parallel_zip1.RData")

load("Z:/NEMO_out/output_none.RData")
zip1 <- m1
plot(zip1)
summary(zip1)
hist(xy$TV)
c <- xy$COUNT
p_nb <- nb$summary.fitted.values$`0.5quant`
p_nb[p_nb > max(c)] <- max(c)
p_zip0 <- zip0$summary.fitted.values$`0.5quant`
p_zip0[p_zip0 > max(c)] <- max(c)
p_zip1 <- zip1$summary.fitted.values$`0.5quant`
p_zip1[p_zip1 > max(c)] <- max(c)


par(mfrow = c(2,1))
plot(c~p_nb, xlim = c(0, max(c)), xlab = "predicted count", ylab = "observed count", main = "negative binomial")
plot(c~p_zip0, xlim = c(0, max(c)), xlab = "predicted count", ylab = "observed count", main = "ZIP")
#plot(c~p_zip1, xlim = c(0, max(c)))

hist(p_nb - c, xlim = c(-50, 50), breaks = 10000, main = "negative binomial", xlab = "residuals", col = "gray")
hist(p_zip0 - c, xlim = c(-50, 50), breaks = 10000, main = "ZIP", xlab = "residuals", col = "gray")
#hist(p_zip1 - c, xlim = c(-50, 50), breaks = 1000)
sum(p_nb - c)
sum(p_zip0 - c)
#sum(p_zip1 - c)

sum(p_nb)
sum(p_zip0)
sum(p_zip1)
sum(c)

t <- inla.posterior.sample(result = m1)
?inla.posterior.sample
t
a <- inla.posterior.sample
test <- poisson_m$model.matrix
test@Dimnames
nb_output <- m1$summary.fitted.values$mean
poisson_output <- poisson_m$summary.fitted.values$mean
plot(xy$COUNT~poisson_output, xlim = c(0,500))
plot(xy$COUNT~m1$summary.fitted.values$`0.5quant`, xlim = c(0,300), ylim = c(0,300))
plot(xy$COUNT~poisson_output, xlim = c(0,500))
cor(xy$COUNT[xy$COUNT>0],poisson_output[xy$COUNT>0])
cor(xy$COUNT[xy$COUNT>0],nb_output[xy$COUNT>0])

plot(order(m1$summary.fitted.values$`0.5quant`) ~ order(xy$COUNT))
cor(order(m1$summary.fitted.values$`0.5quant`), order(xy$COUNT))^2
cor(order(m1$summary.fitted.values$`0.5quant`), order(xy$COUNT))
m1 <- nb
plot(x = seq(length((m1$summary.fitted.values$`0.5quant`))), y = (m1$summary.fitted.values$`0.5quant`)[order(m1$summary.fitted.values$`0.5quant`)], type = "l", ylim = c(0, 40), cex = 3, xlab = "prediction rank", ylab = "observed and predicted value")
points(x = seq(length((m1$summary.fitted.values$`0.5quant`))), (m1$summary.fitted.values$`0.025quant`)[order(m1$summary.fitted.values$`0.5quant`)], col = "lightgray")
points(x = seq(length((m1$summary.fitted.values$`0.5quant`))), (m1$summary.fitted.values$`0.975quant`)[order(m1$summary.fitted.values$`0.5quant`)], col = "lightgray")
points(x = seq(length((m1$summary.fitted.values$`0.5quant`))), xy$COUNT[order(m1$summary.fitted.values$`0.5quant`)], col = "red", pch = "|", cex = .5)

summary(nb_output)
