library(INLA)
library(rgdal)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(spdep)
library(MBA)
library(DescTools)
library(ggthemes)
library(rasterVis)
library(raster)
library(grid)
library(scales)
library(viridis)
library(gridExtra)


rm(list = ls()[which(ls() != "segments")])
#load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_5km_splines.RData")
load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LSO_5km.RData")

preds <- stringi::stri_remove_empty(stringi::stri_replace_all_fixed(strsplit(f, "+ ")[[1]], "+", ""))

invlogit <- function(x) exp(x)/(1 + exp(x))

sample.fixed.params <- function(model = model_test, nsample = 10000, param = "") {
  params <- names(model$marginals.fixed)
  which.param <- grep(param, params)
  
  if(sum(which.param)==0) stop("Check parameter names")
  dist <- matrix(ncol = length(which.param), nrow = nsample)
  for(i in 1:NROW(which.param)) {
    dist[,i] <- sample(model$marginals.fixed[[which.param[i]]][,1], nsample, TRUE,
                       model$marginals.fixed[[which.param[i]]][,2])
    
  } 
  dist <- data.frame(dist)
  names(dist) <- names(model$marginals.fixed)
  return(dist)
}

dist <- sample.fixed.params(model = xy_without_block10_spatial_model, nsample = 500)

#elaborated effect plots  

tf_sheet <- read.csv("Z:/modelling/transform_sheet.csv", stringsAsFactors = F)
scale_sheet <- read.csv("Z:/modelling/scale_sheet.csv", stringsAsFactors = F)
tf_sheet$name[is.na(tf_sheet$name)] <- "NA." -> scale_sheet$name[is.na(scale_sheet$name)]
unscaled_df <- read.csv("Z:/modelling/yxtable.csv")
unscaled_df <- unscaled_df[as.character(unscaled_df$ID) %in% as.character(xy$ID),]

maxs <- apply(unscaled_df[,tf_sheet$name], 2, max)
#maxs <- apply(unscaled_df[,tf_sheet$name], 2, quantile, probs = c(0.9))
#mins <- apply(unscaled_df[,tf_sheet$name], 2, quantile, probs = c(0.025))
mins <- apply(unscaled_df[,tf_sheet$name], 2, min)
original.seqs <- seqs <- data.frame(apply(rbind(mins, maxs), 2, function(x) seq(x[1], x[2], length = 100)))
  
for(i in seq(length(mins))){
  bt_value <- tf_sheet$value[i]
  scale <- scale_sheet$scale[i]
  center <- scale_sheet$center[i]
  name <- scale_sheet$name[i]
  cov <- seqs[,i]
  cov_shift <- cov + abs(mins[i]) + 1e-64
    
    
  if(bt_value == 0) tf_cov <- log(cov_shift)
  if(bt_value != 0) tf_cov <- cov_shift^bt_value  
    
  seqs[,i] <- (tf_cov - center)/scale
}
  
missing_factors <- unique(xy[,preds[!preds %in% names(seqs)]])
  

conditional.plot <- function(cond.vars = seqs[, names(seqs) %in% preds], effect.of = "VD", model = spatial_model, posterior.samples = 500, original.df, with.lines = F, rug = T){
  cond.vars[,-match(effect.of, names(cond.vars))] <- t(replicate(nrow(cond.vars),apply(cond.vars[,-match(effect.of, names(cond.vars))], 2, median), simplify = T))
  
  dist <- sample.fixed.params(model = model, nsample = posterior.samples)
  dist <- dist[,colnames(dist) %in% names(cond.vars)]
  dist <- dist[,match(names(dist), names(cond.vars))]
  
  m.cond.vars <- as.matrix(cond.vars)
  m.dist <- as.matrix(dist)
  
  out <- vector("list", length = nrow(m.dist))
  for(i in seq(nrow(m.dist))) out[[i]] <- as.vector(invlogit(m.cond.vars %*% m.dist[i,]))
  if(with.lines) gg.out.all <- reshape2::melt(data.frame(do.call(cbind, out), pred = original.seqs[,effect.of]), id.vars = "pred")
  gg.out.q <- data.frame(t(apply(do.call(cbind, out), 1, quantile, probs = c(0.025, 0.5, 0.975))), pred = original.seqs[,effect.of])
  names(gg.out.q)[1:3] <- c("lower", "mid", "upper")
  original.df <- data.frame(pred = original.df[,effect.of], zero = 0)
  
  gg <- 
    ggplot()+
      geom_ribbon(data = gg.out.q, aes(x = pred, ymin = lower*100, ymax = upper*100), fill = "gray")+
      geom_line(data = gg.out.q, aes(x = pred, y = mid*100), size = 1.5)+
      xlab(effect.of)+
      ylab("propability [%]")+
      theme_bw() 
  
  if(with.lines) gg <- gg + geom_line(data = gg.out.all, aes(x = pred, y = value * 100, group = variable), col = "red", alpha = .1)
  if(rug) gg <- gg + geom_point(data = original.df, aes(x = pred, y = zero), shape = "|", alpha = .1)
  
  if(!rug){
    gg <- grid.arrange(
      gg + theme(plot.margin=unit(c(-0.08,1,1,1), "cm")), 
      ggplot(original.df, aes(x = pred))+
      geom_histogram(bins = 100, col = "black", fill = "gray")+
      theme_void()+
      ggtitle(paste0("Effect of: ", effect.of))+
      theme(plot.margin=unit(c(.3,1,-0.08,1), "cm")),
      layout_matrix = matrix(c(2,rep(1,7),2,rep(1,7)), nc = 2))
  }
  
  
  return(gg)
}

conditional.plot(effect.of = "VD", model = xy_without_block10_spatial_model, posterior.samples = 50, original.df = unscaled_df, rug = F, with.lines = T)
conditional.plot(effect.of = "VD", model = xy_without_block10_spatial_model, posterior.samples = 50, original.df = unscaled_df)

for(pred in names(seqs[, names(seqs) %in% preds])) assign(pred, conditional.plot(effect.of = pred, model = xy_without_block10_spatial_model,posterior.samples = 50, original.df = unscaled_df, rug = T))

do.call(grid.arrange, mget(names(seqs[, names(seqs) %in% preds])))
