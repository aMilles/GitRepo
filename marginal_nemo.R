library(INLA)
library(gridExtra)
library(parallel)
library(foreach)

load("/home/fr/fr_fr/fr_am595/data/run_this/output_all_complex_binomial_nonspatial_spatial_5km_splines.RData")
setwd("/home/fr/fr_fr/fr_am595/data/bootstrap/")

if(is.numeric(xy$ID)) xy$ID <- segments_new$ID

#select one model
if("spatial_model" %in% ls()) that.model = spatial_model
if(!"spatial_model" %in% ls()) that.model = get(ls(pattern = "_spatial")[1])

#remove everything that is not needed
rm(list = ls()[which(!ls() %in% c("that.model", "xy", "splines", "n.knots", "f", "output_name", "segments"))])
#read file with GEC functions
source("https://raw.githubusercontent.com/aMilles/GitRepo/master/function_file.R")
f <- stringi::stri_replace_all_fixed(f, "\n", "")
#get a vector of predictors used for the model
preds <- stringi::stri_remove_empty(stringi::stri_replace_all_fixed(strsplit(f, "+ ")[[1]], "+", ""))

#read the transformation and scale sheet
tf_sheet <- read.csv("transform_sheet.csv", stringsAsFactors = F)
scale_sheet <- read.csv("scale_sheet.csv", stringsAsFactors = F)
#"NA" predictor gets confused with NA ...
tf_sheet$name[is.na(tf_sheet$name)] <- "NA." -> scale_sheet$name[is.na(scale_sheet$name)]

#read dataset without scale/transformation and read data with scale and transformation
unscaled_df <- read.csv("yxtable.csv")
scaled_transformed <- read.csv("yxtable_scaled_transformed.csv")
#get the matching values of the full unscaled data frame (needed if a subset of the dataset was modelled)
if(nrow(unscaled_df) != nrow(xy)) unscaled_df <- unscaled_df[as.character(unscaled_df$ID) %in% as.character(xy$ID),]

{
  #calulcate the minima for transform shift
  mins <- apply(unscaled_df[,tf_sheet$name], 2, min)
  maxs <- apply(unscaled_df[,tf_sheet$name], 2, max)
  #sample the 0 - 100 % quantile of the predictor ranges with length 100 
  #original.seqs <- seqs <- data.frame(apply(unscaled_df[,tf_sheet$name], 2, quantile, probs = seq(0,1,length = 100)))
  
  original.seqs <- seqs <- data.frame(apply(rbind(mins, maxs), 2, function(x) seq(x[1], x[2], length = 100)))
  
  #transform and scale the quantiles of the unscaled dataset
  for(i in seq(nrow(tf_sheet))){
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
  
  #use the nearest scaled+transformed value for each scaled+transformed quantile to get an estimate of the splines
  for(spline in splines){
    #this.spline <- data.frame((spliner(c(scaled_transformed[,spline], seqs[,spline]), n.knots)[seq(100)+nrow(scaled_transformed), ]))
    close <- sapply(seqs[,spline], function(x) which.min(abs(x - xy[,spline])))
    this.spline <- xy[close,paste0(spline, seq(n.knots))]
    names(this.spline) <- paste0(spline, seq(n.knots))
    seqs <- data.frame(seqs, this.spline)
  }  
  
  final.preds <- unique(unlist(strsplit(preds, ":")))
  missing.preds <- unique(final.preds)[!unique(final.preds) %in% names(seqs)]
  seqs[,missing.preds] <- unique(xy[,missing.preds])[12,]
  
  interactions <- preds[as.logical(stringi::stri_count_fixed(preds, ":"))]
  
  for(interaction in preds[as.logical(stringi::stri_count_fixed(preds, ":"))]){
    seqs[, interaction] <- model.matrix(as.formula(paste0("~", interaction)), seqs)[,2]
  }
  seqs <- cbind("(Intercept)" = 1, seqs)
}


#GENERATE EFFECT PLOTS
selection <- unique(substring(names(seqs)[names(seqs) %in% preds & !names(seqs) %in% missing.preds & !names(seqs) %in% interactions],1,2))
selection[selection == "NA"] <- "NA."

samples <- sample(nrow(unscaled_df), nrow(unscaled_df), replace = F)

cl <- doParallel::makeCluster(20)
doParallel::registerDoParallel(cl)

out <- foreach(i = selection, .combine = "cbind")%dopar%{
  assign(
    predictor, 
    marginal.plot(
      marg.vars = seqs[, names(seqs) %in% c(final.preds)], 
      effect.of = i, 
      model = that.model, 
      posterior.samples = 10, 
      original.df = unscaled_df[1:100,],
      scaled.df = xy[1:100,],
      interactions = interactions, 
      factor.values = missing.preds,
      prob = c(0.025, 0.25, 0.5, 0.75, 0.975),
      original.seqs = original.seqs,
      n.cores = 1,
      silent = T))
} 

doParallel::stopCluster(cl)

save("out", file = paste0("marginals_", output_name))
