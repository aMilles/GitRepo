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

#read Model output from nemo
#load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_5km_splines.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LSO_5km.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_5km_splines.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LOSO_5km_splines.RData")
#load("Z:/NEMO_out/output_ZWE_complex_binomial_nonspatial_spatial_5km_splines.RData")
load("Z:/NEMO_out/output_BW")


if(is.numeric(xy$ID)) xy$ID <- segments_new$ID

#select one model
if("spatial_model" %in% ls()) that.model = spatial_model
if(!"spatial_model" %in% ls()) that.model = get(ls(pattern = "_spatial")[1])

#remove everything that is not needed
rm(list = ls()[which(!ls() %in% c("that.model", "xy", "splines", "n.knots", "f", "output_name", "segments"))])
#read file with GEC functions
source("Z:/GitRepo/function_file.R")
f <- stringi::stri_replace_all_fixed(f, "\n", "")
#get a vector of predictors used for the model
preds <- stringi::stri_remove_empty(stringi::stri_replace_all_fixed(strsplit(f, "+ ")[[1]], "+", ""))

#read the transformation and scale sheet
tf_sheet <- read.csv("Z:/modelling/transform_sheet.csv", stringsAsFactors = F)
scale_sheet <- read.csv("Z:/modelling/scale_sheet.csv", stringsAsFactors = F)
#"NA" predictor gets confused with NA ...
tf_sheet$name[is.na(tf_sheet$name)] <- "NA." -> scale_sheet$name[is.na(scale_sheet$name)]

#read dataset without scale/transformation and read data with scale and transformation
unscaled_df <- read.csv("Z:/modelling/yxtable.csv")
scaled_transformed <- read.csv("Z:/modelling/yxtable_scaled_transformed.csv")
#get the matching values of the full unscaled data frame (needed if a subset of the dataset was modelled)
if(nrow(unscaled_df) != nrow(xy)) unscaled_df <- unscaled_df[as.character(unscaled_df$ID) %in% as.character(xy$ID),]

{
  #calulcate the minima for transform shift
  mins <- apply(unscaled_df[,tf_sheet$name], 2, min)
  #sample the 0 - 100 % quantile of the predictor ranges with length 100 
  original.seqs <- seqs <- data.frame(apply(unscaled_df[,tf_sheet$name], 2, quantile, probs = seq(0,1,length = 100)))
  
  #original.seqs <- seqs <- data.frame(apply(rbind(mins, maxs), 2, function(x) seq(x[1], x[2], length = 100)))
  
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
  
  for(interaction in preds[as.logical(stringi::stri_count_fixed(preds, ":"))]){
    seqs[, interaction] <- model.matrix(as.formula(paste0("~", interaction)), seqs)[,2]
  }
  seqs <- cbind("(Intercept)" = 1, seqs)
}

 final.preds <- unlist(strsplit(preds, ":"))

#all possible combinations of factors used by the model - currently not used, factor effect is set to zero  
missing_factors <- unique(xy[,preds[!preds %in% names(seqs)]])
  

source("Z:/GitRepo/function_file.R")

conditional.plot(cond.vars = seqs[, names(seqs) %in% final.preds], effect.of = "WA", model = that.model, posterior.samples = 100, original.df = unscaled_df, rug = F, with.lines = T, interactions = preds[as.logical(stringi::stri_count_fixed(preds, ":"))])

#GENERATE EFFECT PLOTS
for(pred in unique(substring(names(seqs[, names(seqs) %in% preds]),1,2))){
  if(pred == "NA") pred <- "NA."
  gg <- conditional.plot(cond.vars = seqs[, names(seqs) %in% c("(Intercept)", final.preds)], effect.of = pred, model = that.model, posterior.samples = 10000, original.df = unscaled_df, rug = T, with.lines = F, interactions = preds[as.logical(stringi::stri_count_fixed(preds, ":"))]) +
    scale_y_log10()

  if(pred == "LD") gg <- gg + xlim(c(0,0.3))
  assign(pred, gg)
}


#SAVE THE EFFECT PLOTS
{
  pdf(file = paste0("Z:/Plots/effect_plots/", tools::file_path_sans_ext(output_name), ".pdf"), onefile = T)
  for(pred in unique(substring(names(seqs[, names(seqs) %in% preds]),1,2))){
    if(pred == "NA") pred <- "NA."
    grid.arrange(get(pred))
  } 
  dev.off()
}
