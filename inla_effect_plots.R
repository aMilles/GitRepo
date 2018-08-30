library(INLA)
library(ggplot2)
library(gridExtra)
library(parallel)
library(foreach)
library(doSNOW)
library(doParallel)
library(gridExtra)
library(ranger)

#read Model output from nemo
load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_5km_prec1e-04.RData")

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
f <- stringi::stri_replace_all_regex(f, " ", "")
f <- stringi::stri_replace_all_fixed(f, "+", " + ")
f <- stringi::stri_replace_all_fixed(f, " + f(Site,model=\"iid\")", "")
preds <- as.vector(na.omit(stringi::stri_replace_all_fixed(strsplit(f, "+ ")[[1]], "+", NA)))

#read the transformation and scale sheet
tf_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/transform_sheet.csv", stringsAsFactors = F, row.names = 1)
scale_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/scale_sheet.csv", stringsAsFactors = F)
#"NA" predictor gets confused with NA ...
tf_sheet$name[is.na(tf_sheet$name)] <- "NA." -> scale_sheet$name[is.na(scale_sheet$name)]

#read dataset without scale/transformation and read data with scale and transformation
unscaled_df <- read.csv("C:/Users/amilles/Dropbox/modelling/yxtable.csv")
all_unscaled <- unscaled_df
scaled_transformed <- read.csv("C:/Users/amilles/Dropbox/modelling/yxtable_scaled_transformed.csv")
transformed <- read.csv("C:/Users/amilles/Dropbox/modelling/yxtable_transformed.csv")
#get the matching values of the full unscaled data frame (needed if a subset of the dataset was modelled)
if(nrow(unscaled_df) != nrow(xy)) unscaled_df <- unscaled_df[as.character(unscaled_df$ID) %in% as.character(xy$ID),]

{
  #calulcate the minima for transform shift
  all.mins <- apply(all_unscaled[,tf_sheet$name], 2, min)
  mins <- apply(unscaled_df[,tf_sheet$name], 2, min)
  maxs <- apply(unscaled_df[,tf_sheet$name], 2, max)
  #sample the 0 - 100 % quantile of the predictor ranges with length 100 
  original.seqs <- seqs <- data.frame(apply(unscaled_df[,tf_sheet$name], 2, quantile, probs = seq(0,1,length = 100)))
  check.df <- all_unscaled
  #original.seqs <- seqs <- data.frame(apply(rbind(mins, maxs), 2, function(x) seq(x[1], x[2], length = 100)))
  
  #transform and scale the quantiles of the unscaled dataset
  for(i in seq(ncol(seqs))){
    name <- names(seqs)[i]
    bt_value <- tf_sheet$value[tf_sheet$name == name]
    shift <- tf_sheet$shift[tf_sheet$name == name]
    
    scale <- scale_sheet$scale[scale_sheet$name == name]
    center <- scale_sheet$center[scale_sheet$name == name]
    
    min <- all.mins[i]
    cov <- seqs[,i]
    check <- all_unscaled[,name]

    cov <- cov - shift
    check <- check - shift
    
    if(bt_value == 0) {tf_cov <- log(cov); tf_check <- log(check)}
    if(bt_value != 0) {tf_cov <- cov^bt_value; tf_check <- check^bt_value}
    
    seqs[,i] <- (tf_cov - center)/scale
    check.df[,name] <- (tf_check - center)/scale
  }
  if(!all(round(apply(check.df[,names(seqs)], 2, sd), 5) == 1)) warning("scale or transformation maybe did not work")
  
  #use the nearest scaled+transformed value for each scaled+transformed quantile to get an estimate of the splines
  if(!all(is.na(splines))){
    for(spline in splines){
      #this.spline <- data.frame((spliner(c(scaled_transformed[,spline], seqs[,spline]), n.knots)[seq(100)+nrow(scaled_transformed), ]))
      close <- sapply(seqs[,spline], function(x) which.min(abs(x - xy[,spline])))
      this.spline <- xy[close,paste0(spline, seq(n.knots))]
      names(this.spline) <- paste0(spline, seq(n.knots))
      seqs <- data.frame(seqs, this.spline)
    }   
  }
  
  final.preds <- unique(unlist(strsplit(preds, ":")))
  missing.preds <- unique(final.preds)[!unique(final.preds) %in% names(seqs)]
  seqs[,missing.preds] <- unique(xy[,missing.preds])
  
  interactions <- preds[as.logical(stringi::stri_count_fixed(preds, ":"))]
  
  # for(interaction in preds[as.logical(stringi::stri_count_fixed(preds, ":"))]){
  #   seqs[, interaction] <- model.matrix(as.formula(paste0("~", interaction)), seqs)[,2]
  # }
  seqs <- cbind("(Intercept)" = 1, seqs)
}
apply(seqs, 2, sd)


for(factors in names(Filter(is.factor, seqs))){
  factor.levels <- levels(seqs[,factors])
  seqs[,factors] <- as.factor(sort(rep(factor.levels, length = nrow(seqs))))
  original.seqs[,factors] <- as.factor(sort(rep(factor.levels, length = nrow(seqs))))
}

#GENERATE EFFECT PLOTS
selection <- unique(substring(names(seqs)[names(seqs) %in% preds & !names(seqs) %in% missing.preds & !names(seqs) %in% interactions],1,2))
selection[selection == "NA"] <- "NA."


full.names <- data.frame(matrix(c("WA", "distance to water [m]",
                                  "AI", "proportion of agriculture",
                                  "HD", "human population per area",
                                  "LD", "tropical livestock units / km²",
                                  "NA.", "extractable sodium in soils",
                                  "NB", "slope [percent]",
                                  "PI", "accessibility from roads [s]",
                                  "SC", "4-month precipitation sum [mm]",
                                  "SL", "altitude [m]",
                                  "SM", "soil moisture",
                                  "TC", "temperature [K]",
                                  "TD", "tree cover [percent]",
                                  "TV", "terrain ruggedness",
                                  "CC", "Koppen-Geiger-Climate-Zone",
                                  "PA", "protected areas",
                                  "VD", "EVI"), byrow = T, nc = 2), stringsAsFactors = F)
names(full.names) <- c("pred", "full_name")

for(i in c("VD", "WA", "PI", "NA.", "AI", "SM", "LD")){
  source("Z:/GitRepo/function_file.R")
  dist <- sample.fixed.params(that.model, 10)
  
  check <- marginal.quantiles(
    marg.vars = seqs, 
    effect.of = i, 
    posterior.dist = dist, 
    original.df = unscaled_df,
    scaled.df = xy,
    interactions = interactions, 
    original.seqs = original.seqs,
    silent = F,
    formula = f,
    interact = T)
  
  
  out <- check
  df <- out[[3]]
  out[[5]] <- merge(out[[5]], full.names, by = "pred")
  out[[4]] <- full.names$full_name[which(full.names$pred %in% out[[4]])]
  names(df) <- "pred"
  comb.out <- cbind(out[[1]][,-4], out[[2]])
  names(comb.out) <- c("lower.in", "mid", "upper.in", "lower.out", "upper.out", "pred")
  
  comb.out <- data.frame(comb.out, 
                         interactor = rep(out[[5]]$full_name, each = 100), 
                         interactor.value = rep(paste0(out[[5]]$pred, "_", signif(as.numeric(out[[5]]$value), 3)), each = 100))
  
  
assign(paste0("comb.out_", i), comb.out)
}


for(i in ls(pattern = "comb.out_")){
  comb.out <- get(i)
  print(i)
  predictor <- strsplit(i, "_")[[1]][2]
  effect <- full.names[full.names[,1] == predictor,2]
  df <- data.frame(unscaled_df[, predictor])
  names(df) <- "pred"

  plot <- (ggplot(comb.out, aes(x = pred))+
             geom_ribbon(aes(ymin = lower.out, ymax = upper.out, fill = interactor.value), alpha = .6)+
             #geom_ribbon(aes(ymin = lower.in, ymax = upper.in), alpha = .8)+
             geom_line(aes(y = mid, col = interactor.value, group = interactor.value))+
             scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(comb.out$interactor.value)), "Paired"))+
             scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(comb.out$interactor.value)), "Paired"))+
             geom_point(data = df, aes(x = pred, y = 0), shape = "|", alpha = 0.05)+
             facet_wrap(~ interactor, scales = "free")+
             xlab(effect)+
             ylab("Habitat Suitability")+
             theme_classic()+
             guides(fill = guide_legend(title = "interaction value"), color = F, alpha = F, shape = F)+
             theme(text = element_text(size = 20)))
  
  png(paste0("Z:/Plots/effect_plots/marginals/Marginal_Effect_Interaction_", effect,"_", output_name, ".png"), height = 900, width = 500 * length(unique(comb.out$interactor)))
  print(plot)
  dev.off()
    
  # pdf(paste0("Z:/Plots/effect_plots/marginals/Marginal_Effect_Interaction_", effect,"_", output_name, ".pdf"), height = 6, width = 5 * length(unique(comb.out$interactor)))
  # print(plot)
  # dev.off()
    
}
  


scaled_transformed$obs <- as.factor(as.numeric(scaled_transformed$COUNT > 0))
f.ranger <- as.formula("obs ~ AI + LD + TV + PI + SM + VD + WA + NA.")
rF <- ranger(f.ranger, data = scaled_transformed, probability  = T)

check <- marginal.quantiles(
  marg.vars = seqs, 
  effect.of = "WA", 
  posterior.dist = dist, 
  original.df = unscaled_df,
  scaled.df = xy,
  interactions = interactions, 
  original.seqs = original.seqs,
  silent = F,
  formula = "AI + LD + TV + PI + SM + VD + WA + NA.",
  interact = F,randomForest = T, rF.model = rF)

out <- check
df <- out[[3]]
names(df) <- "pred"

print(effect.plot(quantiles.i = out[[1]],
                  quantiles.o = out[[2]],
                  lines = out[[3]],
                  rug = F, 
                  with.lines = F,
                  df = df,
                  effect.of = out[[4]]))





samples <- sample(nrow(unscaled_df), nrow(unscaled_df), replace = F)

xy$Site <- "ZWE_MAT"
check <- marginal.quantiles(
  marg.vars = seqs, 
  effect.of = "WA", 
  posterior.dist = dist, 
  original.df = unscaled_df,
  scaled.df = xy,
  interactions = interactions, 
  original.seqs = original.seqs,
  silent = F,
  formula = f,
  interact = F, randomForest = F)


out <- check
df <- out[[3]]
names(df) <- "pred"

png(paste0("Z:/Plots/effect_plots/marginals/Marginal_Effect_Interaction_", out[[4]],"_", output_name, ".png"), height = 500, width = 500 )
print(effect.plot(quantiles.i = out[[1]],
            quantiles.o = out[[2]],
            lines = out[[3]],
            rug = F, 
            with.lines = F,
            df = df,
            effect.of = out[[4]]))
dev.off()
          
                

#### SIMPLE ####
gg.xy <- xy
gg.xy[,selection] <- unscaled_df[,selection]
gg.xy <- reshape2::melt(data.frame(gg.xy[,c("obs", selection)], "pred" = that.model$summary.fitted.values$`0.5quant`, "site" = xy$Site), id.vars = c("obs", "pred", "site"))

{
  pdf(file = paste0("Z:/Plots/effect_plots/cheap_effect_", tools::file_path_sans_ext(output_name), ".pdf"), onefile = F)
  ggplot(gg.xy, aes(y = pred*100, x = value))+
    geom_point(alpha = 0.01)+
    geom_smooth()+
    #geom_point(aes(y = min(pred*100), x = value), shape = "|", alpha = 0.01)+
    facet_grid(site~variable, scales = "free")+
    ylim(values = c(0,100))+
    ylab("detection probability [%]")
  dev.off()
}

  
#### CONDITIONAL ####  
  
  conditional.plot(cond.vars = seqs[, names(seqs) %in% c("(Intercept)", final.preds)], effect.of = "VD", model = that.model, posterior.samples = 100, original.df = unscaled_df, rug = F, with.lines = T, interactions = preds[as.logical(stringi::stri_count_fixed(preds, ":"))], factor.values = missing.preds, fix.fun = function(x) quantile(x, c(0.8)))
  warnings()
  
  
  for(pred in selection){
    if(pred == "NA") pred <- "NA."
    gg <- conditional.plot(
      cond.vars = seqs[, names(seqs) %in% c("(Intercept)", final.preds)], 
      effect.of = pred, 
      model = that.model, 
      posterior.samples = 10000, 
      original.df = unscaled_df, 
      rug = T, 
      with.lines = F, 
      interactions = interactions, 
      factor.values = missing.preds,
      prob.inner = c(0.25, 0.5, 0.75),
      prob.outer = c(0.25, 0.5),
      fix.fun = function(x) quantile(x, c(0.8)))
    assign(pred, gg)
  }
  
  #SAVE THE EFFECT PLOTS
  {
    pdf(file = paste0("Z:/Plots/effect_plots/", tools::file_path_sans_ext(output_name), ".pdf"), onefile = T)
    for(pred in selection){
      grid.arrange(get(pred))
    } 
    dev.off()
  }
  
  source("Z:/GitRepo/function_file.R")
  