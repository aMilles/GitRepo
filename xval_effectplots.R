library(INLA)
library(ggplot2)
library(gridExtra)
library(parallel)
library(foreach)
library(doSNOW)
library(doParallel)
library(gridExtra)
library(ranger)


load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_xval_LOSO_6km_prec1e-04.RData")
rm(list = ls(pattern = "nonspatial"))

f <- stringi::stri_replace_all_fixed(f, "\n", "")
#get a vector of predictors used for the model
f <- stringi::stri_replace_all_regex(f, " ", "")
f <- stringi::stri_replace_all_fixed(f, "+", " + ")
f <- stringi::stri_replace_all_fixed(f, "+ f(Site,model=\"iid\") + f(PA,model=\"iid\")", "")
f <- stringi::stri_replace_all_fixed(f, "+ Site + PA", "")
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
  missing.preds <- missing.preds[stringi::stri_count_fixed(missing.preds, "^2") < 1]
  seqs[,missing.preds] <- unique(xy[,missing.preds])
  
  interactions <- preds[as.logical(stringi::stri_count_fixed(preds, ":"))]
  
  # for(interaction in preds[as.logical(stringi::stri_count_fixed(preds, ":"))]){
  #   seqs[, interaction] <- model.matrix(as.formula(paste0("~", interaction)), seqs)[,2]
  # }
  seqs <- cbind("(Intercept)" = 1, seqs)
}
apply(seqs[,-22], 2, median)


for(factors in names(Filter(is.factor, seqs))){
  factor.levels <- levels(seqs[,factors])
  seqs[,factors] <- as.factor(sort(rep(factor.levels, length = nrow(seqs))))
  original.seqs[,factors] <- as.factor(sort(rep(factor.levels, length = nrow(seqs))))
}

#GENERATE EFFECT PLOTS
selection <- unique(substring(names(seqs)[names(seqs) %in% preds & !names(seqs) %in% missing.preds & !names(seqs) %in% interactions],1,2))
selection[selection == "NA"] <- "NA."


full.names <- data.frame(matrix(c("WA", "WA - distance to water [km]",
                                  "AD", "AD - agricultural density [%]",
                                  "HD", "HD - human density [people/km²]",
                                  "LD", "LD - livestock density [TLU/km²]",
                                  "SS", "SS - sodium in soils [cmolc / kg]",
                                  "NB", "NB - slope [%]",
                                  "AR", "AR - accessibility from roads [h]",
                                  "SC", "4-month precipitation sum [mm]",
                                  "SL", "altitude [m]",
                                  "SM", "soil moisture",
                                  "TC", "temperature [K]",
                                  "TD", "TD - tree density [%]",
                                  "TV", "terrain ruggedness",
                                  "CC", "Koppen-Geiger-Climate-Zone",
                                  "PA", "protected areas",
                                  "VD", "VD - Enhanced Vegtation Index"), byrow = T, nc = 2), stringsAsFactors = F)
names(full.names) <- c("pred", "full_name")
#full.names$full_name <- full.names$pred

for(that.model in ls(pattern = "_spatial_model")){
  
  for(iterator in seq(8)){
    i <- c("TD", "VD", "WA", "AR", "SS", "LD", "AD", "NB")[iterator]
    source("Z:/GitRepo/function_file.R")
    dist <- sample.fixed.params(get(that.model), 10)
    dist <- dist[,-c(grep("PA", colnames(dist)), grep("Site", colnames(dist)))]
    
    
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
      interact = c(T, T, T, T, T, F, T, F)[iterator])
    
    
    out <- check
    df <- out[[3]]
    out[[5]] <- merge(out[[5]], full.names, by = "pred")
    out[[4]] <- full.names$full_name[which(full.names$pred %in% out[[4]])]
    names(df) <- "pred"
    comb.out <- cbind(out[[1]][,-4], out[[2]])
    names(comb.out) <- c("lower.in", "mid", "upper.in", "lower.out", "upper.out", "pred")
    
    if(c(T, T, T, T, T, F, T, F)[iterator]){
      comb.out <- data.frame(comb.out, 
                             interactor = rep(out[[5]]$full_name, each = 100), 
                             interactor.value = rep(paste0(out[[5]]$pred, "_", signif(as.numeric(out[[5]]$value), 3)), each = 100))
    }
    
    
    assign(paste0(that.model, "_comb.out_", i), comb.out)
  }
}


for(comb in ls(pattern = "_comb.out_")){
  if(ncol(get(comb)) == 6){
    this_comb <- get(comb)
    this_comb$interactor <- this_comb$interactor.value <- "no interaction" 
    assign(comb, this_comb)
  }
}

ggall <- do.call(rbind, mget(ls(pattern = "_comb.out_")))

ggall$predictor <- NA
ggall$interactor.value <- stringi::stri_replace_all_fixed(ggall$interactor.value, "_", ": ")
ggall$xval <- NA
end <- 0
start <- 1
i = 1
for(i in seq(length(ls(pattern = "comb.out_")))){
  rows <- nrow(get(ls(pattern = "comb.out_")[i]))
  end <- end + rows
  ggall$predictor[start:end] <- rep(do.call(rbind, strsplit(ls(pattern = "_comb.out_"), "_"))[i,8], rows)
  ggall$xval[start:end] <- rep(paste(do.call(rbind, strsplit(ls(pattern = "comb.out_"), "_"))[i,3:4], collapse = "_"), rows)
  start <- start + rows
}

names(full.names) <- c("predictor", "full_name")
ggall$sort <- seq(nrow(ggall))
ggall <- merge(ggall, full.names)
ggall <- ggall[order(ggall$sort), ]

ggall$interactor.value <- factor(ggall$interactor.value, levels = unique(ggall$interactor.value)[c(3, 4, 8, 9, 6, 7, 10, 11, 1, 2, 5)])
ggall <- ggall[ggall$xval != "AGO_Lui",]

effect_plot <- ggplot(ggall, aes(x = pred, y = mid, group = paste0(xval, interactor.value)))+
  geom_ribbon(aes(ymin = lower.in, ymax = upper.in, fill = interactor.value), alpha = .1)+
  #geom_line(aes(color = interactor.value), size = 1)+
  geom_line(aes(color = interactor.value), size = 1)+
  facet_wrap(~full_name, scales = "free_x")+
  scale_fill_manual(values = c("gold1", "gold4",  "brown1", "brown4", "skyblue1", "skyblue4", "springgreen1", "springgreen4", "darkorchid1", "darkorchid4", "black"))+
  scale_color_manual(values = c("gold1", "gold4",  "brown1", "brown4","skyblue1", "skyblue4", "springgreen1", "springgreen4", "darkorchid1", "darkorchid4", "black"))+
  scale_y_log10("habitat suitability")+
  theme_bw()+
  theme(text = element_text(size = 11), legend.position = "top")+
  xlab("predictor value")+
  guides(fill = guide_legend(title = "interaction", nrow = 2),  color = guide_legend(title = "interaction", nrow = 2), alpha = F)

effect_plot


pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/effect_plot_xval.pdf")
effect_plot
dev.off()

which.max(ggall$lower.in)
ggall[1801,]
