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
library(ranger)
library(ggiraphExtra)
rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- readOGR("Z:/GEC/segments_GEE.shp")
if(!"buffer" %in% ls()) buffer <- readOGR("Z:/GEC/buffered_segments.shp")
#source("https://raw.githubusercontent.com/aMilles/GitRepo/master/function_file.R")
source("Z:/GitRepo/function_file.R")

load("Z:/NEMO_out/output_ZWE_MAT_ZWE_ZV_ZWE_SELV_ZWE_SEB_complex_binomial_nonspatial_spatial_xval_LSO_5km_prec1e-04.RData")

dropbox.file = "C:/Users/amilles/Dropbox/"
xy_backup <- xy
#replace scaled and transformed xy values with original values
xy2 <- read.csv(paste0(dropbox.file, "modelling/yxtable.csv"))[,-c(1)]
unique(xy2$Site)
xy2 <- xy2[match(as.character(xy$ID), as.character(xy2$ID)),]
xy[, names(xy2)] <- xy2
if(is.numeric(xy$ID)) xy$ID <- segments_new$ID

if(model.family %in% c("binomial", "zeroinflated.binomial.0", "zeroinflated.binomial.1") & xval){
  link.function <- function(x) exp(x)/(1 + exp(x)) #logit
  if(!"obs" %in% names(xy)) xy$obs <- xy$dnd
}


if(xval){
  identifier <- ifelse(xval.type %in% c("LOSO", "KOSI"), 4, 3)
  spatial_summary <- nonspatial_summary <- vector("list", length(ls(pattern = "_nonspatial_model$")))
  pred.site.names <- NULL
  for(model in ls(pattern = "_nonspatial_model$")){
    pred.site.name <- paste(strsplit(model, "_")[[1]][c(1:identifier)], collapse = "_")
    pred.site <- get(pred.site.name)
    where.na <- which(is.na(pred.site$obs))
    if(xval.type == "KOSI") where.na <- -where.na
    m <- get(model)
    nonspatial_summary[[which(ls(pattern = "_nonspatial_model$") == model)]] <- rbind(summary(m)$fixed[,-7])
    xy[, paste0("nonspatial_pred_", pred.site.name)] <- link.function(m$summary.linear.predictor$`0.5quant`)
    xy[where.na, "nonspatial_pred"] <- link.function(m$summary.linear.predictor$`0.5quant`[where.na])
    pred.site.names <- append(pred.site.names, pred.site.name)
  } 
  for(model in ls(pattern = "_spatial_model$")){
    pred.site.name <- paste(strsplit(model, "_")[[1]][c(1:identifier)], collapse = "_")
    pred.site <- get(pred.site.name)
    where.na <- which(is.na(pred.site$obs))
    if(xval.type == "KOSI") where.na <- -where.na
    m <- get(model)
    spatial_summary[[which(ls(pattern = "_spatial_model$") == model)]] <- rbind(summary(m)$fixed[,-7])
    xy[, paste0("spatial_pred_", pred.site.name)] <- link.function(m$summary.linear.predictor$`0.5quant`)
    xy[where.na, "spatial_pred"] <- link.function(m$summary.linear.predictor$`0.5quant`[where.na])
    pred.site.names <- append(pred.site.names, pred.site.name)
  }
}

if(!xval){
  xy$spatial_pred <- spatial_model$summary.fitted.values$`0.5quant`
  xy$nonspatial_pred <- nonspatial_model$summary.fitted.values$`0.5quant`
  spatial_summary <- list(spatial_model$summary.fixed)
  nonspatial_summary <- list(nonspatial_model$summary.fixed)
  pred.site.names <- as.character(unique(xy$Site))
}


f.inla <- paste0(stringi::stri_replace_all_fixed(f,  "+ f(Site, model = \"iid\")", ""))
f.ranger <- paste0("as.factor(obs) ~ ", paste(unique(strsplit(stringi::stri_replace_all_fixed(stringi::stri_replace_all_fixed(f.inla, ":", " "), "+ ", ""), " ")[[1]]), collapse = " + "))
segments_df$Block <- as.factor(xy$Site)

rF_preds <-  NULL

for(block in unique(xy$Block)){
  rows = which(xy$Block == block)
  rF <- ranger::ranger(f.ranger, xy[-rows,], probability = T, importance = "impurity")
  p <- predict(rF, data = xy[rows,])$predictions[,1]
  rF_preds <- append(rF_preds, p)
}

hist(xy[-rows,]$obs - rF$predictions[,1])

xy$rF_preds <- rF_preds

# for(block in ls(pattern = "xy_without_")[seq(1,length(ls(pattern = "xy_without_")), 3)]) plot(is.na(get(block)$obs))
# m$summary.linear.predictor$`0.5quant`

fit_summary <- read.csv("Z:/residual_analysis/summary.csv", stringsAsFactors = F)[,-1]

#Rough Summary
Efrons <- NULL
for(i in c("spatial_pred", "nonspatial_pred", "rF_preds")){
  print(paste0(i, ": binomial sum ", sum(xy[,i])))
  print(paste0("observed binomial sum ", sum(xy$obs)))
  #plot(xy[,i], xy$obs, pch = "|", main = i)
  Efrons <- append(Efrons, 1 - (1/length(xy$obs) * sum((xy$obs - xy[,i])^2))/(1/length(xy$obs) * sum((xy$obs - mean(xy$obs))^2)))
  print(paste0("Lave/Efron: ", signif(Efrons, 2)))
  print("####################")
}

if(!output_name %in% fit_summary$dataset) fit_summary<- rbind(fit_summary, c(output_name, Efrons, as.character(f.nonspatial)[3]))
write.csv(fit_summary,"Z:/residual_analysis/summary.csv")

spatial_model$summary.fixed
#Fit by Site
out <- NULL
for(i in c("spatial_pred", "nonspatial_pred")){
  for(Site in unique(xy$Site)){
    ss <- xy[xy$Site == Site,]
    out <- append(out, c(1/length(ss$obs) * sum((ss$obs - mean(ss$obs))^2), 
                         1- ((1/length(ss$obs) * sum((ss$obs - ss[,i])^2))/((1/length(ss$obs)) * sum((ss$obs - mean(ss$obs))^2)))
                         )) 
  }
}

fit.df <- data.frame(matrix(out, nc = 2, byrow = T), model = rep(c("spatial", "non-spatial"), each = length(unique(xy$Site))), site = unique(xy$Site))
names(fit.df)[1:2] <- c("detection_ratio", "effron's R-squared")
ggplot(fit.df, aes(x = detection_ratio, y = `effron's R-squared`, group = model, col = site, fill = model, shape = model, color = model))+
  geom_point(size = 5)+
  geom_smooth(method = "lm", se = F)+
  geom_label(aes(label = site))


ggfit <- ggplot(fit.df[fit.df$model == "spatial",], aes(x = detection_ratio, y = `effron's R-squared`))+
  geom_label(aes(label = site),position=position_jitter(width=0.0,height=.1), label.size = 0, alpha = 0)+
  geom_point(size = 5, shape = "+", col = "firebrick")+
  xlab("ratio of detections per subunit")+
  theme_bw()+
  theme(text = element_text(size = 14))
  #geom_smooth(method = "lm", se = F)

pdf(paste0("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/resubstitution_effron", output_name, ".pdf"), width = 9, height = 3)
ggfit
dev.off()
#Variation of coefficients
summaries <- data.frame(rbind(do.call(rbind, nonspatial_summary), do.call(rbind, spatial_summary)))

predictors <- data.frame("predictor" = rep(row.names(nonspatial_summary[[1]]), length(ls(pattern = "nonspatial_model$"))*2), 
           "upper" = summaries$X0.025quant,              
           "value" = summaries$X0.5quant,
           "lower" = summaries$X0.975quant,
           "spatial" = rep(c("non-spatial", "spatial"), each = nrow(spatial_summary[[1]]) * length(ls(pattern = "nonspatial_model$"))),
           "Site" = rep(pred.site.names, each = nrow(spatial_summary[[1]])))

(summary.plot <- 
  ggplot(predictors, aes(x = predictor, y = value, col = Site, ymax = upper, ymin = lower))+
    geom_point()+
    geom_errorbar()+
    facet_wrap(~spatial, nc = 1)+
    ylab("coefficient estimate")+
    ylim(c(-2, 2)))



tf_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/transform_sheet.csv")
reverse <- tf_sheet$name[tf_sheet$value < 0]

for(pred in reverse){
  rows <- (stringi::stri_count_regex(as.character(predictors$predictor), pred)) > 0
    predictors[rows, c(2:4)] <- predictors[rows, c(2:4)] * -1

}
  

pdf(paste0("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/predictor_coefficients_", output_name, ".pdf"), width = 5, height = 8)
(summary.plot <- 
    ggplot(predictors[predictors$spatial == "spatial" & predictors$predictor != "(Intercept)",], aes(x = predictor, y = value, ymax = upper, ymin = lower, col = Site))+
    geom_abline(slope = 0, col = "gray")+
    geom_point()+
    geom_errorbar()+
    ylab("coefficient estimate (95%-interval)")+
    theme_bw()+
    theme(legend.position = "none", text = element_text(size = 14)))+
    coord_flip()
dev.off()


# Spatial analysis
#add observations and predictions to the segments
matchs <- match(as.character(xy$ID), as.character(segments$ID))
summary(matchs)
segments_CC <- segments_new
segments_CC$SC <- NULL
  #segments[match(as.character(xy$ID), as.character(segments$ID)),]

segs <- segments_CC
segs_cent <- geosphere::centroid(segs)
segs <- cbind(long = segs_cent[,1], lat = segs_cent[,2], data.frame(segs))
segments_df <- cbind(data.frame(segs_cent), xy, data.frame(segments_CC))
names(segments_df)[c(1,2)] <- c("long", "lat") 

#create spatial correlogram
spatial_dep <- segments_df
spatial_dep <- subset(spatial_dep, Site == "ZWE_MAT")
spatial_dep_non_spatial_m <- ncf::correlog(x = spatial_dep$long, y = spatial_dep$lat, z = spatial_dep$nonspatial_pred - spatial_dep$obs, increment = 5, resamp = 2, latlon = T)
spatial_dep_spatial_m <- ncf::correlog(x = spatial_dep$long, y = spatial_dep$lat, z = spatial_dep$spatial_pred - spatial_dep$obs, increment = 5, resamp = 2, latlon = T)


source("Z:/GitRepo/function_file.R")


#save those plots in one pdf

hist(segments_df$rF_pred)
summary(segments_df$rF_pred)
segments_df$prediction_log10 <- log10(segments_df$spatial_pred)
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/", output_name, "predictions_rF.pdf"), onefile = T)
for(Site in as.character(unique(segments_df$Site))){
 print(
   map <- map_site(df = segments_df[segments_df$Site == Site,], o = "obs", p = "rF_preds", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-0.01, 1)) +
     geom_point(data = segments_df[segments_df$obs == 1 & segments_df$Site == Site, c("long", "lat", "Block")], aes(x = long, y = lat), shape = "+")+
     guides(fill = guide_colorbar(title = "Habitat Suitability"))+
     ggtitle("")+
     theme(text = element_text(size = 14))
   )
}
dev.off()




samples <- t(as.matrix(sample.fixed.params(xy_without_block10_spatial_model, nsample = 100)))
xy.m <- xy_backup
xy.m$Site[] <- "ZWE_MAT"
m <- model.matrix(as.formula(paste0("obs ~", f.inla)), xy_backup)
prob <- invlogit(m %*% samples)
prob <- apply(data.frame(prob), 1, median)

segments_df$prediction_log10 <- log10(segments_df$spatial_pred_xy_without_ZWE_SELV)
segments_df$inla.pred <- prob 

all.map <- map_site(df = segments_df, o = "obs", p = "inla.pred", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-0.01, 1))+
  geom_point(data = segments_df[segments_df$obs == 1, c("long", "lat", "Block")], aes(x = long, y = lat), shape = 20, col = "white", fill = "black", alpha = 0.8, size = .5)+
  guides(fill = guide_colorbar(title = "Habitat Suitability"))+
  scale_fill_gradientn(colors = viridis::cividis(5))+
  ggtitle("")+
  theme(text = element_text(size = 14))

pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/", output_name, "all_predictions_nositeeffect.pdf"), onefile = T)
all.map
dev.off()



map.preds <- c("AD", "LD", "TV", "AR", "SM", "VD", "WA", "SS", "SC", "TC", "SL")

for(i in seq(10)) dev.off()
{
  pdf(file = paste0("Z:/residual_analysis/",tools::file_path_sans_ext(output_name), ".pdf"), onefile = T,paper = "a4", width = 8.27, height = 11.69)
    for(predictor in map.preds){
      print(map_site(df = segments_df, o = "obs", p = predictor, plot.what = "p", title = predictor, lower = F, plot.range = range(segments_df[, predictor]))+
              geom_point(data = segments_df[segments_df$obs == 1, c("long", "lat", "Block")], aes(x = long, y = lat), shape = "+"))
    }
  dev.off()
}

#plot all covariates
for(i in seq(10)) dev.off()
{
  pdf(file = paste0("Z:/residual_analysis/", selection, "_predictors.pdf"), onefile = T,paper = "a4", width = 8.27, height = 11.69)
  
  for(pred in names(segments_df[,3:21])[which(!names(segments_df[,3:21]) %in% c("ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
    map_site(df = segments_df, o = "obs", p = pred, plot.what = "p", title = pred)
  }
  
  dev.off() 
}


#simple "effect" plots
gg.xy <- cbind(xy_without_ZWE_MAT[, rev(strsplit(f, " +")[[1]][-seq(2, length(strsplit(f, " +")[[1]]), 2)])[-1]], xy[,c("obs", "nonspatial_pred", "spatial_pred")])
gg.xy <- reshape2::melt(gg.xy, id.vars = c("nonspatial_pred", "spatial_pred", "obs"))

  ggplot(gg.xy, aes(x = value, y = obs))+
    geom_smooth()+
    facet_wrap(~variable, scales = "free")

  
#elaborated effect plots  
  if(effect){
    tf_sheet <- read.csv("Z:/modelling/transform_sheet.csv", stringsAsFactors = F)
    scale_sheet <- read.csv("Z:/modelling/scale_sheet.csv", stringsAsFactors = F)
    tf_sheet$name[is.na(tf_sheet$name)] <- "NA." -> scale_sheet$name[is.na(scale_sheet$name)]
    #scale_sheet$name <- tf_sheet$name
    df <- read.csv("Z:/modelling/yxtable_scaled_transformed.csv")
    names(df)
    bt_mins <- us_mins <- mins <- apply(df[,tf_sheet$name], 1, min)
    bt_maxs <- us_maxs <- maxs <- apply(df[,tf_sheet$name], 1, max)
    for(i in length(mins)){
      us_mins[i] <- mins[i] * scale_sheet$scale[i] + scale_sheet$center[i]
      us_maxs[i] <- maxs[i] * scale_sheet$scale[i] + scale_sheet$center[i]
    }
    
  }  


  
  
  sample.fixed.params <- function(model = model_test, nsample = 10000, param = "") {
    params <- names(model$marginals.fixed)
    which.param <- grep(param, params)
    
    if(sum(which.param)==0) stop("Check parameter names")
    dist <- matrix(ncol = length(which.param), nrow = nsample)
    for(i in 1:NROW(which.param)) {
      dist[,i] <- sample(model$marginals.fixed[[which.param[i]]][,1], nsample, TRUE,
                         model$marginals.fixed[[which.param[i]]][,2])
      
    } 
    return(dist)
  }
  sample.fix  

  
  
## ggradar

xy_radar <- read.csv("Z:/modelling/yxtable.csv")  
rm(list = ls(pattern = "radar_"))  
for(Site in c("BWA_NOR", "KEN_LAI", "KEN_TSV", "XWA_TBC", "ZWE_MAT", "ZWE_ZV", "ZWE_SELV")) assign(paste0("radar_", Site), ggiraphExtra::ggRadar(data = xy_radar[xy_radar$Site == Site, c(splines, "HT", "Site")], aes(color = HT), scales = "free", rescale = T))
ggiraphExtra::ggRadar(data = xy_radar[xy_radar$Site == Site, c(splines, "HT")], aes(color = HT), scales = "free", rescale = T, interactive = T)+ggtitle(Site)+theme(legend.position = "none")

do.call(gridExtra::grid.arrange, mget(ls(pattern = "radar_")))
