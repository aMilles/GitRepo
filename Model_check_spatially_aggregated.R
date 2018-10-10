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

#load("Z:/NEMO_out/output_ZWE_MAT_ZWE_ZV_ZWE_SELV_ZWE_SEB_complex_binomial_nonspatial_spatial_xval_LSO_5km_prec1e-04.RData")
load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_6km_prec1e-04.RData")
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


f.inla <- paste0(stringi::stri_replace_all_fixed(f,  "+ f(Site ,model=\"iid\") + f(PA ,model=\"iid\")", ""))
f.ranger <- paste0("as.factor(obs) ~ ", paste(unique(strsplit(stringi::stri_replace_all_fixed(stringi::stri_replace_all_fixed(f.inla, ":", " "), "+ ", ""), " ")[[1]]), collapse = " + "))


rF_preds <-  NULL
if(all(is.na(xy$Block))) xy$Block <- xy$Site
for(block in unique(xy$Block)){
  rows = which(xy$Block == block)
  rF <- ranger::ranger(f.ranger, xy[-rows,], probability = T, importance = "impurity")
  p <- predict(rF, data = xy[rows,])$predictions[,1]
  rF_preds <- append(rF_preds, p)
}

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
  geom_label(aes(label = site), position=position_jitter(width=0.0,height=.1), label.size = 0, alpha = 0)+
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
           "Site" = rep(c("nonspatial", "spatial"), each = nrow(spatial_summary[[1]])))

(summary.plot <- 
  ggplot(predictors, aes(x = predictor, y = value, col = Site, ymax = upper, ymin = lower))+
    geom_point()+
    geom_errorbar()+
    facet_wrap(~spatial, nc = 1)+
    ylab("coefficient estimate")+
    ylim(c(-2, 2)))


tf_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/transform_sheet.csv")
reverse <- tf_sheet$name[tf_sheet$value < 0]
reverse
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


##################################################################
################### SPATIAL PREDICTIONS ##########################
##################################################################

#add observations and predictions to the segments
matchs <- match(as.character(xy$ID), as.character(segments$ID))
segments_CC <- segments_new
segments_CC$SC <- NULL
segs <- segments_CC
segs_cent <- geosphere::centroid(segs)
segs <- cbind(long = segs_cent[,1], lat = segs_cent[,2], data.frame(segs))
segments_df <- cbind(data.frame(segs_cent), xy, data.frame(segments_CC))
segments_df$Block <- as.factor(xy$Site)
names(segments_df)[c(1,2)] <- c("long", "lat") 

source("Z:/GitRepo/function_file.R")


#save those plots in one pdf

### PREDICTIONS WITH RANDOM EFFECTS (SPATIAL AND NON-SPATIAL) ####

#spatial

segments_df$prediction_log10 <- log10(segments_df$spatial_pred)
segments_df$prediction_log10nsp <- log10(segments_df$nonspatial_pred)

pred_map <- map_site(df = segments_df, o = "obs", p = "prediction_log10", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-5, 0))

pred_map2 <- pred_map+
  scale_fill_gradientn(colors = viridis::inferno(5, end = 0.75), limits = c(-4, 0.1), na.value = viridis::inferno(5, end = 0.75)[1])+
  guides(fill = guide_colorbar(title = "habitat suitability (log10)"), color = F)+
  geom_point(data = segments_df[segments_df$obs == 1, c("long", "lat", "Block", "prediction_log10")], aes(x = long, y = lat), shape = "O", size = 0.01, fill = "black", color = "white")+
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 9, angle = 30), axis.text.y = element_text(size = 9, angle = 30), panel.background = element_rect(color = "gray50"), plot.background = element_rect(color = "gray50"))+
  ggtitle("")

pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/predictions_spatial_model.pdf"), onefile = T)
pred_map2
dev.off()

#non-spatial

pred_mapnsp <- map_site(df = segments_df, o = "obs", p = "prediction_log10nsp", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-5, 0))

pred_map2nsp <- pred_mapnsp+
  scale_fill_gradientn(colors = viridis::inferno(5, end = 1), limits = c(-4, 0.1), na.value = viridis::inferno(5, end = 1)[1])+
  guides(fill = guide_colorbar(title = "habitat suitability (log10)"), color = F)+
  geom_point(data = segments_df[segments_df$obs == 1, c("long", "lat", "Block", "prediction_log10nsp")], aes(x = long, y = lat), shape = "O", size = 0.01, fill = "black", color = "white")+
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 9, angle = 30), axis.text.y = element_text(size = 9, angle = 30), panel.background = element_rect(color = "gray50"), plot.background = element_rect(color = "gray50"))+
  ggtitle("")

pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/predictions_spatial_model.pdf"), onefile = T)
pred_map2nsp
dev.off()


### PREDICT WITHOUT RANDOM EFFECTS ####


if(!"df_backup" %in% ls()) df_backup <- segments_df
df_backup -> segments_df

segments_df$Block <- xy$Site 
samples <- t(as.matrix(sample.fixed.params(spatial_model, nsample = 100)))


# PREDICT HS AT THE TIMES OF OBSERVATIONS

xy.pred <- xy_backup
m <- model.matrix(as.formula(paste0("obs ~", f.inla)), xy.pred)
inla.pred <- invlogit(m %*% samples)
inla.pred <- apply(data.frame(inla.pred), 1, median)
segments_df$inla.pred <- inla.pred
segments_df$log10inla.pred <- log10(inla.pred)


segments_df <- merge(segments_df, setNames(aggregate(obs ~ Site, xy, sum), c("Block", "detection_sum")))
names(segments_df)
inla.pred.map <- map_site(df = segments_df, o = "obs", p = "log10inla.pred", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-5, 0))
 
 pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/Map_nonspatial_pred.pdf"))
 inla.pred.map +
   scale_fill_gradientn(colors = viridis::inferno(5, end = 0.75), limits = c(-4, -1), na.value = viridis::inferno(5, end = 0.75)[1])+
   geom_point(data = segments_df[segments_df$obs == 1, c("long", "lat", "Block", "inla.pred", "detection_sum")], aes(x = long, y = lat, size = 1/(detection_sum), alpha = 1/detection_sum^3), fill = "black", color = "gray90", shape = 3)+
   scale_size_continuous(range = c(0.0000001, 2))+
   scale_alpha_continuous(range = c(0.3,1))+
   guides(fill = guide_colorbar(title = "log10(habitat suitability)"), color = F, size = F, alpha = F)+
   theme(text = element_text(size = 10), axis.text.x = element_text(size = 9, angle = 30), axis.text.y = element_text(size = 9, angle = 30), panel.background = element_rect(color = "gray50"), plot.background = element_rect(color = "gray50"))+
   ggtitle("")

dev.off()
###predict at the time of minimum HS  

xy.m <- xy_backup
xy.m$Block <- NULL

#read predictors at the time of minimum HS and replace them
xy_season <- read.csv(paste0(dropbox.file, "modelling/xytable_minHS.csv"))
all(as.character(xy_season$ID) == as.character(xy.m$ID))
xy.m[, c("VD","WA","TC", "SC")] <- xy_season[match(as.character(xy.m$ID), as.character(xy_season$ID)), c("VD","WA","TC", "SC")]
xy_season <- na.omit(xy.m)

#there were some more NAs in the seasonal dataset
segments_df <- segments_df[na.omit(match(as.character(xy_season$ID), as.character(segments_df$ID))), ]

# create predictions

xy.minHS <- xy_season
m <- model.matrix(as.formula(paste0("obs ~", f.inla)), xy.minHS)
minHSprob <- invlogit(m %*% samples)
minHSprob <- apply(data.frame(minHSprob), 1, median)
segments_df$minHSprob <- minHSprob
segments_df$log10minHSprob <- log10(minHSprob)

pred_minHS_minHS_map <- map_site(df = segments_df, o = "obs", p = "log10minHSprob", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-5, 0))
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/Map_noLD.pdf"), onefile = T)
pred_minHS_minHS_map+
  scale_fill_gradientn(colors = viridis::inferno(5, end = 0.75), limits = c(-4, -1), na.value = viridis::inferno(5, end = 0.75)[1])+
  guides(fill = guide_colorbar(title = "log10(habitat suitability)"), color = F)+
  theme(text = element_text(size = 10), axis.text.x = element_text(size = 9, angle = 30), axis.text.y = element_text(size = 9, angle = 30))+
  ggtitle("")
dev.off()


#generate a scenario with environmental effects set to the same level
xy.anthro <- xy_season
xy.anthro$TD[] <-  0
xy.anthro$VD[] <-  0
xy.anthro$WA[] <-  0
xy.anthro$SS[] <- 0
xy.anthro$NB[] <- 0
m <- model.matrix(as.formula(paste0("obs ~", f.inla)), xy.anthro)
anthro <- invlogit(m %*% samples)
anthro <- apply(data.frame(anthro), 1, median)
segments_df$anthro <- anthro
segments_df$log10anthro <- log10(anthro)


#generate a scenario with anthropogenic effects set to the same level
xy.enviro <- xy_season
xy.enviro$HD <- 0
xy.enviro$LD <- 0
xy.enviro$AD <- 0
xy.enviro$AR <- 0
m <- model.matrix(as.formula(paste0("obs ~", f.inla)), xy.enviro)
enviro <- invlogit(m %*% samples)
enviro <- apply(data.frame(enviro), 1, median)
segments_df$enviro <- enviro
segments_df$log10enviro <- log10(enviro)


xy.nono <- xy_season
xy.nono$HD <- min(xy.nono$HD)
xy.nono$LD <- min(xy.nono$LD)
xy.nono$AD <- min(xy.nono$AD)
xy.nono$AR <- min(xy.nono$AR)
xy.nono$SS[] <- 0
xy.nono$NB[] <- 0
xy.nono$TD[] <-  0
xy.nono$VD[] <-  0
xy.nono$WA[] <-  0
m <- model.matrix(as.formula(paste0("obs ~", f.inla)), xy.nono)
nono <- invlogit(m %*% samples)
nono <- apply(data.frame(nono), 1, median)
segments_df$nono <- nono
segments_df$log10nono <- log10(nono)



pred_enviro <- map_site(df = segments_df, o = "obs", p = "enviro", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-5, 0))
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/pred_enviro.pdf"), onefile = T)
pred_enviro+
  scale_fill_gradientn(colors = viridis::inferno(5), limits = c(-3, -1), na.value = viridis::inferno(5)[1])+
  guides(fill = guide_colorbar(title = "habitat suitability"), color = F)+
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 9, angle = 30), axis.text.y = element_text(size = 9, angle = 30))+
  ggtitle("")
dev.off()

pred_anthro <- map_site(df = segments_df, o = "obs", p = "anthro", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-5, 0))
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/pred_anthro.pdf"), onefile = T)
pred_anthro+
  scale_fill_gradientn(colors = viridis::inferno(5), limits = c(-3, -1), na.value = viridis::inferno(5)[1])+
  guides(fill = guide_colorbar(title = "habitat suitability"), color = F)+
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 9, angle = 30), axis.text.y = element_text(size = 9, angle = 30))+
  ggtitle("")
dev.off()



#### aggreagte LAND USE CHANGE ####
comp.LUCCdf <- segments_df[, c("enviro", "PA", "Site", "sampef")]
check <- aggregate(LUCC_perc ~ Site, FUN = spatstat::weighted.median, w = 1/comp.LUCCdf$sampef,  data = comp.LUCCdf)
#comp.LUCC <- aggregate(LUCC_perc ~ PA + Site, FUN = spatstat::weighted.median, w = 1/comp.LUCCdf$sampef, data = comp.LUCCdf)
comp.LUCC <- aggregate(LUCC_perc ~ PA + Site, FUN = median, data = comp.LUCCdf)
names(comp.LUCC) <- c("PA", "Site", "LUCC_perc")

comp.LUCC$Site <- factor(comp.LUCC$Site, levels = check$Site[order(check$LUCC_perc)])

pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/change_habitatsuitability.pdf"), onefile = T)
ggplot(comp.LUCC, aes(x = Site, y = LUCC_perc, fill = as.factor(PA)))+
  geom_bar(stat="identity",position="dodge")+
  coord_flip()+
  scale_fill_manual(name="",
                      breaks=c(0, 1),
                      labels=c( "not protected", "protected"), 
                      values = c("azure3", "darkslategray4"))+
  theme_bw()+
  ylab("median change of habitat suitability [%]")+
  xlab("")+
  theme(legend.position = "top", text = element_text(size = 20), panel.grid.major.y = element_line(size = 1, colour = c("gray90")), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
dev.off()

a <- map_site(df = segments_df, o = "obs", p = "LUCC_perc", plot.what = "p", title = paste0("Prediction: ", Site), buffer.it = T, lower = F, plot.range = c(-5, 0))+
  scale_fill_gradientn(colors = viridis::inferno(5), limits = c(-60, 5), na.value = viridis::inferno(5)[1])+
  guides(fill = guide_colorbar(title = "change in habitat suitability [%]"))+
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 9, angle = 30), axis.text.y = element_text(size = 9, angle = 30))+
  ggtitle("")
  
a+
  scale_fill_gradientn(colors = viridis::inferno(5), limits = c(-60, 5), na.value = viridis::inferno(5)[1])+
  guides(fill = guide_colorbar(title = "change in habitat suitability [%]"))+
  theme(text = element_text(size = 12), axis.text.x = element_text(size = 9, angle = 30), axis.text.y = element_text(size = 9, angle = 30))+
  ggtitle("")
#### Sensitivity / Elasticity ####
full.names <- data.frame(matrix(c("WA", "WA - distance to water [km]",
                                  "AD", "AD - agricultural density [%]",
                                  "HD", "HD - human density [people/km²]",
                                  "LD", "LD - livestock density [TLU/km²]",
                                  "AR", "AR - accessibility from roads [h]",
                                  "VD", "VD - Enhanced Vegtation Index"), byrow = T, nc = 2), stringsAsFactors = F)
names(full.names) <- c("predictor", "full_name")

vary.pred <- function(df, predictor, change, formula, posterior){
  tf_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/transform_sheet.csv", stringsAsFactors = F)
  scale_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/scale_sheet.csv", stringsAsFactors = F)
  
  bt_value <- tf_sheet$value[tf_sheet$name == predictor]
  shift <- tf_sheet$shift[tf_sheet$name == predictor]
   
  scale <- scale_sheet$scale[scale_sheet$name == predictor]
  center <- scale_sheet$center[scale_sheet$name == predictor]
  
  par(mfrow = c(3,1))
  cov <- df[,predictor]
  changer <- (((cov * scale + center)^(1/bt_value)) + shift) * change
  hist(changer, main = predictor)
  boxplot(changer ~ df$Site)
  changer <- (((changer - shift)^bt_value) - center) / scale
  hist(changer, main = predictor)
  df[, predictor] <- changer
  m <- model.matrix(as.formula(paste0("obs ~", formula)), df)
  prob <- apply(data.frame(invlogit(m %*% posterior)), 1, function(x) median(x, na.rm = T))
  return(prob)
}
sampling.effort <- read.csv(paste0(dropbox.file, "modelling/sampling_effort.csv"), row.names = 1)
segments_df$sampef <- sampling.effort[match(segments_df$ID, sampling.effort$ID), 2]


samples <- t(as.matrix(sample.fixed.params(spatial_model, nsample = 50)))
preds <- c("WA", "VD", "AR", "LD", "AD")
change <- c(1.1, 1.1, .9, 1.1, 1.1, 1.1)
dir <- rep(c("up"), length(preds))
out <- vector("list", length(preds))
for(pred in seq(length(preds))) out[[pred]] <- vary.pred(df = xy_season, predictor = preds[pred], change = change[pred], formula = f.inla, posterior = samples)

ref <- vary.pred(df = xy_season, predictor = preds[pred], change = 1, formula = f.inla, posterior = samples)
change <- data.frame(prob = 100 *(do.call(c, out) - ref)/ref, predictor =  rep(preds, each = nrow(xy_season)), direction = "up", site = xy_season$Site)
aggchange <- aggregate(prob ~ predictor + site, FUN = spatstat::weighted.median, w = 1/segments_df$sampef, data = change)
aggchange <- merge(aggchange, full.names)

weird <- scales::trans_new("signed_log",
                           transform=function(x) sign(x)*log10(abs(x)),
                           inverse=function(x) sign(x)*(10^abs(x)))
aggchange
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/change_habitatsuitability_elasticity.pdf"), onefile = T)
ggplot(aggchange, aes(y = prob, x = site, fill = prob > 0))+
  geom_abline(intercept = 0, slope = 0)+
  geom_bar(stat="identity",position="dodge", color = "black")+
  facet_wrap(~full_name, scales = "free_y", nrow = 6)+
  scale_fill_manual(name="",
                    breaks=c(0, 1),
                    labels=c( "", ""), 
                    values = c("gray90", "gray50"))+
  theme_bw()+
  ylab("median change of habitat suitability [%]")+
  xlab("site")+
  theme(legend.position = "none", text = element_text(size = 14), panel.grid.major.y = element_line(size = .1, colour = "darkgray"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 30, hjust = .9))
dev.off()


gg <- data.frame(rbind(segments_df[, c("Site", "nono", "anthro")], setNames(segments_df[, c("Site", "nono", "enviro")], c("Site", "nono", "anthro"))), conditions = rep(c("anthropogenic", "environmental"), each = nrow(segments_df)))
gg <- data.frame(rbind(segments_df[, c("Site", "nono", "anthro", "PA")]), conditions = rep(c("anthropogenic"), each = nrow(segments_df)))

gg$Site <- factor(gg$Site, as.character(unique(gg$Site)[order(aggregate(anthro ~ Site, gg, median)$anthro)]))

pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/anthropogenic_effect.pdf"), width = 8, height = 3.5)
ggplot(gg)+
  geom_boxplot(aes(y = 100 * anthro, x = Site, fill = as.factor(PA)))+
  scale_y_continuous("habitat suitability [%]", breaks = c(0.25, 0.5, 1, 2))+
  theme_bw()+
  scale_fill_manual(name="",
                    breaks=c(0, 1),
                    labels=c( "not protected", "protected"), 
                    values = c("azure3", "darkslategray4"))+
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 30, hjust = 0.9), legend.position = "top")
dev.off()


