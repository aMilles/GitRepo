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
rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- readOGR("Z:/GEC/segments_GEE.shp")
if(!"buffer" %in% ls()) buffer <- readOGR("Z:/GEC/buffered_segments.shp")



#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_5km.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LSO_5km.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LOSO_5km.RData")
#load("Z:/NEMO_out/output_KEN_ZWE_simple_binomial_nonspatial_spatial_xval_LOSO_5km.RData")
#load("Z:/NEMO_out/output_all_simple_binomial_nonspatial_spatial_xval_LSO_5km.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LOSO_5km_splines.RData")
#load("Z:/NEMO_out/output_ZWE_simple_zeroinflated.binomial.0_nonspatial_spatial_xval_LOSO_5km.RData")
#load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_5km_splines.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_5km_splines.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LSO_5km.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_5km_splines.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LOSO_5km_splines.RData")
load("Z:/NEMO_out/output_BWA_NOR_KEN_LAI_KEN_TSV_XWA_TBC_ZWE_MAT_ZWE_ZV_ZWE_SELV_complex_binomial_nonspatial_spatial_5km_splines.RData")

xy_backup <- xy

#replace scaled and transformed xy values with original values
xy2 <- read.csv("Z:/modelling/yxtable.csv")[,-c(1)]
unique(xy2$Site)
xy2 <- xy2[match(as.character(xy$ID), as.character(xy2$ID)),]
xy[, names(xy2)] <- xy2
if(is.numeric(xy$ID)) xy$ID <- segments_new$ID

if(model.family %in% c("binomial", "zeroinflated.binomial.0", "zeroinflated.binomial.1") & xval){
  link.function <- function(x) exp(x)/(1 + exp(x)) #logit
  if(!"obs" %in% names(xy)) xy$obs <- xy$dnd
}


if(xval){
  identifier <- ifelse(xval.type %in% c("Site", "LOSO"), 4, 3)
  spatial_summary <- nonspatial_summary <- vector("list", length(ls(pattern = "_nonspatial_model$")))
  pred.site.names <- NULL
  for(model in ls(pattern = "_nonspatial_model$")){
    pred.site.name <- paste(strsplit(model, "_")[[1]][c(1:identifier)], collapse = "_")
    pred.site <- get(pred.site.name)
    where.na <- which(is.na(pred.site$obs))
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
  spatial_summary <- spatial_model$summary.fixed
  nonspatial_summary <- nonspatial_model$summary.fixed
}


# for(block in ls(pattern = "xy_without_")[seq(1,length(ls(pattern = "xy_without_")), 3)]) plot(is.na(get(block)$obs))
# m$summary.linear.predictor$`0.5quant`

fit_summary <- read.csv("Z:/residual_analysis/summary.csv", stringsAsFactors = F)[,-1]

#Rough Summary
Efrons <- NULL
for(i in c("spatial_pred", "nonspatial_pred")){
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
ggplot(fit.df, aes(x = detection_ratio, y = `effron's R-squared`, group = model, col = site))+
  geom_point(size = 5)+
  geom_smooth(method = "lm")+
  facet_wrap(~model, nc = 1, scales = "free")
fit.df[fit.df$detection_ratio > 0.1,]
dev.off()

#Variation of coefficients
summaries <- data.frame(rbind(do.call(rbind, nonspatial_summary), do.call(rbind, spatial_summary)))

predictors <- data.frame("predictor" = rep(row.names(nonspatial_summary[[1]]), length(ls(pattern = "nonspatial_model$"))*2), 
           "upper" = summaries$X0.025quant,              
           "value" = summaries$X0.5quant,
           "lower" = summaries$X0.975quant,
           "spatial" = rep(c("non-spatial", "spatial"), each = nrow(spatial_summary[[1]]) * length(ls(pattern = "nonspatial_model$"))),
           "Site" = rep(pred.site.names, each = nrow(spatial_summary[[1]])*2))

(summary.plot <- 
  ggplot(predictors, aes(x = predictor, y = value, col = Site, ymax = upper, ymin = lower))+
    geom_point()+
    geom_errorbar()+
    facet_wrap(~spatial, nc = 1)+
    ylab("coefficient estimate"))

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


### FUNCTION ###
# Plot a Map of residuals, observations or predictions..
map_site <- function(df = segments_df, p = "spatial_pred", o = "obs", plot.what = "r", filter = NA, title = NA, interpolate = T, buffer.it = T, buffer.spdf = buffer){
  
  #use unique names for plotting
  df <- df[,match(c("long","lat", "Site", p,o, "HT", "Block"), names(df))]
  names(df)[c(4,5)] <- c("p", "o")
  df$r <- df$o - df$p
  df$input <- df[,match(plot.what, names(df))]
  df$input[sample(nrow(df), 2)] <- c(ifelse(plot.what == "r", -1, 0), 1)

  #create a filter, to exclude extreme values, that may affect color ranges
  if(!all(is.na(filter))){
    df$input[df$input > quantile(df$input, filter[2], na.rm = T)] <- NA
    df$input[df$input < quantile(df$input, filter[1], na.rm = T)] <- NA
  }
  
  name <- ifelse(plot.what == "r", "Residuals", get(plot.what))
  
  #take a subsample of the data for residuals vs. predicted plot
  samples.df <- c(which(df$o == 1)[round(seq(1, length(which(df$o == 1)), length = 1000))],
                  which(df$o == 0)[round(seq(1, length(which(df$o == 1)), length = 1000))])

  
  #create a interpolation surface
  for(site in unique(df$Block)){
    map.subset <- df[segments_df$Site == site, c("long", "lat", "input")]
      
    if(!buffer.it){
      mba <- MBA::mba.surf(map.subset, 300, 300, extend = FALSE)
      dimnames(mba$xyz.est$z) <- list(mba$xyz.est$x, mba$xyz.est$y)
      out <- data.frame(reshape2::melt(mba$xyz.est$z, varnames = c('long', 'lat'), value.name = "input"), Block = site)
    }
      
    if(buffer.it){
      out <- raster(MBA::mba.surf(map.subset, 300, 300, sp = T)$xyz.est)
      masked_mba <- mask(out, buffer.spdf)
      out <- data.frame(as.data.frame(as(masked_mba, "SpatialPixelsDataFrame")), Block = site)
    }
      
    names(out) <- c("input", 'long', 'lat', "Block")
    assign(paste0("ipblock_", site), out)
      
  }
    
  df_map <- cbind(do.call(rbind, mget(ls(pattern = "ipblock_"))))
    
  gg_map <- 
    ggplot()+
    geom_raster(data = df_map, aes(x = long, y = lat, fill = input))+
    coord_equal()+
    scale_fill_viridis()+
    facet_wrap(~Block, scales = "free")+
    theme_light()+
    theme(legend.position="bottom", plot.background = element_rect(color = "gray")) +
    theme(legend.key.width=unit(2, "cm"))+
    guides(fill = guide_legend(title = name))+
    xlab("Longitude [°]")+
    ylab("Latitude [°]")+
    ggtitle(name)
  
    
  #print("create hist")
  if(plot.what == "r"){
    gg_resid_hist <- 
      ggplot()+
      geom_histogram(data = df, aes(x = input, col = o, fill = o), bins = 100)+
      scale_y_log10()+
      xlab("residuals (o - p)")+
      ylab("log10(count)")+
      theme_gray()+
      theme(legend.position = "none")+
      xlim(c(-1,1))
    
    #print("create scatter")
    gg_resid_scatter <- 
      ggplot()+
      geom_point(data = df[samples.df, ], aes(x = r, y = p), shape = "|", alpha = .2)+
      xlab("residuals (o - p)")+
      ylab("predicted")+
      xlim(c(-1,1))+
      theme_gray()
 
     gg_map <- 
       gridExtra::grid.arrange(gg_map, 
                               gg_resid_hist, 
                               gg_resid_scatter,
                               layout_matrix = matrix(c(1,1,1,2,1,1,1,3), nc = 2),
                               top = title)
  }else{
    
    gg_hist <- 
      ggplot(data = df, aes(input, fill = HT, color = HT, alpha = .5))+
        geom_density()+
        guides(fill = guide_legend(title = "Observed"), alpha = F, col = F)+
        scale_fill_manual(values = c("firebrick", "turquoise", "white"))+
        facet_wrap(~Block, scales = "free")

    
    gg_map <- 
      gridExtra::grid.arrange(gg_map, 
                              gg_hist,
                              layout_matrix = matrix(c(1,1,1,2,2,1,1,1,2,2), nc = 2),
                              top = title)

  }
  return(gg_map)
}

#save those plots in one pdf
segments_df$Block <- segments_df$Site
segments_df$WA <- segments_df$WA/1000
map_site(df = segments_df, o = "obs", p = "spatial_pred", plot.what = "p", title = "non-spatial model", buffer.it = T)
map_site(df = segments_df, o = "obs", p = "VD", plot.what = "p", title = "non-spatial model", buffer.it = T)

names(test) <- c("input", "long", "lat", "Block")

ggplot()+
  geom_raster(data = test, aes(x = long, y = lat, fill = input))+
  coord_equal()+
  scale_fill_viridis()+
  facet_wrap(~Block, scales = "free")+
  theme_light()+
  theme(legend.position="bottom", plot.background = element_rect(color = "gray")) +
  theme(legend.key.width=unit(2, "cm"), legend.title = name)+
  xlab("Longitude [°]")+
  ylab("Latitude [°]")


for(i in seq(10)) dev.off()
{
  pdf(file = paste0("Z:/residual_analysis/",tools::file_path_sans_ext(output_name), ".pdf"), onefile = T,paper = "a4", width = 8.27, height = 11.69)
    plot(spatial_dep_non_spatial_m$mean.of.class, spatial_dep_non_spatial_m$correlation, type = "l", lty = "dashed", main = "correlogram, ZWE_MAT")
    lines(x = spatial_dep_spatial_m$mean.of.class, y = spatial_dep_spatial_m$correlation, lty = "solid", col = "green")
    map_site(df = segments_df, o = "obs", p = "spatial_pred", plot.what = "r", title = "Residuals: spat model")
    map_site(df = segments_df, o = "obs", p = "nonspatial_pred", plot.what = "r", title = "Residuals: nspat model")
    for(predictor in rev(strsplit(f, " +")[[1]][-seq(2, length(strsplit(f, " +")[[1]]), 2)])[-1]){
      map_site(df = segments_df, o = "obs", p = predictor, plot.what = "p", title = predictor)
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

