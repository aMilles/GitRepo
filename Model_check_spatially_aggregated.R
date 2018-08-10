library(INLA)
library(rgdal)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(spdep)
library(MBA)

rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- readOGR("Z:/GEC/segments_GEE.shp")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_5km.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LSO_5km.RData")
#load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LOSO_5km.RData")
#load("Z:/NEMO_out/output_KEN_ZWE_simple_binomial_nonspatial_spatial_xval_LOSO_5km.RData")
#load("Z:/NEMO_out/output_all_simple_binomial_nonspatial_spatial_xval_LSO_5km.RData")
load("Z:/NEMO_out/output_ZWE_simple_binomial_nonspatial_spatial_xval_LOSO_5km_splines.RData")
#load("Z:/NEMO_out/output_ZWE_simple_zeroinflated.binomial.0_nonspatial_spatial_xval_LOSO_5km.RData")

xy_backup -> xy

#replace scaled and transformed xy values with original values
xy2 <- read.csv("Z:/modelling/yxtable.csv")[,-c(1)]
unique(xy2$Site)
xy2 <- xy2[match(as.character(xy$ID), as.character(xy2$ID)),]
unique(xy2$Site)
xy[, 1:21] <- xy2
xy <- xy[,1:22]

if(model.family %in% c("binomial", "zeroinflated.binomial.0", "zeroinflated.binomial.1") & xval){
  link.function <- function(x) exp(x)/(1 + exp(x)) #logit
  if(!"obs" %in% names(xy)) xy$obs <- xy$dnd
}
m$logfile
identifier <- ifelse(xval.type %in% c("Site", "LOSO"), 4, 3)
plot(m$summary.linear.predictor$`0.5quant`)
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

# for(block in ls(pattern = "xy_without_")[seq(1,length(ls(pattern = "xy_without_")), 3)]) plot(is.na(get(block)$obs))
# m$summary.linear.predictor$`0.5quant`


#Rough Summary
for(i in c("spatial_pred", "nonspatial_pred")){
  print(paste0(i, ": binomial sum ", sum(xy[,i])))
  print(paste0("observed binomial sum ", sum(xy$obs)))
  plot(xy[,i], xy$obs, pch = "|", main = i)
  print(paste0("Spearman: ",cor(xy[,i], xy$obs, method = "spearman")^2))
  print(paste0("Pearson: ", cor(xy[,i], xy$obs, method = "pearson")^2))
  print("####################")
}
summaries <- data.frame(rbind(do.call(rbind, nonspatial_summary), do.call(rbind, spatial_summary)))

predictors <- data.frame("predictor" = rep(row.names(nonspatial_summary[[1]]), length(ls(pattern = "_nonspatial_model$"))*2), 
           "upper" = summaries$X0.025quant,              
           "value" = summaries$X0.5quant,
           "lower" = summaries$X0.975quant,
           "spatial" = rep(c("non-spatial", "spatial"), each = nrow(spatial_summary[[1]]) * length(ls(pattern = "_nonspatial_model$"))),
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

spatial_dep <- segments_df
spatial_dep <- subset(spatial_dep, Site == "ZWE_MAT")
spatial_dep_non_spatial_m <- ncf::correlog(x = spatial_dep$long, y = spatial_dep$lat, z = spatial_dep$nonspatial_pred - spatial_dep$obs, increment = 5, resamp = 2, latlon = T)
spatial_dep_spatial_m <- ncf::correlog(x = spatial_dep$long, y = spatial_dep$lat, z = spatial_dep$spatial_pred - spatial_dep$obs, increment = 5, resamp = 2, latlon = T)


### FUNCTION ###
# Plot a Map of residuals, observations or predictions..
map_site <- function(df = segments_df, p = "spatial_pred", o = "obs", plot.what = "r", filter = NA, title = NA, interpolate = T){
  
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
  
  samples.df <- c(which(df$o == 1)[round(seq(1, length(which(df$o == 1)), length = 1000))],
                  which(df$o == 0)[round(seq(1, length(which(df$o == 1)), length = 1000))])
  #samples.df <- seq(nrow(df))
  
  if(interpolate){
    for(site in unique(df$Block)){
      test <- df[segments_df$Site == site, c("long", "lat", "input")]
      mba <- mba.surf(test, 300, 300, extend = FALSE)
      dimnames(mba$xyz.est$z) <- list(mba$xyz.est$x, mba$xyz.est$y)
      assign(paste0("ipblock_", site), reshape2::melt(mba$xyz.est$z, varnames = c('long', 'lat'), value.name = "input"))
    }
    
    df_map <- cbind(do.call(rbind, mget(ls(pattern = "ipblock_"))), "Block" = rep(unique(segments_df$Block), each = 300*300))
    
    gg_map <- 
    ggplot(data = df_map, aes(x = long, y = lat)) +
      geom_raster(aes(fill = input), interpolate = T, alpha = .5)+
      geom_point(data = df, aes(x = long, y = lat, col = input), fill = "black", shape = 21, size = .1) +
      scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(5, "Spectral")), guide =  "colourbar", na.value = "white")+
      scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(5, "Spectral")), guide =  "colourbar", na.value = "white")+
      theme_light()+
      guides(fill = guide_legend(title = title), color = F)+
      facet_wrap(~Block, scales = "free")
  }
  
  
  #print("create map")
  
  if(!interpolate){
    gg_map <- 
      ggplot(data = df, aes(x = long, y = lat, col = input)) +
      geom_point(size = .5)
    ggtitle(name)+
      theme_dark()+
      theme(legend.position = "top", legend.box = "vertical")+
      guides(color=guide_legend(title=name), size=guide_legend(title=name))+
      facet_wrap(~Site, scales = "free")+
      scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(5, "Spectral")), guide =  "colourbar")+
      theme_gray()
  }

  
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

segments_df$Block <- segments_df$Site
segments_df$WA <- segments_df$WA/1000
map_site(df = segments_df, o = "obs", p = "WA", plot.what = "p", title = "non-spatial model")

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

for(i in seq(10)) dev.off()
{
  pdf(file = paste0("Z:/residual_analysis/", selection, "_predictors.pdf"), onefile = T,paper = "a4", width = 8.27, height = 11.69)
  
  for(pred in names(segments_df[,3:21])[which(!names(segments_df[,3:21]) %in% c("ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
    map_site(df = segments_df, o = "obs", p = pred, plot.what = "p", title = pred)
  }
  
  dev.off() 
}

gg.xy <- cbind(xy_without_ZWE_MAT[, rev(strsplit(f, " +")[[1]][-seq(2, length(strsplit(f, " +")[[1]]), 2)])[-1]], xy[,c("obs", "nonspatial_pred", "spatial_pred")])
gg.xy <- reshape2::melt(gg.xy, id.vars = c("nonspatial_pred", "spatial_pred", "obs"))

  ggplot(gg.xy, aes(x = value, y = obs))+
    geom_smooth()+
    facet_wrap(~variable, scales = "free")



segm
# 
# bc_bbox <- make_bbox(lat = lat, lon = long, data = spatial_dep)
#  sc_surface <- get_map(location = bc_bbox, source = "google", maptype = "terrain")
#  sc_satellite <- get_map(location = bc_bbox, source = "google", maptype = "satellite")
#  
#  satellite <- 
#    ggmap(sc_satellite) +
#    ggtitle("Site: ZWE_MAT")
#  
#  terrain <- 
#    ggmap(sc_satellite) +
#    ggtitle("Site: ZWE_MAT"))

