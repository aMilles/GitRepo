library(spdep)
library(INLA)



rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- rgdal::readOGR("Z:/GEC/segments.shp")
xy <- read.csv("Z:/modelling/yxtable_scaled_transformed.csv")[,-c(1)]

#create setup
model.type = "simple"
model.family = "binomial"
selection = c("ZWE") #ISO3 Country Code or "all" if non selection should be made, Site Codes if sites are selected ("ZWE_MAT")
selection.level = "Country" #"Country" or "Site"
nb_dist = 5000 #maximum distance considered as a neighbor [m]
xval = F #prepare data for cross validation
xval.type = "LSO" #LSO = leave some out, LOSO = leave one site out, KOSI = keep one site in
LSO.folds = 10 #number of folds
splines <- c("WA", "VD")
n.knots <- 2

output_name <- paste0(paste0(selection, collapse = "_"), "_", model.type, "_", model.family, "_nonspatial_spatial_", ifelse(xval, paste0("xval_", xval.type, "_"), ""), nb_dist/1000, "km", ifelse(all(is.na(splines)), "", "_splines"), ".RData")
source("Z:/GitRepo/function_file.R")



# define simple/complex formula and spatial model
{
  if(model.type == "complex"){
    f <- 
      "AI + HD + LD + NB + PA + PI + SC + SL + SM + TC + TD + VD + WA + Site + 
    WA^2 + PI^2 + 
    SC:WA + SC:VD + TC:WA + TC:VD + PA:PI + Site:PI + LD:VD"
  }
  
  if(model.type == "simple"){
    f <- "LD + NB + SC + VD + WA + Site"
    for(spline in splines) f <- stringi::stri_replace_all_regex(f, spline, paste0(spline, seq(n.knots), collapse = " + "))
    #f <- "LD + NB + SC + VD + WA + f(Site, model =\"iid\")"
  }
  
  if(model.family %in% c("binomial", "zeroinflated.binomial.0", "zeroinflated.binomial.1")){
    xy$obs <- as.numeric(xy$COUNT > 0) 
  }
  
  if(!model.family %in% c("binomial", "zeroinflated.binomial.0", "zeroinflated.binomial.1")){
    xy$obs <- xy$COUNT
  }
  
  
  f.nonspatial <- formula(paste0("obs ~", f))
  #f.spatial <- formula(paste0(paste(as.character(f.nonspatial)[c(2,1,3)], collapse = " "), "+ f(ID, model = \"besag\", graph = ", output_name,"_nb.txt, param = c(0.1, 0.5))"))
  f.spatial <- formula(paste0(paste(as.character(f.nonspatial)[c(2,1,3)], collapse = " "), "+ f(ID, model = \"besag\", graph = \"dinla.txt\", param = c(0.1, 0.5))"))
  }

#create function to run different models
inla.fun <- function(x, formula, family){
  return(inla(formula, data = x, control.predictor = list(compute = T), family = family))
}

# r = inla(..., control.inla = list(
#   
#   correct = TRUE,
#   
#   correct.factor = <VALUE>))


#create a subset of the data if selected
if(!"all" %in% selection) xy <- xy[xy[,selection.level] %in% selection, ]
for(i in c("ID", "CC", "HT", "PA", "Site", "Country", "Transect")){
  xy[, match(i, names(xy))] <- as.factor(as.character(xy[, match(i, names(xy))]))
}

#splines
for(spline in splines){
  splined <- data.frame(spliner(xy[,match(spline, names(xy))], n.knots))
  names(splined) <- paste0(spline, seq(n.knots))
  xy <- cbind(xy, splined)
}


#create the neighborhood
{
  segments_new <- segments[match(as.character(xy$ID), as.character(segments$ID)),]
  metric_segs <- sp::spTransform(segments_new, "+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
  cents <- geosphere::centroid(metric_segs)
  dinla <- dnearneigh(as.matrix(cents), d1 = 1, d2 = nb_dist)
}

#create blocked xval in put
xy$Block <- NA
if(xval){
  if(xval.type %in% c("KOSI", "LOSO")){
    for(Site in as.character(unique(xy$Site))){
      xy_site <- xy
      if(xval.type == "LOSO") which.na <- xy_site$Site == Site
      if(xval.type == "KOSI") which.na <- xy_site$Site != Site
      xy_site[which.na, "obs"] <- NA
      xy$Block[which.na] <- Site
      assign(paste0("xy_without_", Site), xy_site)
    }
  }
  if(xval.type == "LSO"){
    zeroes <- which(xy$obs == 0)
    ones <- which(xy$obs > 0)
    which.zeroes <- rep(1:LSO.folds, length = length(zeroes))
    which.ones <- rep(1:LSO.folds, length = length(ones))
    
    for(i in seq(LSO.folds)){
      xy_fold <- xy
      xy_fold[c(ones[which.ones == i], zeroes[which.zeroes == i]), "obs"] <- NA
      xy$Block[c(ones[which.ones == i], zeroes[which.zeroes == i])] <- paste0("block", i)
      assign(paste0("xy_without_block", i), xy_fold)
    }
  }
} 


#WORK IN PROGRESS... 
if(effect){
  tf_sheet <- read.csv("Z:/modelling/transform_sheet.csv")
  scale_sheet <- read.csv("Z:/modelling/scale_sheet.csv")
  #scale_sheet$name <- tf_sheet$name
  xy_backup -> xy

  bt_mins <- us_mins <- mins <- apply(xy[,tf_sheet$name], 1, min)
  bt_maxs <- us_maxs <- maxs <- apply(xy[,tf_sheet$name], 1, max)
  
  for(i in length(mins)){
    us_mins[i] <- mins[i] * scale_sheet$scale[i] + scale_sheet$center[i]
    us_maxs[i] <- maxs[i] * scale_sheet$scale[i] + scale_sheet$center[i]
    # bt_mins[i] <- 
    # bt_maxs[i] <-
  }
#   seq()
}




for(block in ls(pattern = "^xy_without_")) plot(is.na(get(block)$obs), pch = "|", cex = .5, xlim = c(0,1000), main = "block")

#create an output name, that defines setup characteristics, will be used in NEMO later on.
{
  save(list = ls()[!ls() %in% c("segments", "metric_segs")], file = paste0("Z:/NEMO_in/input_", output_name))
}
output_name
