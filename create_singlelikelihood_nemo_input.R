library(spdep)
library(INLA)



rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- rgdal::readOGR("Z:/GEC/segments.shp")
xy <- read.csv("Z:/modelling/yxtable_scaled_transformed.csv")[,-c(1)]

#create setup
model.type = "complex"
model.family = "binomial"
selection = "ZWE" #ISO3 Country Code or "all" if non selection should be made, Site Codes if sites are selected ("ZWE_MAT") c("BWA_NOR", "KEN_LAI", "KEN_TSV", "XWA_TBC", "ZWE_MAT", "ZWE_ZV", "ZWE_SELV")
selection.level = "Country" #"Country" or "Site"
nb_dist = 5000 #maximum distance considered as a neighbor [m]
xval = F #prepare data for cross validation
xval.type = "LSO" #LSO = leave some out, LOSO = leave one site out, KOSI = keep one site in
LSO.folds = 10 #number of folds
splines = c("WA", "VD", "AI") #transform these predictors to splines c("AI", "LD", "WA", "VD", "PI", "TD", "NA.")
n.knots = 2 #number of knots per spline per predictor
ridge = F #ridge regression
precision = 1e-4 #1e-4 is default
  

output_name <- paste0(paste0(selection, collapse = "_"), "_", model.type, "_", model.family, "_nonspatial_spatial_", ifelse(xval, paste0("xval_", xval.type, "_"), ""), nb_dist/1000, "km", ifelse(all(is.na(splines)), "", "_splines"), ifelse(ridge, "_ridge", paste0("_prec",precision)), ".RData")
source("Z:/GitRepo/function_file.R")

# define simple/complex formula and spatial model
{
  if(model.type == "complex"){
    f <- 
      "AI + HD + LD + NA. + NB + PA + PI + SM + TD + VD + WA + Site +
       SC:WA + SC:VD + TC:WA + CC:WA + TC:VD + LD:VD + PA:PI + Site:PI"
    for(spline in splines) f <- stringi::stri_replace_all_regex(f, paste0(" ",spline), paste0(spline, seq(n.knots), collapse = " + "))
  }
  
  if(model.type == "simple"){
    f <- "LD + NB + SC + VD + WA"
    for(spline in splines) f <- stringi::stri_replace_all_regex(f, paste0(" ",spline), paste0(spline, seq(n.knots), collapse = " + "))
    #f <- "LD + NB + SC + VD + WA + f(Site, model =\"iid\")"
  }
  
  if(model.family %in% c("binomial", "zeroinflated.binomial.0", "zeroinflated.binomial.1")){
    xy$obs <- as.numeric(xy$COUNT > 0) 
  }
  
  if(!model.family %in% c("binomial", "zeroinflated.binomial.0", "zeroinflated.binomial.1")){
    xy$obs <- xy$COUNT
  }
  
  f <- stringi::stri_replace_all_regex(f, "\n", "")
  #if(ridge) f <- paste0("-1 + X +  f(idx.Z, model=\"z\", Z=Z)")
  
  f.nonspatial <- formula(paste0("obs ~", f))
  f.spatial <- formula(paste0(paste(as.character(f.nonspatial)[c(2,1,3)], collapse = " "), "+ f(ID, model = \"besag\", graph = \"", tools::file_path_sans_ext(output_name),"_nb.txt\", param = c(0.1, 0.5))"))
 
  }

#create function to run different models
inla.fun <- function(x, formula, family, ridge, precision){
  #if(!ridge) 
  
  return(inla(formula, data = x, control.predictor = list(compute = T), family = family, control.fixed=list(mean=0,prec=precision)))
  
  # if(ridge){
  #   n <- nrow(x)
  #   X <- matrix(1,nr = n, nc = 1)
  #   Z <- x
  #   return(inla(formula, data = list(y= x$obs, idx.Z = 1:n, X=X), control.predictor = list(compute=TRUE), family = family))
  }
}

#create a subset of the data if selected
if(!"all" %in% selection) xy <- xy[xy[,selection.level] %in% selection, ]
for(i in c("ID", "CC", "HT", "PA", "Site", "Country", "Transect")){
  xy[, match(i, names(xy))] <- as.factor(as.character(xy[, match(i, names(xy))]))
}

#splines
if(!all(is.na(splines))){
  for(spline in splines){
    splined <- data.frame(spliner(xy[,match(spline, names(xy))], n.knots))
    for(i in ncol(splined))  splined[,i] <- scale(splined[,i])
    names(splined) <- paste0(spline, seq(n.knots))
    xy <- cbind(xy, splined)
  }
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


for(block in ls(pattern = "^xy_without_")) plot(is.na(get(block)$obs), pch = "|", cex = .5, xlim = c(0,1000), main = "block")

#create an output name, that defines setup characteristics, will be used in NEMO later on.
{
  save(list = ls()[!ls() %in% c("segments", "metric_segs")], file = paste0("Z:/NEMO_in/input_", output_name))
}

