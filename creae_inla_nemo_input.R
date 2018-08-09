library(spdep)
library(INLA)

segments <- rgdal::readOGR("Z:/GEC/segments_GEE.shp")

# f <- Y ~ 
#   mu.z + mu.y + 
#   AI + HD + LD + NB + PA + PI + SC + SL + SM + TC + TD + VD + WA + Site + 
#   WA^2 + PI^2 + 
#   SC:WA + SC:VD + TC:WA + TC:VD + PA:PI + Site:PI + LD:VD +
#   f(ID, model = "besag", graph = "dinla.txt", param = c(0.1, 0.5))+
#   f(HT, model = "iid")



f <- dnd ~ LD + NB + SC + VD + WA

inla.fun <- function(x, f){
  return(inla(f, data = x, control.predictor = list(compute = T), family = "binomial"))
}

m1 <- inla.fun(xy.adv, f)



#+ f(ID, model = "besag", graph = "dinla.txt", hyper = list(prec = list(prior = "loggamma", param = c(0.1, 1), initial = 0.5)))

xy <- read.csv("Z:/modelling/yxtable_scaled_transformed.csv")[,-c(1)]
table(xy$Country)
table(xy$Site)
selection = "none"
# selection = "BAW_NOR"
 selection = "ZWE"
#selection = "ETH"

if(selection != "none") xy <- subset(xy, Country == selection)

n.z <- as.numeric(rownames(xy))
n.y <- as.numeric(rownames(xy[which(xy$COUNT > 0),]))
xy.sel <- match(c(n.z, n.y), rownames(xy))

xy.adv <- xy[xy.sel,]
xy.adv[seq(length(n.z)), match(c("HT") , names(xy.adv))] <- NA
#xy.adv[(length(n.z) + 1) : length(xy.sel), match(c("NB", "TC") , names(xy.adv))] <- NA
for(i in c("ID", "CC", "HT", "PA", "Site", "Country", "Transect")) {
  xy.adv[, match(i, names(xy.adv))] <- as.factor(as.character(xy.adv[, match(i, names(xy.adv))]))
}

xy.adv <- xy.adv[seq(length(n.z)),]
xy.adv$dnd <- as.numeric(xy.adv$COUNT > 0) 
segments_new <- segments[match(as.character(xy.adv$ID), as.character(segments$ID)),]

for(i in seq(nrow(segments_new))) segments_new@polygons[[i]]@ID <- as.character(i)

#calculate neighborhood
metric_segs <- sp::spTransform(segments_new, "+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
cents <- geosphere::centroid(metric_segs)
dinla <- dnearneigh(as.matrix(cents), d1 = 1, d2 = 40000)


save(list = ls()[!ls() %in% c("segments", "metric_segs")], file = paste0("Z:/NEMO_in/input_", selection, ".RData"))

summary(xy.adv)


# 
# for(i in seq(length(dinla))){
#   if(i > length(n.z)){
#     input <- dinla[[i]][dinla[[i]] > length(n.z)]
#     if(length(input) == 0 | is.null(length(input))) input <- 0
#     dinla[[i]] <- input
#   }else{
#     input <- dinla[[i]][dinla[[i]] <= length(n.z)]
#     if(length(input) == 0 | is.null(length(input))) input <- 0
#     dinla[[i]] <- input
#   }
#   dinla[[i]] <- as.integer(dinla[[i]])
# } 
# 
# 


# 
# dinla[[1]]
# dinla[[length(n.z)]]
# dinla[[length(n.z)]]
# 
# for(i in sample(seq(nrow(segments_new)), 10, F)){
#   plot(segments_new[dinla[[i]],], main = paste0("neighbors of ", segments_new$ID[i]))
#   plot(segments_new[i,], col = "red", add = T)
#   points(cents[dinla[[i]],], pch = "+")
# }

# 
# 
# Y <- cbind(c(as.numeric(xy$COUNT > 0), rep(NA, sum(as.numeric(xy$COUNT > 0)))),
#            c(rep(NA, nrow(xy)), xy$COUNT[xy$COUNT > 0]))
# 
# ldat <- c(
#   list(
#     Y = Y,
#     mu.z = as.numeric(!is.na(Y[,2])),
#     mu.y = as.numeric(!is.na(Y[,1]))),
#   as.list(xy.adv))
# 



#rm(list = (ls()[ls() != "segments"]))
# 
# f
#n = sample(seq(nrow(segments)), 5000)

# metric_segs <- sp::spTransform(segments[c(n.1, n.2),], "+proj=lcc +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
# 
# dinla <- dnearneigh(as.matrix(geosphere::centroid(metric_segs)), d1 = 0.01, d2 = 6000)

# f <- COUNT ~ AI + COUNT + HD + LD + NB + PA + PI + SC + SL + SM + TC + TD + VD + WA +
#              WA^2 + PI^2 + Site + 
#              SC:WA + SC:VD + TC:WA + TC:VD + PA:PI + Site:PI + LD:VD 

# f(ID, model = "besag", graph = "dinla.txt", hyper = list(prec = list(prior = "loggamma", param = c(0.1, 1), initial = 0.5)))

