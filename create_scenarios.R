library(ggplot2)
library(reshape2)
setwd("Z:/IMAGE/")

segments <- readOGR("Z:/GEC/segments_GEE.shp")
load("image_extraction.RData")
xy <- read.csv("Z:/modelling/yxtable.csv")


for(img in c(7,8,12)){
  assign(paste0("y2100_", names(image.ext)[img]), image.ext[[img]]$y_2100)
} 

for(img in c(7,8,12)){
  assign(paste0("y2015_", names(image.ext)[img]), image.ext[[img]]$y_2015)
} 
names(image.ext)

scenario_2015 <- do.call(cbind, mget(ls(pattern = "y2015_")))
names_sc <- NULL
for(i in ls(pattern = "y2015_")) names_sc <- append(names_sc, names(get(i)))
names(scenario_2015)[c(1, 11)] <- c("aridity", "npp")
names(scenario_2015)[-c(1, 11)] <- names_sc[-c(1, 11)]
scenario_2015$ID <- segments$ID

scenario_2100 <- do.call(cbind, mget(ls(pattern = "y2100_")))
names_sc <- NULL
for(i in ls(pattern = "y2100_")) names_sc <- append(names_sc, names(get(i)))
names(scenario_2100)[c(1, 11)] <- c("aridity", "npp")
names(scenario_2100)[-c(1, 11)] <- names_sc[-c(1, 11)]
scenario_2100$ID <- segments$ID

write.csv(scenario_2100, "Z:/cc_model/IMAGE_2100.csv")
write.csv(scenario_2015, "Z:/cc_model/IMAGE_2015.csv")