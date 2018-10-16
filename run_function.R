library(spdep)
library(INLA)

setwd("/home/fr/fr_fr/fr_am595/data/run_this")

#always load the new input file
print(paste0("load ", list.files(pattern = "input_")[which.max(base::file.info(list.files(pattern = "input_"))$ctime)]))
load(list.files(pattern = "input_")[which.max(base::file.info(list.files(pattern = "input_"))$ctime)])

summary <- data.frame(start = Sys.time(), end = NA, name = output_name)

nb2INLA(paste0(tools::file_path_sans_ext(output_name),"_nb.txt"), dinla)
print(paste0("created nb: ", paste0(tools::file_path_sans_ext(output_name),"_nb.txt")))
#nb2INLA("dinla.txt", dinla)

if(xval){
  for(site_input in ls(pattern = "^xy_without")){
    df <- get(site_input)
    df$ID <- 1:nrow(df)
    assign(paste0(site_input, "_spatial_model"), inla.fun(df, f.spatial, family = model.family, precision = precision))
    assign(paste0(site_input, "_nonspatial_model"), inla.fun(df, f.nonspatial, family = model.family, precision = precision))
  }
}

if(!xval){
  xy.model <- xy
  xy.model$ID <- 1:nrow(xy)
  nonspatial_model <- inla.fun(xy.model, f.nonspatial, family = model.family, precision = precision)
  spatial_model <- inla.fun(xy.model, f.spatial, family = model.family, precision = precision)
}

summary$end <- Sys.time()

all.summaries <- read.csv("summary.csv", stringsAsFactors = F)[,-1]
all.summaries <- rbind(all.summaries, summary)
write.csv(all.summaries, "summary.csv")

save.image(paste0("output_", output_name))