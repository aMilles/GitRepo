library(ggplot2)
library(gridExtra)
load("Z:/NEMO_out/marginals_all_complex_binomial_nonspatial_spatial_5km_splines.RData")

source("Z:/GitRepo/function_file.R")
for(i in seq(length(out)/4) - 1){
  selection <- seq(4) + i * 4
  print(selection)
  for(j in selection) assign(paste0("out", which(selection == j)), out[[j]])
  names(out3) <- "pred"
  effect.plot(quantiles.i = out1,
              quantiles.o = out2,
              df = out3,
              effect.of = out4,
              rug = F,
              with.lines = F)
}



f
