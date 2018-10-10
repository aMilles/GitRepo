library(INLA)
library(ggplot2)
library(gridExtra)
library(parallel)
library(foreach)
library(doSNOW)
library(doParallel)
library(gridExtra)
library(ranger)

#read Model output from nemo
load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_6km_prec1e-04.RData")
load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_xval_LOSO_6km_prec1e-04.RData")

rm(list = ls(pattern = "nonspatial"))

model <- ls(pattern = "spatial_model")[1]
out <- vector("list", length = length(ls(pattern = "spatial_model")))
for(i in seq(length(ls(pattern = "spatial_model")))){
  model <- ls(pattern = "spatial_model")[i]
  sum <- summary(get(model))$fixed[,3:5]
  site.rows <- c(grep("Site", rownames(sum)), grep("PA", rownames(sum)))
  if(length(site.rows) == 0) site.rows <- nrow(sum) + 1
  out[[i]] <- sum[-c(1, site.rows), ]
}



coefs <- data.frame(do.call(rbind, out))
coefs$predictor <- row.names(out[[1]])
coefs$model <- rep(ls(pattern = "spatial_model"), each = nrow(out[[1]]))
names(coefs)[1:3] <- c("low", "mid", "up")
coefs$midrel <- (coefs[1:nrow(out[[1]]),2] - coefs[,2]) 
coefs$region <- region <- apply(do.call(rbind, strsplit(coefs$model, "_"))[,c(3,4)], 1, paste, collapse = "_")
coefs$region[coefs$region == "spatial"] <- NA


tf_sheet <- read.csv("C:/Users/amilles/Dropbox/modelling/transform_sheet.csv")
reverse <- tf_sheet$name[tf_sheet$value < 0]
reverse
for(pred in reverse){
  rows <- (stringi::stri_count_regex(as.character(coefs$predictor), pred)) > 0
  coefs[rows, c(1:3)] <- coefs[rows, c(1:3)] * -1
  
}

coefs$region[coefs$region == "spatial_model"] <- NA 
coefs$region[!coefs$region %in% c(NA, "BWA_NOR", "COD_GAR", "KEN_TSV", "ZWE_SELV", "XWA_TBC")] <- "other (n = 11)"
coefs$region <- factor(coefs$region, levels = c("BWA_NOR", "COD_GAR", "KEN_TSV", "ZWE_SELV", "XWA_TBC", "other (n = 11)"))
unique(coefs$region)

colors <- RColorBrewer::brewer.pal(6, "Accent")
colors[4] <- "gray70"
colors
levels(coefs$region)
xval <- ggplot(coefs[-(1:nrow(out[[1]])),], aes(x = predictor, y = mid, fill = region))+
  geom_abline(aes(intercept = 0, slope = 0), col = "navyblue", size =2)+
  geom_crossbar(aes(ymin = low, ymax = up, y = mid), color = NA, alpha = 1)+
  geom_point(shape = "|", size = 4.5, color = "white")+
  geom_crossbar(data = coefs[1:nrow(out[[1]]),], aes(y = mid, ymin = low, ymax = up), color = "black", alpha = 0, size = 0.5)+
  theme_bw()+
  coord_flip()+
  scale_fill_manual(values = colors)+
  #scale_fill_brewer(type = "qual", palette =  "Set2")+
  theme(legend.position = "none", text = element_text(size = 14), panel.grid.major.y = element_line(color = c("gray60", "white")))+
  ylab("fixed effects, 95-% interval")
xval <- xval+theme(legend.position = "top")+guides(fill = guide_legend(title = "site"))
xval
pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/xval.pdf", width = 8, height = 5)
xval
dev.off()



