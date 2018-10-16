library(rgdal)
library(zoo)
library(stringi)
library(pscl)
library(MASS)
library(corrplot)

#this is the path to directory with a "modelling"-folder
dropbox.file = "C:/Users/amilles/Dropbox/"

#read csv.files
xy_season <- read.csv(paste0(dropbox.file, "modelling/yxtable_season.csv"), row.names = 1)
names(xy)
xy <- read.csv(paste0(dropbox.file, "modelling/yxtable.csv"), row.names = 1)

xy_seasonbackup -> xy_season
names(xy_season)[match(c("AI", "NA.", "PI"), names(xy_season))] <- c("AD", "SS", "AR")
xy_season$AR <- xy$AR
xy_season$HD <- xy$HD
xy_season$VD <- xy_season$VD / 10000
xy_season$WA <- xy_season$WA / 1000
xy_season$TC <- xy_season$TC - 273.15
summary(xy_season$AR)
xy_season <- na.omit(xy_season) 

#combine seasons and xy to transform and scale them in one run.
all.preds <- rbind(xy, xy_season[,names(xy)])
transformed <- all.preds

#apply box cox transformation, store the transformation values for later backtransformation
bc_values <- NULL
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/", as.Date(Sys.time()),"_transform.pdf"), onefile = T, paper = "a4")
for(col in names(all.preds)[which(!names(all.preds) %in% c("REPS", "ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
  cov <- all.preds[,col]
  if(min(cov) <= 0) cov <- cov + abs(min(cov))
  cov <- cov + (min(cov[cov > 0]) / 2)
  par(mfrow = c(1,1))
  #hist(cov, main = col)
  fm <- lm(cov ~ 1)
  
  lambdas <- c(0, tan(seq(-1.2,1.2, length = 200))) # length can be set even higher, increase range if bc_values are at max/min values
  bc <- boxcox(fm, lambda = lambdas, plotit = F)
  bc_value <- bc$x[which.max(bc$y)]
  if(bc_value == 0){
    transformed[,col] <- log(cov)
  }else{
    transformed[,col] <- cov^bc_value
  }
  print(bc_value)
  bc_values <- append(bc_values, bc_value)
  print(bc_values)
  
  par(mfrow = c(2,2))
  qqnorm(transformed[,col], main = paste0(col, " transformed"))
  qqline(transformed[,col])
  hist(transformed[,col])
  
  qqnorm(all.preds[,col], main = paste0(col, " not transformed"))
  qqline(all.preds[,col])
  hist(all.preds[,col])
}
dev.off()
#save and check transformation - plots should look the same for original values and backtransformation

transform_sheet <- 
  data.frame(name = names(all.preds)[which(!names(all.preds) %in% c("REPS", "ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))],
             value = bc_values, shift = 0)

for(pred in names(all.preds)[which(!names(all.preds) %in% c("REPS", "ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
  orig <- all.preds[, pred]
  
  bt_value = transform_sheet$value[which(transform_sheet$name == pred)]
  if(bt_value == 0) bt <- exp(transformed[,pred])
  if(bt_value != 0) bt <- transformed[,pred]^(1/bt_value)
  
  shift <- median(orig - bt)
  transform_sheet$shift[which(transform_sheet$name == pred)] <- shift
  par(mfrow = c(1,2))
  hist(orig, main = paste0("original ", pred))
  hist((bt + shift), main = paste0("backtransformed ", pred))
} 


backtransform <- c("SC", "SL", "SM", "TC", "VD")
transformed[,backtransform] <- all.preds[,backtransform]

transform_sheet$value[transform_sheet$name %in% backtransform] <- 1
transform_sheet$shift[transform_sheet$name %in% backtransform] <- 0


for(pred in names(all.preds)[which(!names(all.preds) %in% c("REPS", "ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
  orig <- all.preds[, pred]
  
  bt_value = transform_sheet$value[which(transform_sheet$name == pred)]
  if(bt_value == 0) bt <- exp(transformed[,pred])
  if(bt_value != 0) bt <- transformed[,pred]^(1/bt_value)
  
  shift <- median(orig - bt)
  transform_sheet$shift[which(transform_sheet$name == pred)] <- shift
  par(mfrow = c(1,2))
  hist(orig, main = paste0("original ", pred))
  hist((bt + shift), main = paste0("backtransformed ", pred))
} 



write.csv(transform_sheet, paste0(dropbox.file, "modelling/transform_sheet.csv"))
write.csv(transformed[seq(nrow(xy)), ], paste0(dropbox.file, "/modelling/yxtable_transformed.csv"))
write.csv(transformed[seq(nrow(xy)+1, nrow(transformed)), ], paste0(dropbox.file, "/modelling/yxtable_season_transformed.csv"))

#scale non-factor predictors
scale_transformed <- transformed

centers <- NULL
scales <- NULL
names <- NULL
for(col in which(!names(transformed) %in% c("REPS", "ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))){
  scale_transformed[,col] <- scale(transformed[,col])
  centers <- append(centers, attr(scale_transformed[,col], "scaled:center"))
  scales <- append(scales, attr(scale_transformed[,col], "scaled:scale"))
  names <- append(names, col)
} 

#save and check scaling parameters - plots should look the same

scale_sheet <- data.frame("center" = centers, "scale" = scales, "name" = names(all.preds)[which(!names(all.preds) %in% c("REPS", "ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))])

for(pred in names(all.preds)[which(!names(all.preds) %in% c("REPS", "ID", "CC", "COUNT", "HT", "PA", "Site", "Transect", "Country"))]){
  orig <- transformed[, pred]
  scaled <- scale_transformed[, pred]
  
  center_value = scale_sheet$center[which(scale_sheet$name == pred)]
  scale_value = scale_sheet$scale[which(scale_sheet$name == pred)]
  
  
  us <- scaled * scale_value + center_value
  
  hist(orig, main = paste0("original ", pred))
  hist(us, main = paste0("unscaled ", pred))
}

write.csv(scale_sheet, paste0(dropbox.file, "modelling/scale_sheet.csv"))

#check how the scaled, transformed values are distributed in comparison to the original values
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/", as.Date(Sys.time()),"_transform_scale_hist.pdf"), onefile = T, paper = "a4")
par(mfrow = c(2,2))
for(i in which(!names(transformed) %in% c("REPS", "ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))){
  hist(scale_transformed[,i], main = names(scale_transformed)[i], breaks = 1000)
  hist(all.preds[,i], main = names(scale_transformed)[i], breaks = 1000)
}
dev.off()

scale_transformed$CC <- as.factor(scale_transformed$CC)


source("Z:/GitRepo/modified_cormat.R")

#check for correlation
#correlation between values before and after transforamtion + scaling, HD old and HD new are highly NEGATIVELY correlated (just to keep it in mind)

rquery.cormat.spearman(cbind(all.preds[seq(nrow(xy)), which(!names(scale_transformed) %in% c("REPS", "ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))],
                             scale_transformed[seq(nrow(xy)), which(!names(scale_transformed) %in% c("REPS", "ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))]))$r


#correlation between transformed + scaled values
spearman <- rquery.cormat.spearman(scale_transformed[seq(nrow(xy)), which(!names(scale_transformed) %in% c("REPS", "ID", "CC", "HT", "PA", "COUNT", "Country", "Transect", "Site"))])$r
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/corplot.pdf"), width = 4, height = 4)
spearman
dev.off()
spearman <- data.frame(apply(spearman, 1, function(x) as.numeric(as.character(x))))
cornames <- names(spearman) 
spearman$names <- cornames

cors <- na.omit(reshape2::melt(spearman))
cors[which(abs(cors$value) > 0.5 & cors$names != as.character(cors$variable)), ]


write.csv(scale_transformed[seq(nrow(xy)), ], paste0(dropbox.file, "/modelling/yxtable_scaled_transformed.csv"))
write.csv(scale_transformed[seq(nrow(xy)+1, nrow(transformed)), ], paste0(dropbox.file, "/modelling/yxtable_season_scaled_transformed.csv"))
names(scale_transformed)
summary(scale_transformed)