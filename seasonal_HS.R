library(ggplot2)
library(ranger)
library(rgdal)
#dropbox.file = "C:/Users/Alexander Milles/Dropbox/"
dropbox.file = "C:/Users/amilles/Dropbox/"
source("https://raw.githubusercontent.com/aMilles/GitRepo/master/function_file.R")

#read model output
load("Z:/NEMO_out/output_all_complex_binomial_nonspatial_spatial_6km_prec1e-04.RData")

segments <- readOGR(paste0(dropbox.file, "modelling/segments_GEE.shp"))
segments <- data.frame(segments, stringsAsFactors = F)
segments$time <- as.POSIXct(segments$time * 1000000000, origin = "1970-01-01")
backup -> segments
segments$ID <- as.character(segments$ID)

xy_season <- read.csv(paste0(dropbox.file, "/modelling/yxtable_season_scaled_transformed.csv"))
season_backup -> xy_season

#prepare seasonal dataset for prediction with model output
times <- as.Date(as.numeric(strftime(segments$time, format = "%j")), origin = "1970-01-01")
xy$time <- NA

months <- seq.Date(as.Date("1970-01-01"), as.Date("1970-12-01"), "month")
xy <- xy[match(as.character(xy_season$ID), as.character(xy$ID)),]

#replace predictors at the time of observation with seasonal predictors
xy[, c("VD", "WA", "TC", "SC")] <- xy_season[, c("VD", "WA", "TC", "SC")]
xy_season <- xy
xy_season$month <- rep(months, each = nrow(xy_season)/12)
xy_season$Site <- as.factor(as.character(xy_season$Site))
sites <- xy_season$Site
xy_season$Site <- "ZWE_MAT"
xy_season$Site <- factor(xy_season$Site, levels = levels(xy$Site))
xy_season$time <- times[match(xy_season$ID, segments$ID)]
time <- xy_season$time
xy_season$CC <- factor(xy_season$CC[1], levels = levels(xy$CC))

names(xy_season)

#create model formula without random effects
f.inla <- paste0(stringi::stri_replace_all_fixed(f,  "+ f(Site ,model=\"iid\") + f(PA ,model=\"iid\")", ""))
m <- model.matrix(as.formula(paste0("~", f.inla)), xy_season)
dist <- as.matrix(sample.fixed.params(spatial_model, 10))
out <- invlogit(m %*% t(dist))
HS <- t(apply(data.frame(out), 1, quantile, probs = c(0.025, 0.5, 0.975)))
HS <- data.frame(HS)
names(HS) <- c("lower", "mid", "upper")


gg_season <- cbind(xy_season[,c("Country", "Site", "ID", "month", "time")], HS)
sites -> gg_season$Site
#gg_season <- reshape2::melt(gg_season, id.vars = c("Country", "Site", "ID", "month"))
#gg_season$ranger <- test$predictions[,2]

gg_season.agg <- aggregate(gg_season$mid, by = list(gg_season$Site, gg_season$month, gg_season$Country), FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
gg_season.agg <- cbind(gg_season.agg[,-4], data.frame(gg_season.agg$x))
names(gg_season.agg) <- c("Site", "month", "Country", "lower", "mid", "upper")
gg_season.agg$ID <- NA
gg_season.agg$time <- as.Date(0, origin = "1970-01-01")

#calculate accumulated habitat suitability (3-months)
gg_season.agg$mm <- NA
for(Site in unique(gg_season.agg$Site)){
  ss <- gg_season.agg[gg_season.agg$Site == Site,]
  ss$mm <- (pracma::movavg(c(ss$upper[11:12], ss$upper), 3, "s")[-c(1,2)])
  gg_season.agg[gg_season.agg$Site == Site,] <- ss
}
SC_agg <- aggregate(SC ~ Site + month, xy_season, mean)
names(SC_agg) <- c("Site", "month", "SC")
gg_season.agg <- merge(gg_season.agg, SC_agg)
row.nr <- nrow(gg_season.agg)
gg_season$mm <- NA
gg_season$SC <- NA
gg <- rbind(gg_season.agg, gg_season)
get.mins <- aggregate(mm ~ Site, data = gg[seq(row.nr),], which.min)
get.mins$month <- months[get.mins$mm]

n.sites<- table(unlist(strsplit(as.character(unique(gg_season$Site)), "_"))[seq(1,length(unique(gg_season$Site))*2, 2)])

colors <- NULL
for(i in n.sites) colors <- append(colors, viridis::viridis(i, begin = .1, end = .8))

#create plot of seasonal variation of habitat suitability
hs_season <- ggplot(gg)+
  geom_ribbon(data = gg[seq(row.nr),], aes(x = month, ymin = lower, ymax = upper, group = Site), fill = "red", alpha = .5)+
  #scale_fill_manual(values = colors)+
  scale_color_manual(values = c("black", "turquoise"))+
  geom_line(data = gg[seq(row.nr),], aes(x = month, y = mid, group = Site, col = "median HS"), size = 1)+  geom_line(data = gg[seq(row.nr),], aes(x = month, y = scales::rescale(SC, to  = c(0, 0.09)), group = Site, col = "4-month precipitation"), size = 1)+
  #geom_point(data = gg[seq(row.nr),], aes(x = month, y = mid, group = Site))+
  geom_rug(data = gg[seq(nrow(xy_season)/12)+row.nr,], aes(x = time), alpha = 0.05, color = "navyblue")+
  facet_wrap(~Site)+
  ylab("habitat suitability / scaled 4-month precipitation")+
  theme_bw()+
  geom_vline(aes(xintercept = month), data = get.mins, linetype = "dashed", color = "black")+
  scale_x_date(date_breaks = "4 months", date_labels = "%b")+
  scale_color_manual(values = c("skyblue", "black"))+
  guides(color = guide_legend(title = ""))+
  theme(text = element_text(size = 14), legend.position = "top")
hs_season
pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/HS_season.pdf"))
hs_season
dev.off()

xy <- read.csv(paste0(dropbox.file, "modelling/yxtable.csv"), row.names = 1)

xy_seasonbackup -> xy_season
xy_season_us <- read.csv(paste0(dropbox.file, "/modelling/yxtable_season.csv"))
xy_season_us$month <- as.numeric(factor(as.character(xy_season_us$month), unique(xy_season_us$month)))
names(xy_season_us)[c(4, 10, 13)] <- c("AD","SS","AR")
xy_season_us$AR <- xy$AR
xy_season_us$HD <- xy$HD
for(i in seq(16)){
  Site = get.mins$Site[i]
  month = as.numeric(strftime(get.mins$month, "%m"))[i]
  assign(paste0("var_", Site), apply(xy_season_us[xy_season_us$Site == Site & xy_season_us$month == month, c("AD", "AR", "HD", "LD", "NB", "SS", "TD", "VD", "WA")], 2, mean, na.rm = T))
}

table <- do.call(rbind, mget(ls(pattern = "var_")))
rownames(table) <- apply(do.call(rbind, strsplit(rownames(table), "_"))[,2:3], 1, paste, collapse = "_")
table <- data.frame(table, row.names = rownames(table))
table$VD <- table$VD / 10000
table$WA <- table$WA / 1000
signif(table, 2)
psych::df2latex(table)

paste0(as.numeric(as.factor(get.mins$Site, gg$Site[seq(row.nr)])))
                  
                  
##############################
# SECTION BELOW WAS NOT USED #
##############################

if(F){
  get.mins
  
  #Select month of minimum suitability
  xy_season$Site <- sites
  summary(xy_season)
  xy_minHS <- xy_season[seq(nrow(xy_season)/12),]
  xy_minHS[, c("VD", "WA", "SC", "TC")] <- NA
  
  xy_season$Site <- sites
  
  for(i in seq(length(get.mins$Site))){
    rows <- which(xy_season$month == get.mins$month[i] & as.character(xy_season$Site) == as.character(get.mins$Site[i]))
    these.id <- xy$ID[match(as.character(xy_season$ID[rows]), as.character(xy$ID))]
    #if(!all(as.character(xy_season$ID[rows]) == as.character(xy_minHS$ID[xy_minHS$Site == get.mins$Site[i]]))) warning('Something is wrong!!!!')
    
    xy_minHS[match(these.id, xy_minHS$ID), ] <- xy_season[rows, ]
  }
  write.csv(xy_minHS, paste0(dropbox.file, "modelling/xytable_minHS.csv"))
  warnings()
  summary(xy_minHS)
  get.maxs <- aggregate(mid ~ Site, data = gg[seq(row.nr),], which.max)
  get.maxs$month <- months[get.maxs$mid]
  
  xy_maxHS <- xy
  xy_maxHS[, c("VD", "WA", "SC", "TC")] <- NA
  xy_maxHS$Block <- NULL
  xy_season$Site <- sites
  for(i in seq(length(get.maxs$Site))){
    rows <- which(xy_season$month == get.maxs$month[i] & as.character(xy_season$Site) == as.character(get.maxs$Site[i]))
    these.id <- xy$ID[match(as.character(xy_season$ID[rows]), as.character(xy$ID))]
    #if(!all(as.character(xy_season$ID[rows]) == as.character(xy_maxHS$ID[xy_maxHS$Site == get.maxs$Site[i]]))) warning('Something is wrong!!!!')
    
    xy_maxHS[match(these.id, xy_maxHS$ID), ] <- xy_season[rows, -1]
  }
  summary(xy_maxHS)
  write.csv(xy_maxHS, paste0(dropbox.file, "modelling/xytable_maxHS.csv"))
  
  
  
  
  subset <- gg_season[1:2,]
  
  for(i in seq(nrow(xy_minHS))){
    x <- gg_season[gg_season$Site == xy_minHS$Site[i] & gg_season$month == xy_minHS$month[i],]
    subset <- rbind(subset, x)
  }
  
  
  count <- aggregate(xy_season$obs, by = list(xy_season$Site), FUN = sum)
  count$x <- count$x/table(xy_season$Site)
  names(count) <- c("Site", "density")
  subset <- merge(gg[seq(row.nr),], count)
  subset$density <- as.vector(subset$density)
  best.sites <- c("BWA_NOR", "KEN_TSV", "ZWE_MAT", "ZWE_SELV", "ZWE_ZV", "KEN_LAI")
  subset$good <- subset$Site %in% best.sites
  summary(subset$mm)
  agg <- aggregate(subset$mm, by = list(subset$Site), FUN = median)
  agg$density <- as.vector(count$density)
  names(agg)[1:2] <- c("Site", "mid")
  agg$good <- agg$Site %in% best.sites
  #subset$Site <- as.character(subset$Site)
  
  pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/minHS_bottleneck.pdf"))
  ggplot(agg, aes(y = density, x = mid))+
    geom_point(size = 5, shape = "+")+
    geom_text(aes(label = Site), size =3, nudge_y = .02, check_overlap = T)+
    theme_bw()+
    theme(legend.position = "bottom", text = element_text(size = 14))+
    xlab("median minimum habitat suitability [%]")+
    ylab("detection rate per subunit")+
    xlim(0, 0.02)
  dev.off()
  
  
  agg.estimate <- agg
  
  agg.estimate$density[agg$Site == "KEN_LAM"] <- agg$density[agg$Site == "KEN_LAM"] * 700
  agg.estimate$density[agg$Site == "ZWE_SEB"] <- agg$density[agg$Site == "ZWE_SEB"] * 1/.26
  agg.estimate$density[agg$Site == "TCD_ZAK"] <- agg$density[agg$Site == "TCD_ZAK"] * 8 #https://news.nationalgeographic.com/2017/01/wildlife-watch-chad-zakouma-elephants-poaching/
  agg.estimate$density[agg$Site == "ZWE_MAT"] <- agg$density[agg$Site == "ZWE_MAT"] / 3
  agg.estimate$density[agg$Site == "COD_GAR"] <- agg$density[agg$Site == "COD_GAR"] * 227/12 # 22,700 in the mid-'70s to less than 1,200 today. https://www.gq.com/story/inside-the-ivory-wars
  agg.estimate$density[agg$Site == "COD_VIR"] <- agg$density[agg$Site == "COD_VIR"] * 5000/300 #https://www.reuters.com/article/us-congo-democratic-elephants/ivory-poachers-decimate-congo-elephant-population-idUSLM40286220080822 
  ggplot(agg.estimate, aes(y = density, x = mid, col = good))+
    geom_point()+
    geom_text(aes(label = Site), size =3, nudge_y = .02)+
    theme_bw()+
    theme(legend.position = "bottom")+
    guides(col = guide_legend(title = "equlibrium assumed"))+
    xlab("median minimum habitat suitability")+
    ylab("detection rate per subunit")
  
  
  aggregate(all.xy$COUNT, by = list(all.xy$Site), FUN = mean)
  geom_point(shape = "|", alpha = .01)
  
  
  all.xy <- read.csv(paste0(dropbox.file, "/modelling/yxtable.csv"))  
  f
  gg.all <- (all.xy[,c("AD", "HD", "LD", "SS", "AR", "TD", "VD", "WA", "TV", "TC","SC", "Site")])
  
  
  pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/", output_name, "site_predictors_all.pdf"), onefile = T)
  for(i in c("AD", "HD", "LD", "SS", "AR", "TD", "VD", "WA", "TV", "TC","SC")){
    print(ggplot(gg.all, aes(y = gg.all[,i], x = Site))+
            geom_boxplot()+
            coord_flip()+
            ylab(i))
  }
  dev.off()
  
  
  plot(minHS)
  hs_season
  pdf(paste0(dropbox.file, "Master/Umweltwissenschaften/Masterarbeit/figures/", output_name, "_HS_season.pdf"))
  hs_season
  dev.off()
}