install.packages('ggplot2')
library(ggplot2)

splines <- c("AI", "LD", "PI", "WA", "VD", "NA.")

xy_radar <- read.csv("Z:/modelling/yxtable.csv")  
rm(list = ls(pattern = "radar_"))  
for(Site in c("BWA_NOR", "KEN_LAI", "KEN_TSV", "XWA_TBC", "ZWE_MAT", "ZWE_ZV", "ZWE_SELV")) assign(paste0("radar_", Site), ggiraphExtra::ggRadar(data = xy_radar[xy_radar$Site == Site, c(splines, "HT", "Site")], aes(color = HT), scales = "free", rescale = T))

do.call(gridExtra::grid.arrange, mget(ls(pattern = "radar_")))


ggiraphExtra::ggRadar(data = xy_radar[, c(splines, "HT", "Site")], aes(color = HT), scales = "free", rescale = T)

xy_boxplot <- reshape2::melt(xy_radar[xy_radar$Site %in% c("BWA_NOR", "KEN_LAI", "KEN_TSV", "XWA_TBC", "ZWE_MAT", "ZWE_ZV", "ZWE_SELV"), c(splines, "HT", "Site")], id.vars = c("HT", "Site"))
ggplot(xy_boxplot, aes(x = HT, y = value))+
  geom_boxplot(outlier.alpha = 0)+
  facet_wrap(~variable, scales = "free")


agg <- aggregate(xy_boxplot$value, by = list(xy_boxplot$HT, xy_boxplot$variable, xy_boxplot$Site), FUN = function(x) quantile(x, probs = c(0.5)))
#for(i in seq(4)) agg[,i] <- as.factor(as.character(agg[,i]))
agg <- cbind(agg[,1:3], data.frame(agg[,4]))
agg <- reshape2::melt(agg, id.vars = c("Group.1", "Group.2", "Group.3"))

pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/predictor_boxplot.pdf", width = 7.27, height = 3.5)
ggplot(agg, aes(x = Group.1, y = value))+
  geom_boxplot()+
  geom_point(aes(shape = Group.3), size = 2)+
  facet_wrap(~Group.2, scales = "free")+
  theme_classic()+
  xlab("Herd Type")+
  ylab("Predictor Value")+
  scale_shape_manual(values = seq(8)+1) +
  guides(shape = guide_legend(title = "Site"))+
  theme(text = element_text(size=14))
dev.off()  
