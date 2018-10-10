library(rgdal)
library(ggplot2)
library(ggmap)
library(MBA)
library(geosphere)
library(sp)
library(raster)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(raster)
library(rasterVis)
library(rgdal)
rm(list = ls()[which(ls() != "segments")])
if(!"segments" %in% ls()) segments <- readOGR("Z:/GEC/segments_GEE.shp")
if(!"africa" %in% ls()) africa <- readOGR("Z:/africa/africa_countries.shp")

centroids <- rgeos::gCentroid(segments, byid = T)
coords <- coordinates(centroids)
plot(centroids, pch = ".", cex = .01)
lines(coords[chull(coords),], col = "red")
chull(centroids)

sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords[chull(coords),])), ID=1)))
sp_poly_df <- SpatialPolygonsDataFrame(sp_poly, data=data.frame(ID=1))
extent(sp_poly)+1
r <- raster(extent(sp_poly) + 1)
r@crs <- crs(segments)
res(r) <- 0.1
r[] <- seq(ncell(r))
plot(r)
points(centroids)
library(velox)
segments
vr <- velox::velox(r)
what_pixel <- vr$extract(segments,small = T)

pixel_ext <- data.frame(pixel.id =  do.call(c, what_pixel), count = 1)
pixel_ext_agg <- aggregate(count ~ pixel.id, pixel_ext, sum)

r[r < max(pixel_ext_agg$count)] <- 0

for(i in seq(length(pixel_ext_agg$pixel.id))){
  id = pixel_ext_agg$pixel.id[i]
  print(i)
  r[r == id] <- pixel_ext_agg$count[i]
}
r[r > max(pixel_ext_agg$count)] <- 0
plot(r)
ggplot(r)
r2 <- r
r2[r2 == 0] <- NA
test_spdf <- as(r2, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

map <- ggplot() +  
  coord_equal()+
  geom_polygon(data=africa, aes(x=long, y=lat, group=group), 
             fill="grey90", color="grey80", size=0.1, alpha = .5)+
  geom_tile(data=test_df, aes(x=x, y=y, fill=value/100)) + 
  scale_fill_gradientn(colors =  viridis::viridis(10, end = .8)) +
  theme_bw()+
  xlab("longitude [°E]")+
  ylab("latitude [°N]")+
  xlim(c(-20, 52))+
  ylim(c(-35, 15))+
  theme(text = element_text(size = 11), legend.position = "bottom")+
  guides(fill = guide_colorbar(title = "sampling effort [subunits/km²]", barwidth = 10))

pdf("C:/Users/amilles/Dropbox/Master/Umweltwissenschaften/Masterarbeit/figures/sampling_effort_map.pdf")  
map  
dev.off()

vr <- velox::velox(r)
what_pixel <- vr$extract(segments,small = T, fun = mean)
sampling.effort <- data.frame(ID = segments$ID, sampling.effort = what_pixel[,1])

write.csv(sampling.effort, "C:/Users/amilles/Dropbox/modelling/sampling_effort.csv")
summary(sampling.effort$sampling.effort)

