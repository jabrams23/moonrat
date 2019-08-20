library(raster)
library(rgdal)
library(rgeos)

wd <- ""
setwd(wd)

# Loading camera trap locations
#camera_traps <- readOGR(wd, "VN_CTs")
#rasterize(camera_traps, cities, filename = "C:/Users/tilker/Dropbox (ScreenForBio)/andrew+jesse/Andrew_covariates/input_files/shapefiles/camera_traps_r.tif", options = "COMPRESS = IZW", progress = "text")

grid_preds <- raster("distance_landcover_plantation.tif")

#############################################
# Extracting distances from Da Nang and Hue #
#############################################

roads_v <- shapefile("Roads_Deramakot-Tangkulap-Kuamut_latest_03072017.shp")
roads_v = as(roads_v, "SpatialPointsDataFrame")
roads_v <- spTransform(roads_v, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

road_distance <- distanceFromPoints(grid_preds, roads_v, progress = "text")

road_distance_clipped <- overlay(road_distance, grid_preds, fun = function(x, y) {
  x[is.na(y[])] <- NA
  return(x)
})


writeRaster(road_distance_clipped, "roads_distance", format = "GTiff")

road_distance_clipped <- raster("roads_distance.tif")


png(filename = "road_distance_OVERLAID.png", width=6*600, height=6*485, 
    res=600, bg="white")
plot(road_distance_clipped)
plot(roads_v, add = T, col = "black", cex = 1, pch = 19)
title(main="Distance from cities")
dev.off()
