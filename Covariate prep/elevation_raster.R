
library(raster)
library(rgeos)
library(RQGIS)
library(rgdal)
library(vegan)
library(igraph)

setwd("~/Dropbox (ScreenForBio)/Projects/moonrat")

#Sabah
wd_landcov_sabah <- "~/Dropbox (ScreenForBio)/Projects/moonrat"
filename_sabah <- "Sabah_25_08_2017.tif"
filepath_landcov_sabah <- file.path(wd_landcov_sabah, filename_sabah)
landcov <- raster(filepath_landcov_sabah)

TFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/sg pinangahfr east with Tangkulap.shp")
DFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/deramakot_fr.shp")
KFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/kuamut fr_west.shp")

TFR <- spTransform(TFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
DFR <- spTransform(DFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
KFR <- spTransform(KFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

combined <- TFR+DFR+KFR
combined<- spTransform(combined, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

r.dem30 <- raster("~/Dropbox (ScreenForBio)/Projects/moonrat/n05_e117_1arc_v3.tif")    # change to the correct file (see email 2017-02-02)

fs <- raster('~/Dropbox (ScreenForBio)/Projects/moonrat/Prediction_rasters/FS_agg_100m_cell_clipped.tif')

elev <- projectRaster(r.dem30, fs, res=100, crs="+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
r5 <- crop(elev, extent(combined))
elev_clipped <- mask(r5, combined)
writeRaster(elev_clipped, filename = "elev_100mcell_clipped.tif", format = "GTiff", datatype = "FLT4S")
