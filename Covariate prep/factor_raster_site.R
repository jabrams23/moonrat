
library(raster)
library(rgeos)
library(RQGIS)
library(rgdal)
library(vegan)
library(igraph)

setwd("~/Dropbox (ScreenForBio)/Projects/moonrat/Prediction_rasters")

TFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/sg pinangahfr east with Tangkulap.shp")
DFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/deramakot_fr.shp")
KFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/kuamut fr_west.shp")

TFR <- spTransform(TFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
DFR <- spTransform(DFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
KFR <- spTransform(KFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

fs <- raster('FS_agg_100m_cell_clipped.tif')

r <- fs #raster(nrow=1110, ncol=1112)
res(r)<-100
r[DFR] <- 1
r[KFR] <- 2
r[TFR] <- 3
r <- ratify(r)

rat <- levels(r)[[1]]
rat$site <- c('Deramakot','Kuamut','Tangkulap_Pinangah')
levels(r) <- rat
r

# extract values for some cells
i <- extract(r, c(1,2, 25,100))
i

# get the attribute values for these cells
factorValues(r, i)

is.factor(r)

# write to file:
rr <- writeRaster(r, 'site_raster_clipped.tif', overwrite=TRUE)
rr
