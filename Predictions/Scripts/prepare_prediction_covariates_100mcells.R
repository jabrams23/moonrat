
library(raster)
library(rgeos)
library(RQGIS)
library(rgdal)
library(vegan)
library(igraph)

# tell RQGIS where to find QGIS
#dir_QGIS <- "C:/Program Files/QGIS 2.18"
#dir_QGIS <- "/Applications/QGIS.app"

#set_env(root = dir_QGIS)

######################################################################
# load land cover map

setwd("~/Dropbox (ScreenForBio)/Projects/moonrat")

#Sabah
wd_landcov_sabah <- "~/Dropbox (ScreenForBio)/Projects/moonrat"
#wd_landcov_sabah <- "C:/Dropbox (ScreenForBio)/Projects/moonrat"
filename_sabah <- "Sabah_25_08_2017.tif"
filepath_landcov_sabah <- file.path(wd_landcov_sabah, filename_sabah)
landcov <- raster(filepath_landcov_sabah)

#TFR <- shapefile("C:/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/sg pinangahfr east with Tangkulap.shp")
#DFR <- shapefile("C:/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/deramakot_fr.shp")

TFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/sg pinangahfr east with Tangkulap.shp")
DFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/deramakot_fr.shp")
KFR <- shapefile("~/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/kuamut fr_west.shp")

TFR <- spTransform(TFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
DFR <- spTransform(DFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
KFR <- spTransform(KFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

combined <- TFR+DFR+KFR
combined<- spTransform(combined, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

landcov_cropped <- crop(landcov, extent(combined))

# set new pixel resolution and calculate aggregation factor
new.resolution <- 100
agg.factor <- new.resolution / unique(res(landcov_cropped))

# forest score: if more than 50% of 5m cells is NA, make whole aggregated cell NA
NA_perc_threshold <- 0.5  

# make forest score data frame. CHECK THIS!!! value = "is"; new_value = "becomes"
forest_score <- data.frame(value = c(11, 12, 22, 23, 31, 51),
                           name = c("Dense Forest", "Forest", "Plantation", "Degraded Areas", "Bare Ground", "Water"),
                           new_value = c(3, 2, 0, 1, 0, NA))

# reclassify land cover to forest score weights
fs_weights <- reclassify(x = landcov_cropped,
                         rcl = as.matrix(forest_score[,c(1,3)]))

# aggregate forest score weights to 100m pixel size (mean per pixel)
fs_agg <- aggregate(fs_weights, fact = agg.factor, fun = mean, na.rm = TRUE)
NA_percent_raster <- aggregate(fs_weights, fact = agg.factor, fun = function(x, ...){mean(is.na(x))})

# apply some threshold so only aggregated cells with a minimum of e.g. 50% actual values are assigned forest scores
values(fs_agg) [values(NA_percent_raster) >= NA_perc_threshold] <- NA   # set values above threshold of percent NA values in aggregated cell = NA

####################################################################
# function for calculating distances from certain land cover classes and aggregating to lower resolution with tempate
# requires RQGIS and raster packages
# landcover is in UTM so there is no issues with imprecision
# maybe think about another solution later. It is good enough for testing though
####################################################################

distanceFromLandcoverRaster <- function(infile,                   # land cover raster
                                        outfile_proximity,        # distance to land cover classes raster
                                        values,                   # pixel values from which to compute distance
                                        raster_template          # raster for which to extract values - an aggregated raster in lower resolution
                                        #  pixelsize_landcov = 5    # the pixel size of the land cover raster
){
  
  alg_proximity <- "gdalogr:proximity"   # name of the WGISalgorithm used
  # get_usage(alg = alg_proximity)
  
  # define function parameters
  params <- get_args_man(alg = alg_proximity)
  
  values_tmp <- paste('"\"', paste(values, collapse  = ","), '\""', sep = "")
  
  params$INPUT <- infile
  params$VALUES <- values_tmp   # value from which to calculate distances (1 = stream)
  params$UNITS <- 0    # GEO = meters
  params$MAX_DIST <- 1000
  params$RTYPE <- 5    # Float32
  params$OUTPUT <- outfile_proximity    # comment this line if proximity raster is not needed on disk. then only temporary file will be created
  
  # run proximity tool
  prox_m <- run_qgis(alg = alg_proximity,
                     params = params,
                     load_output = TRUE)
  
  # get coordinates of cell centres of template raster
  output_raster_xy <- xyFromCell(raster_template, cell = which(!is.na(values(raster_template))), spatial = TRUE)
  
  # extract values of aggregated pixels from the stream shapefile
  distance_landcov <- extract(x = prox_m, y = output_raster_xy)
  
  # write these values into a new raster
  values(raster_template)[which(!is.na(values(raster_template)))] <- distance_landcov
  
  #plot(r_distance_stream)
  return(list(proximity_raster = prox_m,
              aggregated_distance = raster_template))
}



#####################
# land cover example
outfile_tmp <- "distance_landcov_test1.tif"
landcov_value_plantation <- 22  # this is plantation    # can also be several values such as c(11,12)

# distance to plantation
distance_plantation <- distanceFromLandcoverRaster (infile = landcov_cropped,        # land cover raster  (in utm)
                                                    outfile_proximity = outfile_tmp,        # distance to land cover classes raster
                                                    values = landcov_value_plantation,                      # pixel values from which to compute distance
                                                    raster_template = fs_agg                # raster for which to extract values - an aggregated raster in lower resolution
)



###############
# this gives us three covariate layers which can be used for predictions. Water has three layers, so you can choose one and play around with it
#plot(fs_agg)
#plot(distance_plantation[[2]])

# write rasters
#writeRaster(fs_agg, filename = "FS_agg_100m_cell.tif", format = "GTiff", datatype = "FLT4S")
#writeRaster(distance_plantation[[2]], filename = "distance_landcover_plantation_100m_cell.tif", format = "GTiff", datatype = "FLT4S")

#het <- aggregate(landcov_cropped, fact = agg.factor, fun= function(x,...){
#  if(mean(!is.na(x) >= NA_perc_threshold)){
#    tmp <- table(x)
#    out <- diversity(tmp) / log(length(tmp))
#    return(out)
#  } else NA
#})

r3 <- crop(fs_agg, extent(combined))
fs_agg_clipped <- mask(r3, combined)
r4 <- crop(distance_plantation[[2]], extent(combined))
distance_plantation_clipped <- mask(r4, combined)
#r5 <- crop(het, extent(combined))
#het_clipped <- mask(r5, combined)


TRI <- raster('TRI_not_smoothed_proper_scale_150m.tif')

TRI <- projectRaster(TRI, fs, res=100, crs="+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

#new.resolution <- 10
#agg.factor <- new.resolution / unique(res(TRI))

# aggregate forest score weights to 100m pixel size (mean per pixel)
#TRI_agg <- TRI # aggregate(TRI, fact = agg.factor, fun = mean, na.rm = TRUE)
#combined <- spTransform(combined, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
r5 <- crop(TRI, extent(combined))
TRI_clipped <- mask(r5, combined)

# compare extent etc of covariate rasters
stopifnot(compareRaster(fs_agg_clipped, distance_plantation_clipped, TRI_clipped))


plot(fs_agg_clipped)
plot(distance_plantation_clipped)
plot(TRI_clipped)

writeRaster(fs_agg_clipped, filename = "FS_agg_100m_cell_clipped.tif", format = "GTiff", datatype = "FLT4S")
writeRaster(distance_plantation_clipped, filename = "distance_landcover_plantation_100m_cell_clipped.tif", format = "GTiff", datatype = "FLT4S")
writeRaster(TRI_clipped, filename = "TRI_not_smoothed_proper_scale_150m_100mcell_clipped.tif", format = "GTiff", datatype = "FLT4S")
