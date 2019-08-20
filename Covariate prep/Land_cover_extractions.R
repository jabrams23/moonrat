#################################################################################
#### land cover class extraction  -------------------############################
#################################################################################


# To do: - check if projection argument is necessary

extractLandCover <- function(landCoverRaster,                 # land cover raster
                             pointsLayer,                     # spatialPointsDataFrame
                             extractionRadii,                 # in meters (one or more)
                             pointsLayerIDcolumnName,         # station ID in pointsLayer
                             landCoverRasterClasses,          # optional vector with the land cover value names in landCoverRaster
                             dir_csv,                         # directory to save csv table to
                             writeShapefile,                  # create shapefiles of extraction radii?
                             dir_extraction_radii_shapefiles, # directory to save shapefiles
                             forestScoreWeights)              # data frame contains at least 2 columns: value and weight. value is land cover value, weight is associated weight for computation of forest score (0...3, or NA)
{
  
  if(projection(landCoverRaster) != projection(pointsLayer)) warning("landCoverRaster and pointsLayer need to have same projection. Check crs() or projection() of both",
                                                                     immediate. = TRUE, call. = FALSE)
  stopifnot(dir.exists(dir_csv))
  if(isTRUE(writeShapefile)){
    stopifnot(dir.exists(dir_extraction_radii_shapefiles))
    if(length(list.files(dir_extraction_radii_shapefiles )) != 0) stop("there's stuff in dir_extraction_radii_shapefiles. Please provide emtpy directory", call. = FALSE)
  }
  
  
  extract.mask <- extract.mask.df <- layer.names <- list()
  
  for(i in 1:length(extractionRadii)){
    extract.mask[[i]] <- gBuffer(pointsLayer, byid = TRUE,
                                 width = extractionRadii[i], quadsegs=10)
    
    extract.mask.df[[i]] <- SpatialPolygonsDataFrame(extract.mask[[i]], data = pointsLayer@data)       # create spatial polygon
  }
  names(extract.mask.df) <- extractionRadii
  
  
  
  if(any(is.na(extract(landCoverRaster, pointsLayer)))) stop ("at least one point is outside landCoverRaster")    # check that all ct points have background
  
  landcov.extract <- list()
  
  for(i in 1:length(extract.mask)){
    message(extractionRadii[i])
    landcov.extract[[i]] <- extract(landCoverRaster, extract.mask[[i]])
  }
  
  table_landcover <- table(values(landCoverRaster))
  n_landcov_classes <- length(table_landcover)
  values_landcov_classes <- names(table_landcover)
  
  
  # prepare table
  
  d <- data.frame(matrix(data = 0,
                         ncol = n_landcov_classes ,
                         nrow = length(pointsLayer) * length(extractionRadii)))
  
  colnames(d) <- values_landcov_classes
  
  d.info <- expand.grid(extractionRadii, pointsLayer@data[,pointsLayerIDcolumnName])
  colnames(d.info)[1] <- "radius"
  colnames(d.info)[2] <- pointsLayerIDcolumnName
  
  
  # percentages of each class
  
  for(i in 1:length(landcov.extract)){    # loop through distances
    for(j in 1:length(pointsLayer)){      # loop through sites
      
      tmp <- table(unlist(landcov.extract[[i]][j]))                          # get pixel values of station / distance as vector
      
      index.val <- intersect(which(d.info$radius ==  extractionRadii[i]),    # get row index of data frame d
                             which(d.info[,pointsLayerIDcolumnName] ==  pointsLayer@data[j, pointsLayerIDcolumnName]))
      
      d[index.val, match(as.integer(names(tmp)), names(d))] <- tmp
    }
  }
  
  n_cells <-  rowSums(d)
  
  d_perc <- d / rowSums(d)
  
  # calculate Forest Score
  
  if(nrow(forestScoreWeights) != n_landcov_classes) stop("mismatch between number of values in forestScoreWeights and number of classes in landcoverRaster")
  if(any(forestScoreWeights$weight < 0, na.rm = TRUE) | any(forestScoreWeights$weight > 3, na.rm = TRUE)) stop("forestScoreWeights$weight must be between 0 and 3, or NA")
  if(!all(colnames(d_perc) == as.character(forestScoreWeights$value))) stop("please sort forestScoreWeights by forestScoreWeights$value")
  
  # if any weight = NA (= omit class), correct percentage data frame
  d_perc2 <- d_perc
  d_perc2[,which(!is.na(forestScoreWeights$weight))] <- d[,which(!is.na(forestScoreWeights$weight))] / rowSums(d[,which(!is.na(forestScoreWeights$weight))])
  d_perc2[,which(is.na(forestScoreWeights$weight))]  <- 0
  
  
  fs <- rep(NA, times = nrow(d_perc))
  for(m in 1:nrow(d_perc)){
    #fs[m] <- sum(d_perc[m,] * forestScoreWeights$weight, na.rm = TRUE)  # this does not respect land cover classes with weight = NA (which are to be ignored)
    fs[m] <- sum(d_perc2[m,] * forestScoreWeights$weight, na.rm = TRUE) # this one does
  }
  
  
  if(hasArg(landCoverRasterClasses)) colnames(d) <- colnames(d_perc) <- landCoverRasterClasses
  
  
  # heterogeneity: computed as Pielou's evenness
  # high evenness = even distibution of classes/no dominant class = heterogeneous habitat
  
  H <- diversity(d)           # Shannon index
  J <- H/log(specnumber(d))   # Pielou's evenness
  
  
  
  # richness (of habitat types)
  spec_num <- specnumber(d)
  
  d3 <-data.frame(d.info,
                  round(d_perc, digits = 3),
                  n_cells = n_cells,
                  forest_score = fs,
                  heterogeneity = round(J, 2),
                  n_habitat = spec_num
  )
  
  write.csv(d3, file = file.path(dir_csv, paste("Landcover_extraction_", Sys.Date(), ".csv")))
  
  # append data to shapefiles
  
  if(isTRUE(writeShapefile)) {
    
    for(n in 1:length(extract.mask.df)){
      
      if(!all(extract.mask.df[[n]]@data[,pointsLayerIDcolumnName] == d3[d3$radius == names(extract.mask.df)[n],pointsLayerIDcolumnName])) stop("weirdass bug. shouldnt have happened. Something wrong in order of stations in tables")
      extract.mask.df[[n]]@data <- cbind(extract.mask.df[[n]]@data, d3[d3$radius == names(extract.mask.df)[n],-2])
      
      
      writeOGR(extract.mask.df[[n]], dsn = dir_extraction_radii_shapefiles, layer = paste("extraction_mask_circle_", extractionRadii[n], "m", sep = ""),      # save extraction mask for use in gis
               driver = "ESRI Shapefile")
      
    }
  }
  
  return(d3)
}

##################################################################################################
#####  loading sample data and running functions
##################################################################################################

library(raster)
library(rgdal)
library(rgeos)
library(vegan)
library(igraph)

setwd("/Users/jesse/Dropbox (ScreenForBio)/Projects/moonrat")

Sabah_fs_100m <- raster("FS_agg_100m_cell.tif")

Sabah_final <- raster("Sabah_25_08_2017.tif")


TFR <- shapefile("/Users/jesse/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/sg pinangahfr east with Tangkulap.shp")
DFR <- shapefile("/Users/jesse/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/deramakot_fr.shp")

TFR <- spTransform(TFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
DFR <- spTransform(DFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

combined <- TFR+DFR
combined<- spTransform(combined, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

plot(Sabah_fs_100m)
plot(combined,add=T)

r2 <- crop(Sabah_fs_100m, extent(combined))
r3 <- mask(r2, combined)

plot(r3)

bb <- bbox(TFR)
cs <- c(3.28084, 3.28084)*100  # cell size 6km x 6km (for illustration)
# 1 ft = 3.28084 m
cc <- bb[, 1] + (cs/2)  # cell offset
cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
grd

sp_grd <- SpatialGridDataFrame(grd,
                               data=data.frame(id=1:prod(cd)),
                               proj4string=CRS(proj4string(TFR)))

############################################################################################

forestScoreWeightsDF <- data.frame(value = c(11, 12, 22, 23, 31, 51),
                           name = c("Dense Forest", "Forest", "Plantation", "Degraded Areas", "Bare Ground", "Water"),
                           new_value = c(3, 2, 0, 1, 0, 0))

LCC_extraction_Vietnam <- extractLandCover(landCoverRaster = Sabah_final,    # land cover raster
                                  pointsLayer = TFR,
                                  extractionRadii = c(50),    # in meters
                                  pointsLayerIDcolumnName = "centroids",
                                  forestScoreWeights = forestScoreWeightsDF,
                                  landCoverRasterClasses = forestScoreWeightsDF$name,    # optional vector with the land cover value names in landCoverRaster
                                  dir_csv = ".",
                                  writeShapefile = TRUE,
                                  dir_extraction_radii_shapefiles = "/Users/jesse/Dropbox (ScreenForBio)/Projects/moonrat/extraction")


test <- aggregate(Sabah_final, fact = 20, fun= function(x){
  tmp <- table(values(x))
  out <- diversity(tmp) / log(length(tmp))
  return(out)
}
)


agg.factor=10

het<- aggregate(landcov, fact = agg.factor, fun= function(x,...){
  if(mean(!is.na(x) >= NA_perc_threshold)){
    tmp <- table(x)
    out <- diversity(tmp) / log(length(tmp))
    return(out)
  } else NA
})
