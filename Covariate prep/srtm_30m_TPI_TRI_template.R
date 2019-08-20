TerrainMeasuresAtVariousScales <- function(dem,
                                           matrix_sizes = c(seq(3,13,2), 23, 33, 43, 53, 67),   # size of extraction matrix in cells (3x3 = 90*90m, 67x67 = 2km*2km at 30m resolution)
                                           measure = c("TRI", "TPI", "roughness"),
                                           dir_outout){

  # guess cell size
  if(round(res(dem)[1], 4) == 3e-04){
    dem_cellsize <- 30
  } else {
    if(round(res(dem)[1], 4) == 8e-04) {
      dem_cellsize <- 90
    } else {
      stop("cannot guess dem cell size. I want a 30m or 90m DEM!")
    }
  }

  measure <- match.arg(measure, choices = c("TRI", "TPI", "roughness")) #, choices = c("TRI", "TPI", "roughness"))

  for(matrix_size in matrix_sizes){  # for neighbourhoods of 3:67 cells (90x90m to 2010x2010m)

    message(paste(matrix_size, "x", matrix_size))

    #weight_matrix <- matrix(c(rep(1, times = floor((matrix_size^2)/2)), 0, rep(1, times = floor((matrix_size^2)/2))), nrow = matrix_size, ncol = matrix_size)
    weight_matrix <- matrix(1, nrow = matrix_size, ncol = matrix_size)

    which_matrix_0 <- ceiling((matrix_size^2) / 2)
    sum_weights <- (sum(weight_matrix) - 1)

    if(measure == "TRI") outRaster <- focal(x = dem, w=weight_matrix, fun=function(x){sum(abs(x[-which_matrix_0]-x[which_matrix_0]))/sum_weights})
    if(measure == "TPI") outRaster <- focal(x = dem, w=weight_matrix, fun=function(x){x[which_matrix_0] - mean(x[-which_matrix_0])})
    if(measure == "roughness") outRaster  <- focal(x = dem, w=weight_matrix, fun=function(x){max(x) - min(x)})

    try(writeRaster(outRaster, filename = file.path(dir_outout, paste(measure, "_not_smoothed_proper_scale_",  matrix_size * dem_cellsize, "m.tif", sep = ""))))

    rm(outRaster, which_matrix_0)

  }
}


# # example calls:
#
r.dem30 <- raster("~/Dropbox (ScreenForBio)/Projects/moonrat/n05_e117_1arc_v3.tif")    # change to the correct file (see email 2017-02-02)
dir_tpi_tri <- "~/Dropbox (ScreenForBio)/Projects/moonrat"   # define the directory to save TRI and TPI rasters to!

TerrainMeasuresAtVariousScales(dem = r.dem30, measure = "TRI", dir_outout = dir_tpi_tri, matrix_sizes = c(5,9,13))
