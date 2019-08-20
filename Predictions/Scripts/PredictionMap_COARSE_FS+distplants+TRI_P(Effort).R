####
# PREDICTION OF OCCURRENCE MAPS
####
# COARSE P(Effort+Road): 2014-16 COARSE, forest score + dist from plantations

library(raster)
library(rgdal)
library(rgeos)

# SITE COVARIATES

covariates <- read.csv(file="file:///C:/Users/guest/Dropbox (ScreenForBio)/Projects/moonrat/Modeling results/Moonrat/Coarse grid only/Data/Covariates_MASTERLIST_ALL_2014-2016_COARSE.csv", row.names=1, stringsAsFactors = F)

str(covariates)

# SORT covariates based on station ID
covariates <- covariates[order(covariates$station),]

## NEED TO STANDARDIZE ALL CONTINUOUS COVARIATES HERE ##

covariates$Elevation..DEM.in.m. <- as.numeric(covariates$Elevation..DEM.in.m.)
covariates$roughness_not_smoothed_proper_scale_150m <- as.numeric(covariates$roughness_not_smoothed_proper_scale_150m)
covariates$roughness_not_smoothed_proper_scale_270m <- as.numeric(covariates$roughness_not_smoothed_proper_scale_270m)
covariates$roughness_not_smoothed_proper_scale_390m <- as.numeric(covariates$roughness_not_smoothed_proper_scale_390m)


covariates$site <- relevel(as.factor(covariates$site),ref="Deramakot")

covariates$tri150meters <- covariates$TRI_not_smoothed_proper_scale_150m
covariates$tpi150meters <- covariates$TPI_not_smoothed_proper_scale_150m
covariates$roughness150meters <- covariates$roughness_not_smoothed_proper_scale_150m


str(covariates)
# Use scale-function to scale the continuous variables. 
covariates1 <- as.data.frame(scale(covariates[ , c(4:38)]))

# Put "scale." in front of the column names, allowing later for nice plotting
names(covariates1) <- paste0("scale.", names(covariates1))

# Combine these new variables in covariate table
covariates <- cbind (covariates, covariates1)
str(covariates)

# OBSERVATION COVARIATES
effort <- as.matrix(read.csv("file:///C:/Users/guest/Dropbox (ScreenForBio)/Projects/moonrat/Modeling results/Moonrat/Coarse grid only/Data/EFFORT_Moonrat_COARSE_2014-2016.csv", row.names=1))
road <- as.matrix(read.csv("file:///C:/Users/guest/Dropbox (ScreenForBio)/Projects/moonrat/Modeling results/Moonrat/Coarse grid only/Data/ROADSETUP_Moonrat_COARSE_2014-2016.csv", row.names=1))
site <- as.matrix(read.csv("file:///C:/Users/guest/Dropbox (ScreenForBio)/Projects/moonrat/Modeling results/Moonrat/Coarse grid only/Data/SITE_Moonrat_COARSE_2014-2016.csv", row.names=1))
#survey <- as.matrix....

#SORT those too based on station ID (which is row names here). The model will look at the ORDER of observations! Not the station ID itself!
effort <- effort[order(row.names(effort)),]
road <- road[order(row.names(road)),]
site <- site[order(row.names(site)),]

# COVARIATE DATA MUST HAVE SAME TRAPPING SITES AS THE DETECTION HISTORY! NOT MORE!

# DETECTION HISTORY

detecthist <- as.matrix(read.csv(file="file:///C:/Users/guest/Dropbox (ScreenForBio)/Projects/moonrat/Modeling results/Moonrat/Coarse grid only/Data/DETECTHIST_Moonrat_COARSE_2014-2016.csv", stringsAsFactors = FALSE, row.names = 1))

#SORT
detecthist <- detecthist[order(row.names(detecthist)),]


# UNMARKED #
library(unmarked)

# Create unmarked data frame for occupancy analysis
umf <- unmarkedFrameOccu(y = detecthist,                            
                         siteCovs =  covariates,
                         obsCovs = list(Effort = as.matrix(effort), 
                                        Road = as.matrix(road),
                                        Site = as.matrix(site)
                         )
                         #we put siteCovs as dat so that ALL of them are available for the analysis and we can pick and choose subsets from this 
                         #as.data.frame(siteCovs.tmp)                           
                         #obsCovs = list(effort = scale.effort))
                         #obsCovs = list(effort = as.matrix(scale.effort)) #this also seems to work! 
)


#RUN BEST MODEL
mod <- occu(~ Effort
            ~ scale.forest_score_50m_radius+scale.CT_distance_from_plantations+scale.TRI_not_smoothed_proper_scale_150m, 
            data = umf)

library(raster)
#?raster
#setwd("C:/Users/guest/Desktop/RB/Moonrat/Data Analysis/R Project - Moonrat/Data/Raster files/")
setwd(".")
fsraster <- raster("C:/Users/guest/Dropbox (ScreenForBio)/Projects/moonrat/Prediction_rasters/FS_agg_100m_cell_clipped.tif")
plantraster <- raster ("C:/Users/guest/Dropbox (ScreenForBio)/Projects/moonrat/Prediction_rasters/distance_landcover_plantation_100m_cell_clipped2.tif")
TRIraster <- raster("C:/Users/guest/Dropbox (ScreenForBio)/Projects/moonrat/Prediction_rasters/TRI_not_smoothed_proper_scale_150m_100mcell_clipped.tif")

#CLIP the rasters to the extent we need

#TFR <- shapefile("./Data/Rasters and shapefiles/reserve shapefiles/sg pinangahfr east with Tangkulap.shp")
#DFR <- shapefile("./Data/Rasters and shapefiles/reserve shapefiles/deramakot_fr.shp")

#TFR <- spTransform(TFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#DFR <- spTransform(DFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

#combined <- TFR+DFR
#combined<- spTransform(combined, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

#plot(fsraster)
#plot(combined,add=T)

#r2 <- crop(fsraster, extent(combined))
#r3 <- mask(r2, combined)

#plot(r3)

## FIND means and SDs of the initial covariates (not the scaled ones). 
# Either do manually from initial covs or transform scaled ones back
#Since we standardized the covariates during the model fitting process, we need to transform the countrywide
# data using the same values. Note, we don’t want to use the mean and SD of the rasters themselves, we
# want to use the mean and SD of the original covariates used to fit the models
covariates11 <- scale(covariates[ , c(4:38)])
attr(covariates11, "scaled:center") # for mean
attr(covariates11, "scaled:scale") # for SD

fs.s <- (fsraster-2.1402233)/0.4780219 #subtract mean from raster and divide by SD
plant.s <- (plantraster-6110.9541455)/3341.9843782
TRI.s <- (TRIraster-7.3678241)/4.8369076

ef <- stack(fs.s,plant.s, TRI.s)
names(ef) <- c("scale.forest_score_50m_radius", "scale.CT_distance_from_plantations", "scale.TRI_not_smoothed_proper_scale_150m")


plot(ef, col=terrain.colors(100))

#Get coefficients from the occu model
(beta <- coef(mod, type="state"))

#Make new object. Multiply coefficients with the respective raster objects.
logit.psi <- beta[1] + beta[2]*fs.s + beta[3]*plant.s
#I think this is to backtransform the coefficients
psi <- exp(logit.psi) / (1 + exp(logit.psi))
plot(psi, col=terrain.colors(100))

#png(filename = "roads_predictions.png", width=6*600, height=6*485, 
#    res=600, bg="white")
#plot(psi, col=terrain.colors(100))
#plot(TFR, add = T)
#plot(DFR, add = T)
#plot(roads_v, col="blue", add = T)
#dev.off()

print(spplot(psi, col.regions=terrain.colors(100)))

spplot(psi, col.regions=terrain.colors(100))


E.psi <- predict(mod, type="state", newdata=ef)

plot(E.psi, axes=FALSE, col=terrain.colors(100))
##Plot. Set zlim to sensbile min and max psi values across all maps (same across all maps)
plot(E.psi,zlim=c(0,1), axes=FALSE,
     #col= rev(rainbow(100, start=0.0, end =0.2))
     col=colorRampPalette(c("grey","yellow", "orange", 
                            "red","darkred"
     ))(255))

save(E.psi, file="PredictionResults_COARSE_fs+disttoplants+TRI_P(Effort).RData")


