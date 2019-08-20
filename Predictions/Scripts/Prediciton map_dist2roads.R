####
# PREDICTION OF OCCURRENCE MAPS
####

wd <- "/Users/jesse/Desktop/robert"

setwd(wd)

library(raster)
library(rgdal)
library(rgeos)

# SITE COVARIATES

covariates <- read.csv(file="/Users/jesse/Desktop/robert/Covariates_MASTERLIST_oDFR+oTFR_2008-2010.csv", row.names=1)

## NEED TO STANDARDIZE ALL CONTINUOUS COVARIATES HERE ##

covariates$site <- relevel(as.factor(covariates$site),ref="Deramakot")


# Use scale-function to scale the columns. Their values are transformed to a scale of the variable.
covariates1 <- as.data.frame(scale(covariates[ , c(4:36)]))

# Put "scale." in front of the column names, allowing later for nice plotting
names(covariates1) <- paste0("scale.", names(covariates1))

# Combine these new variables in covariate table
covariates <- cbind (covariates, covariates1)


# OBSERVATION COVARIATES
effort <- as.matrix(read.csv("/Users/jesse/Desktop/robert/EFFORT_Moonrat_oDFR+oTFR_2008-2010.csv", row.names=1))
road <- as.matrix(read.csv("/Users/jesse/Desktop/robert/ROADSETUP_Moonrat_oDFR+oTFR_2008-2010.csv", row.names=1))
site <- as.matrix(read.csv("/Users/jesse/Desktop/robert/SITE_Moonrat_oDFR+oTFR_2008-2010.csv", row.names=1))


# COVARIATE DATA MUST HAVE SAME TRAPPING SITES AS THE DETECTION HISTORY! NOT MORE!

# DETECTION HISTORY

detecthist <- as.matrix(read.csv(file="/Users/jesse/Desktop/robert/DETECTHIST_Moonrat_oDFR+oTFR_2008-2010.csv", stringsAsFactors = FALSE, row.names = 1))

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
mod <- occu(~ Effort + Road 
            ~ scale.CT_distance_from_roads..m., 
            data = umf)

library(raster)
#?raster
#setwd("C:/Users/guest/Desktop/RB/Moonrat/Data Analysis/R Project - Moonrat/Data/Raster files/")

roadraster <- raster("roads_distance.tif")

## FIND means and SDs of the initial covariates (not the scaled ones). 
# Either do manually from initial covs or transform scaled ones back
#Since we standardized the covariates during the model fitting process, we need to transform the countrywide
# data using the same values. Note, we don’t want to use the mean and SD of the rasters themselves, we
# want to use the mean and SD of the original covariates used to fit the models
covariates11 <- scale(covariates[ , c(4:36)])
attr(covariates11, "scaled:center") # for mean
attr(covariates11, "scaled:scale") # for SD

road.s <- (roadraster-1173.07086)/1045.9296 #subtract mean from raster and divide by SD

plot(road.s, col=terrain.colors(100))

(beta <- coef(mod, type="state"))

logit.psi <- beta[1] + beta[2]*road.s
psi <- exp(logit.psi) / (1 + exp(logit.psi))
plot(psi, col=terrain.colors(100))

TFR <- shapefile("/Users/jesse/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/sg pinangahfr east with Tangkulap.shp")
DFR <- shapefile("/Users/jesse/Dropbox (ScreenForBio)/Projects/community_occupancy/Com_occ_model/new_sabah/reserve shapefiles/deramakot_fr.shp")

TFR <- spTransform(TFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
DFR <- spTransform(DFR, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

roads_v <- shapefile("Roads_Deramakot-Tangkulap-Kuamut_latest_03072017.shp")
roads_v <- spTransform(roads_v, CRS("+proj=utm +zone=50 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

png(filename = "roads_predictions.png", width=6*600, height=6*485, 
    res=600, bg="white")
plot(psi, col=terrain.colors(100))
plot(TFR, add = T)
plot(DFR, add = T)
plot(roads_v, col="blue", add = T)
dev.off()

print(spplot(psi, col.regions=terrain.colors(100)))

spplot(psi, col.regions=terrain.colors(100))

ef <- stack(road.s)
names(ef) <- c("scale.CT_distance_from_roads..m.")

E.psi <- predict(mod, type="state", newdata=ef)
plot(E.psi, axes=FALSE, col=terrain.colors(100))


