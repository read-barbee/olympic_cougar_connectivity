#### Merge Covariate Stack Tiles from Google Earth Engine ####

# Author: Read Barbee

# Date:2023-05-26 

# Purpose: Merge covariate stack tiles from google earth engine into single multiband raster for use in models-- Not currently working


###############################################################################
#### Libraries ####
library(terra)
library(future)
#library(gdalUtilities)

### Import all raster tiles from Google Earth Engine ###

#get list of file paths for each tile
rastlist <- list.files(path = "data/Habitat_Covariates/puma_cov_stack_v1", pattern='.tif$', all.files= T, full.names= T)


#import each tile using terra
tiles <- lapply(rastlist, rast)

#merge the tiles together into single image (very slow)
cov_stack <- do.call(merge, tiles)

#write to raster (also very slow)
#writeRaster(cov_stack, filename = "cov_stack1.tif")



# Set up parallel processing
plan(multicore)  # Use multiple cores for parallel processing


# Plot the raster using parallel processing
plot(cov_stack, parallel = TRUE)


#gdalbuildvrt(rastlist)





