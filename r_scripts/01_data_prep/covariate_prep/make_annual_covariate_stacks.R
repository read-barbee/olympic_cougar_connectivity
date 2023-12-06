#### Create annual covariate stacks ####

# Author: Read Barbee

# Date:2023-11-07

# Purpose: Create stacks of each covariate for each year to use in SSF and occupancy analysis


###############################################################################
#### Libraries ####
library(tidyverse)
library(terra)


#define years for annual stacks to output
years <- 2010:2023

# Define an output directory
output_dir <- "data/Habitat_Covariates/annual_covariates/annual_stacks"

output_dir2 <- "data/Habitat_Covariates/annual_covariates/annual_stacks_by_cov"


################################ Helper functions #################################

#function to rename image bands by year
band_rename <- function(image, base_name, start_year, end_year){
  old_names <- names(image)
  years <- as.character(seq(from = start_year, to = end_year, by = 1))
  
  new_names <- vector()
  for(i in 1:length(old_names)){
    new_names[i] <- paste0(base_name, "_", years[i])
  }
  
  names(image) <- new_names
  
  return(image)
}

#assign same projection to all images (takes care of PROJCRS vs ENGCRS problems from gdal_merge)
crs_fun <- function(x){crs(x) <- "EPSG:5070"
return(x)}

#collect images from each cov stack for a single year. If the year isn't in the stack, use the nearest one in time
make_annual_stack <- function(image_list, year){
  
  annual_stack <- list()
  
  for(i in 1:length(cov_list)){
    
    img_names <- names(image_list[[i]])
    
    img_start_year <- as.numeric(str_sub(first(img_names), start = -4))
    
    img_end_year <- as.numeric(str_sub(last(img_names), start = -4))
    
    if(any(str_detect(img_names, coll(as.character(year))))==T){
      
      annual_stack[[i]] <- terra::subset(image_list[[i]], str_detect(img_names, coll(as.character(year))))
    }else if(year < img_start_year){
      
      annual_stack[[i]]<- terra::subset(image_list[[i]], first(img_names))
    }else if(year > img_end_year){
      
      annual_stack[[i]] <- terra::subset(image_list[[i]], last(img_names))
    }
    
    print(paste0(i, "/", length(cov_list)))
  }
  
  annual_stack_c <- rast(annual_stack)
  
  return(annual_stack_c)
}

######################################################################
##
## 1. Import raster stacks for each covariate
##
##########################################################################

### Static Covariates ###

elevation <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/static/elevation.tif")

slope <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/static/slope.tif")

aspect <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/static/aspect.tif")

mtpi <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/static/mtpi.tif")

tri <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/static/tri.tif")

dist_water <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/static/dist_water.tif")

#make stack of static covariates
static_stack <- c(elevation, slope, aspect, mtpi, tri, dist_water)

#replace -999 no data values with NA
static_stack <- mask(static_stack, static_stack, maskvalues = -999, updatevalue=NA)

#calculate northness and eastness layers from aspect 
static_stack$aspect_northness <- cos((pi*static_stack$aspect)/180)
static_stack$aspect_eastness <- sin((pi*static_stack$aspect)/180)

#reorder layers
static_stack <- static_stack[[c('elevation', 'slope', 'aspect', 'aspect_northness', 'aspect_eastness', 'mtpi', 'tri', 'distance_water')]]


### Annual Covariates ##

ndvi_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/modis_ndvi_annual_2010_2023.tif") %>% 
  band_rename("ndvi_annual", 2010, 2023) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

evi_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/modis_evi_annual_2010_2023.tif") %>% 
  band_rename("evi_annual", 2010, 2023) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

gpp_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/gpp_annual_2010_2022.tif") %>% 
  band_rename("gpp_annual", 2010, 2022) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

npp_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/npp_annual_2010_2022.tif") %>% 
  band_rename("npp_annual", 2010, 2022) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

perc_tree_cov_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/perc_tree_cov_annual_2010_2020.tif") %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

perc_nontree_veg_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/perc_nontree_veg_annual_2010_2020.tif") %>%
  mask(., ., maskvalues = -999, updatevalue=NA)

perc_nonveg_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/perc_nonveg_annual_2010_2020.tif") %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

precip_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/precip_annual_2010_2023.tif") %>% 
  band_rename("precip_annual", 2010, 2023) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

land_cover_usfs_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/land_cover_usfs_annual_2010_2022.tif") %>% 
  band_rename("land_cover_usfs_annual", 2010, 2022) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

land_use_usfs_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/land_use_usfs_annual_2010_2022.tif") %>% 
  band_rename("land_use_usfs_annual", 2010, 2022) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

land_use_change_usfs_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/land_use_change_usfs_annual_2010_2022.tif") %>% 
  band_rename("land_use_change_usfs_annual", 2010, 2022) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

popdens_hii_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/popdens_hii_annual_2010_2020.tif") %>% 
  band_rename("popdens_hii_annual", 2010, 2020) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)  * 0.01

landuse_hii_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/landuse_hii_annual_2010_2020.tif") %>% 
  band_rename("landuse_hii_annual", 2010, 2020) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)  * 0.01

infra_hii_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/infra_hii_annual_2010_2020.tif") %>% 
  band_rename("infra_hii_annual", 2010, 2020) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)  * 0.01

roads_hii_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/roads_hii_annual_2010_2022.tif") %>% 
  subset(12, negate = TRUE) %>% #remove 2022 layer that the other hii indices don't have
  band_rename("roads_hii_annual", 2010, 2020) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)  * 0.01

hii_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/hii_annual_2010_2020.tif") %>% 
  band_rename("hii_annual", 2010, 2020) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)  * 0.01

dist_major_roads_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/dist_major_roads_annual_2014_2023.tif") %>% 
  band_rename("dist_major_roads_annual", 2014, 2023) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

dist_minor_roads_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/dist_minor_roads_annual_2014_2023.tif") %>% 
  band_rename("dist_minor_roads_annual", 2014, 2023) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

dist_very_small_roads_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/dist_very_small_roads_annual_2014_2023.tif") %>% 
  band_rename("dist_very_small_roads_annual", 2014, 2023) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

dist_non_motorized_roads_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/dist_non_motorized_roads_annual_2014_2023.tif") %>% 
  band_rename("dist_non_motorized_roads_annual", 2014, 2023) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

dist_all_roads_annual <- rast("data/Habitat_Covariates/annual_covariates/raw_collections/annual/dist_all_roads_annual_2014_2023.tif") %>% 
  band_rename("dist_all_roads_annual", 2014, 2023) %>% 
  mask(., ., maskvalues = -999, updatevalue=NA)

######################################################################
##
## 2. Make annual stacks
##
##########################################################################

#make list of all individual covariate stacks to iterate over
cov_list <- list(ndvi_annual,
                 evi_annual,
                 gpp_annual,
                 npp_annual,
                 perc_tree_cov_annual,
                 perc_nontree_veg_annual,
                 perc_nonveg_annual,
                 precip_annual,
                 land_cover_usfs_annual,
                 land_use_usfs_annual,
                 land_use_change_usfs_annual,
                 popdens_hii_annual,
                 landuse_hii_annual,
                 infra_hii_annual,
                 roads_hii_annual,
                 hii_annual,
                 dist_major_roads_annual,
                 dist_minor_roads_annual,
                 dist_very_small_roads_annual,
                 dist_non_motorized_roads_annual,
                 dist_all_roads_annual
                 )

names(cov_list) <- c("ndvi_annual",
                     "evi_annual",
                     "gpp_annual",
                     "npp_annual",
                     "perc_tree_cov_annual",
                     "perc_nontree_veg_annual",
                     "perc_nonveg_annual",
                     "precip_annual",
                     "land_cover_usfs_annual",
                     "land_use_usfs_annual",
                     "land_use_change_usfs_annual",
                     "popdens_hii_annual",
                     "landuse_hii_annual",
                     "infra_hii_annual",
                     "roads_hii_annual",
                     "hii_annual",
                     "dist_major_roads_annual",
                     "dist_minor_roads_annual",
                     "dist_very_small_roads_annual",
                     "dist_non_motorized_roads_annual",
                     "dist_all_roads_annual")


#make sure all images are in the same CRS
cov_list<- map(cov_list, crs_fun)


#Make an annual stack for each year in the "years" vector
annual_stacks <- list()

for(i in years){
  annual_stacks[[paste0("stack_", i)]]<- make_annual_stack(cov_list, i)
  
  print(i)
}

######################################################################
##
## 3. Export stacks as tif files
##
##########################################################################

# Loop through the list and export each spatraster
for (i in 1:length(annual_stacks)) {
  raster <- annual_stacks[[i]]
  file_name <- paste0("cov_",names(annual_stacks)[i], ".tif")  # Define the file name
  output_path <- file.path(output_dir, file_name)  # Full output path

  # Export the spatraster as a TIFF file
  writeRaster(raster, output_path, filetype = "GTiff", overwrite = TRUE)

  cat("Exported:", output_path, "\n")  # Print the exported file path
}



# Output formatted annual stacks by covariate
for (i in 1:length(cov_list)) {
  raster <-cov_list[[i]]
  file_name <- paste0(names(cov_list)[i], ".tif")  # Define the file name
  output_path <- file.path(output_dir2, file_name)  # Full output path
  
  # Export the spatraster as a TIFF file
  writeRaster(raster, output_path, filetype = "GTiff", overwrite = TRUE)
  
  cat("Exported:", output_path, "\n")  # Print the exported file path
}

#output static stack

#writeRaster(static_stack, file.path(output_dir2, "static_stack.tif"), filetype = "GTiff", overwrite = TRUE)


