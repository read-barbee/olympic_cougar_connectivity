#### OCP iSSF Module_01: GPS Data Screening ####

# Author: Read Barbee

# Date:2023-06-01

# Inputs:
#   •	GPS0: Raw GPS locations
# 
# Outputs:
#   •	GPS1: Screened GPS locations
# 
# Steps:
# •	Calculate global fix rate 
# •	Check for missing data
# •	Remove first 24hr of fixes for each animal to avoid capture effects
# •	Remove 2D fixes with dop score > 5
# •	Remove locations outside extent of study area raster 




################################ Libraries #################################
library(tidyverse)
library(terra)
library(sf)
library(lubridate)
library(amt)
library(janitor)

################################ User-Defined Parameters #################################


#########################################################################
##
## 1. Import and format location data
##
##########################################################################
# Mountain lion location data (April 2023)
locs_raw <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_6-29-2023.csv",  col_types = list(fix_type = col_character()))

# Mountain lion deployments  (April 2023)
deployments <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Deployments/collar_deployments_master_5-11-2023.csv")

# Dispersal inventory from Teams (3-13-2023)
dispersals <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_5-11-2023.csv")

#raster file for study area mask
mask <- rast("data/Habitat_Covariates/homerange_habitat_layers_JR_9-14-22/assets/puma_forest_op.tif")

extent <- as.vector(ext(mask))


#########################################################################
##
## 2. Data Diagnostics
##
##########################################################################

#calculate global fix rate
global_fix_rate <- function(data){
  
  success <- data %>% 
  filter(is.na(lat_wgs84)==FALSE) %>% 
  count() %>% pull()
  
  total = nrow(data)
  
  fix_rate = success/total
  
  return(fix_rate)
  
}

fix_success_rate <- global_fix_rate(locs_raw)

fix_success_rate

#Make sure no essential fields are missing data
summary(locs_raw %>% filter(!is.na(lat_wgs84)))

#make sure no locations with coordinates are missing dop scores. Get the names of the individuals that are if any
locs_raw %>% 
  filter(!is.na(lat_wgs84) & is.na(dop)) %>% 
  distinct(animal_id)

#########################################################################
##
## 3. Remove missing fixes and locations with missing DOP scores
##
##########################################################################

#remove missing locations
locs_raw_no_na <- locs_raw %>% filter(!is.na(lat_wgs84))

#remove the 22 locations that are missing dop scores

locs_raw_filt <- locs_raw_no_na %>% 
  filter(!is.na(lat_wgs84) & !is.na(dop))

#########################################################################
##
## 4. Remove capture effects and low quality locations
##
##########################################################################


######## Remove first 24 hours for each animal to avoid capture effects (removes 1,781 points). Remove 2D fixes with dop > 5 (removes 27,397 points)


cap_eff <- locs_raw_filt %>% 
  group_by(animal_id) %>% 
  filter(date_time_local >= (min(date_time_local) + hours(24))) %>% 
  ungroup() %>% 
  filter(!is.na(lat_wgs84))
  
  
dop_filt <- cap_eff %>%  filter(dop<=5 | fix_type == "3D") %>%
  mutate(unique_id= 1:nrow(.), .before = deployment_id) 

#%>% filter(!(unique_id %in% c(44803, 214496, 182176, 132748, 213380, 78225)))

#compare to minimum date times from original data frame to make sure it worked
locs_raw %>%
  group_by(animal_id) %>%
  summarize(min=min(date_time_local))

cap_eff %>%
  group_by(animal_id) %>%
  summarize(min=min(date_time_local))

#Remove locations outside of study area (removes 299 locations)
geo_filt <- dop_filt %>% sf::st_as_sf(coords=c("lon_wgs84", "lat_wgs84"), crs=4326, remove=FALSE) %>% 
  sf::st_crop(extent) 

#convert geo-sreened sf object back to dataframe for comparison with df before geo-screenieng
geo_screened_df <- geo_filt %>% 
  as.data.frame() %>% 
  select(-geometry)

#get dataframe of differences before and after geo-screening
geo_removed <- setdiff(dop_filt, geo_screened_df)

#geo screening removes 298 points, most of which were in Berlin for Charlotte_28358 in December 2017. Ike also had a number of points in the same location in December 2018. Seems to be an uncorrectable sytematic gps error. Approximately 7 points are due to standard GPS error.
geo_removed %>% sf::st_as_sf(coords=c("lon_wgs84", "lat_wgs84"), crs=4326, remove=FALSE) %>% 
  mapview::mapview()

location_attempts <- nrow(locs_raw_no_na)
missing_locations <- nrow(locs_raw) - nrow(locs_raw_no_na)
unfiltered <- nrow(locs_raw_no_na)
missing_dop <- nrow(locs_raw_no_na) - nrow(locs_raw_filt)
capture_effects <- nrow(locs_raw_filt) - nrow(cap_eff)
dop_filter <- nrow(cap_eff) - nrow(dop_filt)
geo_filter <- nrow(dop_filt) - nrow(geo_screened_df)
remaining <- nrow(geo_screened_df)

summary_table <- tibble(Step = c("Location Attempts",
                                 "Missing Locations",
                                 "Success Rate",
                                 "Total Locations", 
                                 "Missing DOP", 
                                 "Capture effects", 
                                 "DOP filter", 
                                 "Geo filter", 
                                 "Remaining"),
                        Locations_Removed = c(location_attempts,
                                              missing_locations,
                                              fix_success_rate,
                                              unfiltered,
                                              missing_dop,
                                              capture_effects,
                                              dop_filter,
                                              geo_filter,
                                              remaining))

summary_table


#write_csv(geo_screened_df, "data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_7-3-2023.csv")

