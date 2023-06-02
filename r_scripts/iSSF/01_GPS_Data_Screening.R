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


############################### Import Location Data #################################
# Mountain lion location data (November 2022)
locs_raw <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_5-16-2023.csv",  col_types = list(fix_type = col_character()))

# Mountain lion deployments  (September 2022)
deployments <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/collar_deployments_master_5-11-2023.csv")

# Dispersal inventory from Teams (3-13-2023)
dispersals <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/dispersals_master_5-11-2023.csv")

#raster file for study area mask
mask <- rast("data/Habitat_Covariates/homerange_habitat_layers_JR_9-14-22/assets/puma_forest_op.tif")

extent <- as.vector(ext(mask))


############################### Data Diagnostics #################################

#calculate global fix rate
global_fix_rate <- function(data){
  
  success <- data %>% 
  filter(is.na(latitude)==FALSE) %>% 
  count() %>% pull()
  
  total = nrow(data)
  
  fix_rate = success/total
  
  return(fix_rate)
  
}

fix_success_rate <- global_fix_rate(locs_raw)

fix_success_rate

#Make sure no essential fields are missing data
summary(locs_raw)

#make sure no locations with coordinates are missing dop scores. Get the names of the individuals that are if any
locs_raw %>% 
  filter(!is.na(latitude) & is.na(dop)) %>% 
  distinct(animal_id)


######## Remove first 24 hours for each animal to avoid capture effects (removes 1,781 points). Remove 2D fixes with dop > 5 (removes 27,397 points)


locs_screened <- locs_raw %>% 
  group_by(animal_id) %>% 
  filter(date_time_local >= (min(date_time_local) + hours(24))) %>% 
  ungroup() %>% 
  filter(!is.na(latitude)) %>% 
  filter(dop<=5 | fix_type == "3D") %>%
  mutate(unique_id= 1:nrow(.), .before = deployment_id) 

#%>% filter(!(unique_id %in% c(44803, 214496, 182176, 132748, 213380, 78225)))

#compare to minimum date times from original data frame to make sure it worked
locs_raw %>%
  group_by(animal_id) %>%
  summarize(min=min(date_time_local))

locs_screened %>%
  group_by(animal_id) %>%
  summarize(min=min(date_time_local))

#Remove locations outside of study area (removes 299 locations)
locs_screened2 <- locs_screened %>% sf::st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=FALSE) %>% 
  sf::st_crop(extent) 

#convert geo-sreened sf object back to dataframe for comparison with df before geo-screenieng
screened_df <- locs_screened2 %>% 
  as.data.frame() %>% 
  select(-geometry)

#get dataframe of differences before and after geo-screening
geo_removed <- setdiff(locs_screened, screened_df)

#geo screening removes 300 points, most of which were in Berlin for Charlotte_28358 in December 2017. Ike also had a number of points in the same location in December 2018. Seems to be an uncorrectable sytematic gps error. Approximately 7 points are due to standard GPS error.
geo_removed %>% sf::st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=FALSE) %>% 
  mapview::mapview()



#write_csv(screened_df, "data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_6-1-2023.csv")

