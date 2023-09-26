#### Makah Detection History Formatting 2021 ####

# Author: Read Barbee

# Date:2023-09-26 

# Note: this file already has camera activity and cougar detections combined at the station level


################################ Libraries #################################
library(tidyverse)
library(janitor)


#########################################################################
##
## 1. Import and format camera activity file
##
##########################################################################


cam_act <-read_csv("data/Camera_Data/2022/Makah_2022/Makah_CougarDetections_2022.csv") %>% 
  select(-c(Study)) %>% 
  rename(station_id = Station,
         lon = Longitude,
         lat = Latitude,
         camera_id = CameraID)



#########################################################################
##
## 2. Extract raster cell number for each camera station
##
##########################################################################
# Load the terra package
library(terra)
library(sf)

# Load your raster data
raster <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

temp <- raster[[1]]

# Load your sf object with points for station locations
sf_points <- cam_act %>% st_as_sf(coords=c("lon", "lat"), crs = 4326, remove=FALSE) %>% 
  st_transform(crs = 5070)

# Extract cell numbers for each station location
xy <- sf_points %>% st_coordinates()


#append raster cell id to each station location
cam_act <- cam_act %>% mutate(cell_id = cellFromXY(temp, xy), .after = lon)

get_dupes(cam_act, cell_id)
get_dupes(cam_act, station_id)


#########################################################################
##
## 3. Last bit of formatting
##
##########################################################################

det_hist_binary <- cam_act %>% 
  pivot_longer(cols = c(-c(station_id:camera_id)), names_to = "date", values_to = "value") %>% 
  mutate(date = mdy(date)) %>% 
  pivot_wider(names_from = date, values_from = value) %>% 
  mutate(grid_id = "MAKH",
         year = "2022", .before=station_id)



#write_csv(det_hist_binary, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Camera_Data/2022/Makah_2022/makh_2022_det_hist.csv")




