
#### Bjornerass Screening Test ####

# Author: Read Barbee

# Date:2023-05-16 

# Purpose: Practice running Bjornerass 2010 error screening algorithm. Not currently working



################################ Set screening parameters #################################
delta =100000
mu = 10000
alpha = 1500
theta = (-0.97)

################################ Load libraries and import data #################################
library(tidyverse)
library(sf)


data_raw <- read_csv("data/Location Data/Source Files/locations_master/gps_locs_master_5-16-2023.csv") %>% 
  filter(!is.na(latitude))

#load bjornerass script as function
source("r_scripts/iSSF/Bjornerass_2010_GPS_screening_function.R")

################################ Prepare data for screening #################################

#set initial coordinate system (wgs84)
data_sf <- st_as_sf(data_raw, coords = c("longitude", "latitude"), crs = 4326, remove=FALSE)

#transform to projected coordinate system (Albers)
data_sf <- data_sf %>% st_transform(5070)

#prepare data for the screening function
dat_form <- data_raw %>% 
  mutate(long_5070 = st_coordinates(data_sf)[,1],
         lat_5070 = st_coordinates(data_sf)[,2],
         date_time_utc = mdy_hms(date_time_utc, truncated=3)) %>% 
  rename(lat_wgs84 = latitude,
         long_wgs84 = longitude) %>% 
select(deployment_id:date_time_local, lat_wgs84, long_wgs84, lat_5070, long_5070, everything())


################################ Screening #################################

#run the screening function
bjornerass_screening(dat_form$animal_id, 
                     dat_form$long_5070, 
                     dat_form$lat_5070,
                     dat_form$date_time_utc, delta, mu, alpha, theta)
