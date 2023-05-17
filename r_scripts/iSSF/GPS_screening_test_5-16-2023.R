
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
select(deployment_id:date_time_local, lat_wgs84, long_wgs84, lat_5070, long_5070, everything()) %>% 
  filter(animal_id!="Belle") #Belle doesn't have enough locations for adehabitat



################################ Screening #################################

#run the screening function
screened <- bjornerass_screening(dat_form$animal_id, 
                     dat_form$long_5070, 
                     dat_form$lat_5070,
                     dat_form$date_time_utc, delta, mu, alpha, theta) 
  

#default parameters remove 949 locations
screened <- screened %>% 
  rename(animal_id = id,
         long_5070 = x,
         lat_5070 = y,
         date_time_utc = date)


#both of these semi-join calls seem to do the same thing, but I'll use the one that's more specific
#semi_join(dat_form, screened)

locs_error_screened <- semi_join(dat_form, screened, by=join_by(animal_id, long_5070, lat_5070, date_time_utc))



################################ Debugging #################################

# animal_ids <- dat_form %>% distinct(animal_id) %>% pull()
# 
# 
# # Apply function to each element of the list with error handling
# for (i in seq_along(animal_ids)) {
#   tryCatch({
#     dat2 <- dat_form %>% filter(animal_id == animal_ids[i])
#     # Call your function on the current element
#     result <- bjornerass_screening(dat2$animal_id, 
#                                    dat2$long_5070, 
#                                    dat2$lat_5070,
#                                    dat2$date_time_utc, delta, mu, alpha, theta)
#     
#     # Further processing or analysis
#     # ...
#   }, error = function(err) {
#     # Print or store information about the element causing the error
#     print(paste("Error in element:", animal_ids[i]))
#   })
# }






