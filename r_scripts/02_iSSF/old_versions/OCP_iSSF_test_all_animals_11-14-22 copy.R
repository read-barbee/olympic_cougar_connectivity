### OCP iSSF Test Script -- All Animals ###

# Read Barbee
# November 14, 2022


#Load libraries
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)



# Import all location data
locs_raw <- read_csv("data/Location Data/Source Files/locations_master/all_locations_trimmed_2022-11-06.csv")

#import deployments for individual-specific metadata
deployments <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/OCP_Cougar_Deployments_9-30-22.csv")

#filter out duplicate rows for individuals
first_deps <- deployments %>% 
filter(!duplicated(Name))


#remove missing locations
locs <- locs_raw %>% 
  filter(!is.na(latitude))

#nest locations into list columns by animal_id
locs_nested <- locs %>% nest_by(animal_id)

#add column for sex from deployments dataframe   
locs_nested$sex <-  first_deps$Sex

#join <- left_join(locs_nested, first_deps, by=("animal_id" = "Name")) join or merge better than direct assign

#Create track for each animal

multi_track <- function(d){
  make_track(d, longitude, latitude, date_time_gmt, crs = 4326) %>% 
    transform_coords(crs_to = 5070)
}

locs_nested$tracks <- map(locs_nested$data, multi_track)


#inspect sampling rate for each individual
locs_nested$sampling_rate <- map(locs_nested$tracks, summarize_sampling_rate)

sampling_rates <- locs_nested %>% 
select(animal_id, sampling_rate) %>% 
  unnest(cols = c(sampling_rate))

range(sampling_rates$median)

#filter out individuals with sampling rate >= 4 hours (n = 4)
keepers <- sampling_rates %>% 
  filter(median < 4)

locs_nested %>% 
  filter(animal_id %in% keepers$animal_id)


### Import habitat covariate data

elev <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/puma_elev_op.tif")

ndvi <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/puma_ndvi_op.tif")


#reproject habitat layers to NAD83 / Conus Albers: 5070

elev_reproj <- project(x=elev,
                       y= "epsg:5070")

ndvi_reproj <- project(x=ndvi,
                       y= "epsg:5070")

stack <- c(elev_reproj,
           ndvi_reproj)


#function for resampling and fitting iSSFs for each individual: incomplete--not tested yet
model_fun <- function(x){ 
  x %>% 
    track_resample(rate = hours(3),
                   tolerance = minutes(10)) %>% 
    filter_min_n_burst(min_n = 3) %>% 
    steps_by_burst() %>% 
    random_steps(n = 10)
    
}


deployments2 <- deployments %>% 
  clean_names() %>% 
  filter(!duplicated(name)) %>% 
  mutate(deployment_date = mdy(deployment_date),
         end_date = mdy(end_date))




#dispersal inventory from Teams 3-13-2023
dispersals <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Dispersals_03-13-2023.csv")

#count dispersers
dispersals %>% 
  clean_names() %>% 
  filter(!is.na(name)) %>% 
  distinct(name) %>% 
  count()

dispersals %>% 
  clean_names() %>% 
  filter(!is.na(name)) %>% 
  group_by(sex) %>% 
  distinct(name) %>% 
  count()

#Resident Males: 53
#Disperser Males: 22

#Resident Females: 56
#Disperser Females: 13


