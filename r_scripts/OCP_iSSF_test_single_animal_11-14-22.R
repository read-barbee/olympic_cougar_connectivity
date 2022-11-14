### OCP iSSF Test Script -- Single Animal ###

# Read Barbee
# November 14, 2022


#Load libraries
library(tidyverse)
library(terra)
library(lubridate)
library(amt)




# Import all location data
locs_raw <- read_csv("data/Location Data/Source Files/locations_master/all_locations_trimmed_2022-11-06.csv")

#remove missing locations
locs <- locs_raw %>% 
  filter(!is.na(latitude))


#test just with Al
locs_al <- locs %>% 
  filter(animal_id=="Al")

#make track, specify initial crs as WGS84: 4326; transform to NAD83 / Conus Albers: 5070
al_track <- make_track(locs_al, 
                     longitude,
                     latitude,
                     date_time_gmt,
                     all_cols=TRUE,
                     crs = 4326) %>% 
          transform_coords(crs_to = 5070)


#sumarize distribution of time intervals between locations
summarize_sampling_rate(al_track)

#resample track to standarize to 2 hour sampling rate with tolerance of 2 minutes, filter to bursts with minimum of 3 successive fixes to calculate turn angles, and convert from point to step representation with step lengths and turn angles

al_steps <- track_resample(al_track, 
               rate = hours(2),
               tolerance = minutes(2)) %>% 
  filter_min_n_burst(min_n = 3) %>% 
  steps_by_burst()

#check structure of al_steps
str(al_steps)


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

#extract covariate values for the beginning of each step--Note: requires raster class, not SpatRaster



 
full_dat <- al_steps %>% random_steps(n = 10) %>% #generate 10 random steps for every actual step
  extract_covariates(raster(elev_reproj), #extract covariate values at end of each step
                     where = "end") %>% 
  extract_covariates(raster(ndvi_reproj),
                     where = "end") %>%
  mutate(log_sl_ = log(sl_)) %>% # add column for log of step length
  rename(elev = puma_elev_op, #rename covariates for modeling
         ndvi = puma_ndvi_op)

#Fit iSSF
m1 <- full_dat %>% 
  fit_issf(case_ ~ elev + ndvi + log_sl_ + strata(step_id_))
  



