### OCP iSSF Test Script -- All Animals ###

# Read Barbee
# November 14, 2022


################################ Libraries #################################
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)

################################ User-Defined Parameters #################################

project_crs <- "epsg:5070"

resample_int_hours <- 6

tolerance <- minutes(15)

rand_steps <- 10


################################ Data Import #################################
# Mountain lion location data (November 2022)
locs_raw <- read_csv("data/Location Data/Source Files/locations_master/all_locations_trimmed_2022-11-06.csv")

# Mountain lion deployments  (September 2022)
deployments <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/OCP_Cougar_Deployments_9-30-22.csv")

# Dispersal inventory from Teams (3-13-2023)
dispersals <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Dispersals_03-13-2023.csv")


################################ Data Cleaning/Formatting #################################


#Filter deploment list to get one row per individual
first_deps <- deployments %>% 
  clean_names() %>% 
  distinct(name, .keep_all = TRUE) %>% 
  mutate(name = trimws(name))

#Get vector of disperser names
disperser_names <- dispersals %>% 
  clean_names() %>% 
  filter(!is.na(name)) %>% 
  distinct(name) %>% 
  pull(name)


#Add dispersal status to deployment list
dem_cats <- first_deps %>% 
  mutate(dispersal_status = case_when(name %in% disperser_names == TRUE ~ "disperser",
                                      name %in% disperser_names == FALSE ~ "resident")) %>% 
  select(name,
         sex,
         dispersal_status)


#remove missing locations
locs <- locs_raw %>% 
  filter(!is.na(latitude))

#nest locations into list columns by animal_id
locs_nested <- locs %>% nest_by(animal_id)

#add columns for sex and dispersal status to nested locations 
locs_nested$sex <-  dem_cats$sex
locs_nested$dispersal_status <-  dem_cats$dispersal_status

#reorder columns
locs_nested <- locs_nested %>% 
  select(animal_id, sex, dispersal_status, data)


################################ Convert data to amt tracks #################################

#Create track for each animal

multi_track <- function(d){
  make_track(d, longitude, latitude, date_time_gmt, crs = 4326) %>% 
    transform_coords(crs_to = 5070)
}

locs_nested$tracks <- map(locs_nested$data, multi_track)


################################ Examine Sampling Rates #################################

#inspect sampling rate for each individual

locs_nested$sr <-  map(locs_nested$tracks, summarize_sampling_rate)

#sampling rate statistics by individual
sampling_rates <- locs_nested %>% 
  select(-c(data, tracks)) %>% 
  unnest(cols = c(sr))

summary(sampling_rates)

#the largest minimum fix interval across all individuals is 4 hours


################################ Identify optimum resampling interval #################################

#function to calculate bursts for each individual---increment by hours to find optimum
test_fun1 <- function(x) {
  x %>% 
    amt::track_resample(rate = hours(6), tolerance = minutes(15)) %>% #resample to 2 hours
    amt::filter_min_n_burst() #divide into bursts of min 3 pts
}

locs_nested_test <- locs_nested

locs_nested_test$steps <- map(locs_nested$tracks, test_fun1)


locs_nested_test %>% 
  filter(nrow(steps)==0) %>% 
  select(steps) %>% View()

indiv_to_remove <- locs_nested_test %>% 
  filter(nrow(steps)==0) %>% 
  pull(animal_id)


#1 hour removes 69 individuals (n=40)
#2 hour removes 32 individuals (n=77)
#3 hours removes 32 individuals (n=77)
#4 hours removes 33 individuals (n=76)
#5 hours removes 63 individuals (n=46)
#6 hours removes 8 individuals (n=101)

#6 hours retains the most individuals


################################ Resample Tracks/Generate Steps #################################

steps_6h <- function(x) {
  x %>% 
    amt::track_resample(rate = hours(6), tolerance = minutes(15)) %>% #resample to 2 hours
    amt::filter_min_n_burst() %>% #divide into bursts of min 3 pts
    amt::steps_by_burst() %>% 
    amt::random_steps() %>% #generate 10 random steps per used step
    mutate(unique_step = paste(burst_,step_id_,sep="_")) #add column for unique step id 
}

amt_locs <- locs_nested %>% 
  filter(!(animal_id %in% indiv_to_remove)) %>% 
  filter(animal_id != "Cato")

amt_locs$steps <- map(amt_locs$tracks, steps_6h)

amt_steps <- amt_locs %>% 
  select(animal_id:dispersal_status,
         steps) %>% 
  unnest(cols=c(steps))


#### By Min: Retains 100 individuals for 6 hour sampling: confirms Tom's method #####
# 
# #filter sampling rate summary frame to get vector of individual names to resample
# indiv_for_resampling_min <- sampling_rates %>% 
#   mutate(min=round(min, digits = 0)) %>% 
#   filter(min < 4) %>% 
#   pull(animal_id)
# 
# #hist(indiv_for_resampling_min$min)
# 
# #filter locations based on sampling rate critera. Remove additional problematic individuals  
# test_min <-  locs_nested %>%
#   filter(animal_id %in% indiv_for_resampling_min) %>% 
#   filter(!(animal_id %in% c("Butch", "Cato", "Crash", "OtookTom", "Promise"))) %>% 
#   select(-c(data, sr)) 
# 
# 
# #resample tracks and generate steps
# test_min$steps <-  map(test_min$tracks, steps_6h)

#### By Mean: Only retains 63 individuals ##
# 
# #filter sampling rate summary frame to get vector of individual names to resample
# indiv_for_resampling <- sampling_rates %>% 
#   mutate(mean=round(mean, digits = 0)) %>% 
#   filter(mean %in% c(1, 2, 3, 6)) %>% 
#   pull(animal_id)
# 
# #hist(indiv_for_resampling$mean)
# 
# #filter locations based on sampling rate critera. Remove additional problematic individuals  
# test <-  locs_nested %>%
#   filter(animal_id %in% indiv_for_resampling) %>% 
#   filter(!(animal_id %in% c("Cato", "Crash"))) %>% 
#   select(-c(data, sr)) 
# 
# 
# #resample tracks and generate steps
# test$steps <-  map(test$tracks, steps_6h)











################################ Import Covariates #################################

#all layers are in EPSG: 4326 with resolution (0.0002694946, 0.0002694946)
elev <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/assets/puma_elev_op.tif")

ndvi <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/assets/puma_ndvi_op.tif")

dist_water <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/assets/puma_waterDist_op.tif")

roads_hii <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/assets/puma_hii_roads_op.tif")

forest <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/assets/puma_forest_op.tif")

landuse_hii <- rast("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/homerange_habitat_layers_JR_9-14-22/assets/puma_hii_landuse_op.tif")

rast_reproj <- function(raster){
  project(x=raster, y=project_crs )
}

#very slow
#map(list(elev, ndvi, dist_water, roads_hii, forest, landuse_hii), rast_reproj)


#reproject habitat layers to NAD83 / Conus Albers: 5070

elev_reproj <- project(x=elev,
                       y= "epsg:5070")

ndvi_reproj <- project(x=ndvi,
                       y= "epsg:5070")

stack <- c(elev_reproj,
           ndvi_reproj)



