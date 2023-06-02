#### OCP iSSF Module_02: Prepare Screened GPS Data in amt ####

# Author: Read Barbee

# Date:2023-06-02 

# Purpose: Prepare screened GPS data in amt to fit iSSFs.

# Inputs:
#   •	Screened GPS Locations (corrected for error and habitat bias)
#   •	Deployment list
#   •	Disperser list
#   •	Covariate Stack
#
# Outputs:
#   •	amt_step dataframe with fields for all necessary covariates

#Steps
# •	Import location and covariate data
# •	Set time zones
# •	Formatting: Nest locations by individual and join columns for sex and dispersal status 
# •	Generate amt tracks for each individual
# •	Examine sampling rates and identify resampling interval to retain the most individuals
# •	Remove individuals not compatible with that resampling interval
# •	Resample tracks, generate used and random steps, and extract covariate values for end of each step
# •	Filter steps to remove outliers and stationary steps (100 m < sl_ < 20,000 m )
# •	Remove individuals with < 30steps and/or < 30 days of steps????***
# •	Transform aspect to continuous easting and northing values
# •	Define USFS landcover and landuse categories
# •	Derive covariates for climatic season, hunting season, and calving season
# •	Round all floating hii values to integers
# •	reorder columns and export to csv

#

################################ Libraries #################################
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)
library(DataExplorer)
library(future)
# library(sf)
# library(mapview)

################################ User-Defined Parameters #################################

project_crs <- "epsg:5070" #NAD83 Albers Equal Area Projection

resample_int_hours <- 6

tolerance_mins <- 15

rand_steps <- 20

############################### Import Location Data #################################
# Mountain lion location data (May 2023)
locs_raw <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_6-1-2023.csv", col_types = list(fix_type = col_character()))

# Mountain lion deployments  (September 2022)
deployments <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/collar_deployments_master_5-11-2023.csv")

# Dispersal inventory from Teams (3-13-2023)
dispersals <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/dispersals_master_5-11-2023.csv")

############################### Import Covariate Data #################################

cov_stack <- rast("data/Habitat_Covariates/puma_cov_stack_v1/cov_stack1.tif")

#rename covariate bands
names(cov_stack) <- c("tree_cover_hansen",
                      "gpp",
                      "infra_hii",
                      "landuse_hii",
                      "land_cover_usfs",
                      "land_use_usfs",
                      "npp",
                      "popdens_hii",
                      "power_hii",
                      "precip",
                      "rails_hii",
                      "roads_hii",
                      "elevation",
                      "slope",
                      "aspect",
                      "tri",
                      "tpi",
                      "perc_tree_cover",
                      "perc_nontree_veg",
                      "perc_nonveg",
                      "ndvi",
                      "evi",
                      "dist_water")

#inspect covariate stack--plotting takes a long time. Parallel doesn't seem to be working
# plan(multicore)
# terra::plot(cov_stack, parallel = TRUE)

################################ Data Cleaning/Formatting #################################

#Make sure R recognizes the local time column in the correct time zone
tz(locs_raw$date_time_local) <- "US/Pacific"

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
         dispersal_status) %>% 
  rename(animal_id=name)


#Nest locations by animal_id and add columns for sex and dispersal status
locs_nested <- locs_raw %>% 
  nest_by(animal_id) %>% 
  left_join(dem_cats, by = join_by(animal_id)) %>% 
  select(animal_id, sex, dispersal_status, data)


################################ Convert data to amt tracks #################################

#Create track for each animal

multi_track <- function(d){
  make_track(d, longitude, latitude, date_time_local, crs = 4326) %>% 
    transform_coords(crs_to = 5070)
}

locs_nested$tracks <- map(locs_nested$data, multi_track)


#subset only disperser locations
dispersers <- locs_nested %>% 
  filter(dispersal_status=="disperser")


################################ Examine Sampling Rates #################################

#inspect sampling rate for each individual

locs_nested$sr <-  map(locs_nested$tracks, summarize_sampling_rate)

#sampling rate statistics by individual
sampling_rates <- locs_nested %>% 
  select(-c(data, tracks)) %>% 
  unnest(cols = c(sr))

summary(sampling_rates)


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
  select(steps)

#make list of individuals to remove (because they can't be resampled at 6 hours)
indiv_to_remove <- locs_nested_test %>% 
  filter(nrow(steps)==0) %>% 
  pull(animal_id)


#1 hour removes 67 individuals 
#2 hour removes 31 individuals 
#3 hours removes 29 individuals 
#4 hours removes 32 individuals 
#5 hours removes 64 individuals 
#6 hours removes 10 individuals 

#6 hours retains the most individuals


################################ Resample Tracks/Generate Steps/Extract covariate values #################################

#function to convert locatons to 6hr steps and extract covariate values
steps_calc <- function(x) {
  x %>% 
    amt::track_resample(rate = hours(resample_int_hours), tolerance = minutes(tolerance_mins)) %>% #resample to 2 hours
    amt::filter_min_n_burst() %>% #divide into bursts of min 3 pts
    amt::steps_by_burst() %>%  #convert to steps
    amt::random_steps(n_control = rand_steps) %>% #generate random steps per used step
    amt::extract_covariates(cov_stack, where = "end") %>% #extract covariates at start and end
    #amt::time_of_day(include.crepuscule = FALSE) %>% #calculate time of day for each step--not working
    mutate(unique_step = paste(burst_,step_id_,sep="_")) #add column for unique step id 
}

#remove individuals not compatible with the 6 hour sampling interval
amt_locs <- locs_nested %>% 
  filter(!(animal_id %in% indiv_to_remove)) 


#map the step function over each individual and append as a nested dataframe
amt_locs$steps <- map(amt_locs$tracks, steps_calc)

#select the relevant columns and unnest the steps column
amt_steps <- amt_locs %>% 
  select(animal_id:dispersal_status,
         steps) %>% 
  unnest(cols=c(steps)) %>% 
  relocate(unique_step, .after = step_id_) %>% 
  ungroup()

################################ Secondary Screening #################################

#Check global step length and turn angle distributions
#don't seem to be moving more than 10 km in 6 hours
amt_steps %>% 
  #filter(case_==TRUE) %>% 
  filter(sl_ >100 & sl_<20000 ) %>% 
  pull(sl_) %>%
  hist()

hist(amt_steps$ta_)

#Remove individuals with < 30steps and/or < 30 days of steps: removes 4 individuals. not sure yet if this is necessary

# removal_list <- amt_steps %>%
#   group_by(animal_id) %>%
#   mutate(date_range=interval(start=min(t1_), end=max(t2_)),
#          step_days = as.duration(date_range)/ddays(1)) %>% 
#   distinct(unique_step, .keep_all = TRUE) %>% 
#   summarize(n=n(), step_days = round(first(step_days), 0)) %>%
#   filter(n<30 & step_days <= 30) %>% 
# pull(animal_id)
# 
# amt_steps <- amt_steps %>% 
#   filter(!(animal_id %in% removal_list))



#filter steps to remove outliers and stationary steps
steps_unscaled <- amt_steps %>% 
  filter(sl_ >100 & sl_<20000 ) %>% 
  terra::na.omit()


#make sure no data are missing after omitting steps with missing raster values
#plot_missing(steps_unscaled)


################################ Covaraite Transformations #################################

### Aspect #####

#convert aspect from degrees to radians and calculate northing and easting variables with cos and sin transformations respectively
steps_unscaled <- steps_unscaled %>% 
  mutate(aspect_rad = (pi*aspect)/180, .after=aspect) %>%
  mutate(northing = cos(aspect_rad),
         easting = sin(aspect_rad), .after=aspect_rad) %>% 
  rename(aspect_deg = aspect)
  

### USFS Land Cover #####

#round mean values for landuse and landcover to nearest integer
steps_unscaled <- steps_unscaled %>% 
  mutate(land_cover_usfs = round(land_cover_usfs),
         land_use_usfs = round(land_use_usfs))

#define land cover categories
steps_unscaled <- steps_unscaled %>% 
  mutate(land_cover_usfs = case_when(land_cover_usfs == 1 ~ "trees",
                                     land_cover_usfs == 2 ~ "tall_trees_shrubs",
                                     land_cover_usfs == 3 ~ "tree_shrub_mix",
                                     land_cover_usfs == 4 ~ "gfh_tree_mix",
                                     land_cover_usfs == 5 ~ "barren_tree_mix",
                                     land_cover_usfs == 6 ~ "tall_shrubs",
                                     land_cover_usfs == 7 ~ "shrubs",
                                     land_cover_usfs == 8 ~ "gfh_shrub_mix",
                                     land_cover_usfs == 9 ~ "barren_shrub_mix",
                                     land_cover_usfs == 10 ~ "gfh",
                                     land_cover_usfs == 11 ~ "barren_gfh_mix",
                                     land_cover_usfs == 12 ~ "barren_impervious",
                                     land_cover_usfs == 13 ~ "snow_ice",
                                     land_cover_usfs == 14 ~ "water",
                                     land_cover_usfs == 15 ~ NA_character_))


### USFS Land Use #####

#define land use categories
steps_unscaled <- steps_unscaled %>% 
  mutate(land_use_usfs = case_when(land_use_usfs == 1 ~ "agriculture",
                                   land_use_usfs == 2 ~ "developed",
                                   land_use_usfs == 3 ~ "forest",
                                   land_use_usfs == 4 ~ "non_forest_wetland",
                                   land_use_usfs == 5 ~ "other",
                                   land_use_usfs == 6 ~ "rangeland_pasture",
                                   land_use_usfs == 7 ~ NA_character_))


################################ Add Season Covariates #################################

#hunting season (deer and elk): Sep 15 - Nov 15 (WDFW)
#calving season (deer and elk): May 15 - July 1 (WDFW)
#wet season: October - April (en.climate data.org)
#dry season: May - September (en.climate data.org)

steps_unscaled <- steps_unscaled %>% 
  mutate(season = case_when(month(t2_) >= 10 | month(t2_) <= 4 ~ "wet",
                            month(t2_) < 10 & month(t2_) > 4 ~ "dry"),
         hunting_season = ifelse(yday(t2_) %in% c(258:319), "yes", "no"),
         calving_season = ifelse(yday(t2_) %in% c(135:182), "yes", "no"))


# check to make sure categorization worked correctly
# filter(year(t2_)== 2020) %>% 
#   group_by(season) %>% summarize(first = min(t2_), last = max(t2_))



#add unique id and rearrange fields to final format
steps_final <- steps_unscaled %>% 
  mutate(unique_id = 1:nrow(steps_unscaled), .before= animal_id) %>% 
  select(unique_id:burst_,
         step_id_,
         unique_step,
         case_,
         x1_:dt_,
         gpp,
         npp,
         ndvi,
         evi,
         tree_cover_hansen,
         perc_tree_cover:perc_nonveg,
         land_cover_usfs,
         precip,
         dist_water,
         elevation:tpi,
         land_use_usfs,
         roads_hii,
         popdens_hii,
         landuse_hii,
         infra_hii,
         rails_hii,
         power_hii,
         season:calving_season
         ) %>% 
  mutate(across(roads_hii:power_hii, round))


#scale covariates (optional)
# amt_steps_scaled <- amt_steps %>%
#   mutate(across(elev_start:landuse_hii_end, scale)) %>%
#   mutate(across(elev_start:landuse_hii_end, as.numeric))


#write_csv(steps_final, "data/Location_Data/Steps/6h_steps_unscaled_6-02-2023.csv")




