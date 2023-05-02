#### OCP iSSF Module_02: Prepare Screened GPS Data in amt ####

# Author: Read Barbee

# Date:2023-04-28 

# Purpose: Prepare screened GPS data in amt to fit iSSFs. Note, data are not yet actually screened.


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

############################### Import Location Data #################################
# Mountain lion location data (November 2022)
locs_raw <- read_csv("data/Location Data/Source Files/locations_master/all_locations_trimmed_2022-11-06.csv")

# Mountain lion deployments  (September 2022)
deployments <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/OCP_Cougar_Deployments_9-30-22.csv")

# Dispersal inventory from Teams (3-13-2023)
dispersals <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Dispersals_03-13-2023.csv")

############################### Import Covariate Data #################################

cov_stack <- rast("data/homerange_habitat_layers_JR_9-14-22/cov_stack1_jr_4-28-2023.tif")

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
         dispersal_status) %>% 
  rename(animal_id=name)


#remove missing locations and convert GMT timestamp to local time
locs <- locs_raw %>% 
  mutate(date_time_local = with_tz(date_time_gmt,"US/Pacific"), 
         .after=date_time_gmt) %>% 
  filter(!is.na(latitude))

#nest locations into list columns by animal_id
locs_nested <- locs %>% nest_by(animal_id)



#add columns for sex and dispersal status to nested locations 
locs_nested <- locs_nested %>% 
  left_join(dem_cats, by = join_by(animal_id))

#reorder columns
locs_nested <- locs_nested %>% 
  select(animal_id, sex, dispersal_status, data)


################################ Convert data to amt tracks #################################

#Create track for each animal

multi_track <- function(d){
  make_track(d, longitude, latitude, date_time_local, crs = 4326) %>% 
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


################################ Resample Tracks/Generate Steps/Extract covariate values #################################

#function to convert locatons to 6hr steps and extract covariate values
steps_6h <- function(x) {
  x %>% 
    amt::track_resample(rate = hours(6), tolerance = minutes(15)) %>% #resample to 2 hours
    amt::filter_min_n_burst() %>% #divide into bursts of min 3 pts
    amt::steps_by_burst() %>%  #convert to steps
    amt::random_steps() %>% #generate 10 random steps per used step
    amt::extract_covariates(cov_stack, where = "both") %>% #extract covariates at start and end
    #amt::time_of_day(include.crepuscule = FALSE) %>% #calculate time of day for each step--not working
    mutate(unique_step = paste(burst_,step_id_,sep="_")) #add column for unique step id 
}

#remove cato and the other individuals not compatible with the 6 hour sampling interval
amt_locs <- locs_nested %>% 
  filter(!(animal_id %in% indiv_to_remove)) %>% 
  filter(animal_id != "Cato")

#map the step function over each individual and append as a nested dataframe
amt_locs$steps <- map(amt_locs$tracks, steps_6h)

#select the relevant columns and unnest the steps column
amt_steps <- amt_locs %>% 
  select(animal_id:dispersal_status,
         steps) %>% 
  unnest(cols=c(steps))

#scale covariates (optional)
amt_steps_scaled <- amt_steps %>%
  mutate(across(elev_start:landuse_hii_end, scale)) %>%
  mutate(across(elev_start:landuse_hii_end, as.numeric))



# amt_steps_scaled %>% unnest(cols=steps) %>% 
#write_csv(amt_steps_scaled, "6h_steps_scaled_cov_4-28-2023.csv")




