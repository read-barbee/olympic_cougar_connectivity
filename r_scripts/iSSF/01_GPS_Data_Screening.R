#### OCP iSSF Module_01: GPS Data Screening ####

# Author: Read Barbee

# Date:2023-05-05

# Inputs:
#   •	GPS0: Raw GPS locations
# 
# Outputs:
#   •	GPS1: Screened GPS locations
# 
# Steps- TBD
# •	Screen for error locations 
# o	state-space model?
#   o	DOP?
#   •	Screen for habitat bias




################################ Libraries #################################
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)

################################ User-Defined Parameters #################################


############################### Import Location Data #################################
# Mountain lion location data (November 2022)
locs_raw <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_5-16-2023.csv")

# Mountain lion deployments  (September 2022)
deployments <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/collar_deployments_master_5-11-2023.csv")

# Dispersal inventory from Teams (3-13-2023)
dispersals <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/dispersals_master_5-11-2023.csv")


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


#Look at DOP scores: missing scores for 24,145 locations. Missing for the collar downloads?
summary(locs_raw)

locs_raw %>% 
  filter(is.na(latitude)==FALSE) %>% 
  filter(is.na(dop)==TRUE) %>% 
  distinct(animal_id)

locs_raw %>% 
  filter(is.na(latitude)==FALSE) %>% 
  filter(dop<=5 | is.na(dop)==TRUE)



######## Remove first 24 hours for each animal to avoid capture effects. Remove 2D fixes with dop > 5 (optional, kind of removes a lot of data)


locs_screened <- locs_raw%>% 
  group_by(animal_id) %>% 
  filter(date_time_gmt >= (min(date_time_gmt) + hours(24))) %>% 
  ungroup() %>% 
  filter(!is.na(latitude)) %>% 
  filter(!(dop>5 & fix_type %in% c("2-D least-squares", "2D Fix", "2D"))) %>% 
  mutate(unique_id= 1:nrow(.), .before = deployment_id) %>% 
  filter(!(unique_id %in% c(44803, 214496, 182176, 132748, 213380, 78225)))

#compare to minimum date times from original data frame to make sure it worked
# locs_raw%>% 
#   group_by(animal_id) %>% 
#   summarize(min=min(date_time_gmt))

#Check mapview for locations outside of study area
locs_screened %>% sf::st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>% mapview::mapview()

write_csv(locs_screened, "data/Location Data/Source Files/locations_master/locations_screened_5-9-2023.csv")

