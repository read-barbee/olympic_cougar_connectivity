#   Messing around in CTMM                                                         ####

# Author: Read Barbee

# Date:2022-11-06 

library(tidyverse)
library(ctmm)
library(janitor)
library(lubridate)

#import data
cougars_raw <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/locations_master/all_locations_trimmed_2022-11-06.csv")


#make sure all times are in order by individual--just creates more errors so far
# cougars_raw <- cougars_raw %>% 
#   group_by(deployment_id) %>% 
#   arrange(date_time_gmt, .by_group = TRUE) %>% 
#   ungroup()

# cougars_raw %>% 
#   filter(animal_id=="Elwha")

#get into movebank format
cougars_mb <- cougars_raw %>% 
  select(animal_id,
         collar_id,
         date_time_gmt,
         latitude,
         longitude) %>% 
  rename(individual.local.identifier = animal_id,
         tag.local.identifier = collar_id,
         timestamp = date_time_gmt,
         location.long = longitude,
         location.lat = latitude)
 

#check for duplicate time stamps for each individual
get_dupes(cougars_raw, animal_id, date_time_gmt)

#import to ctmm
cougars <- as.telemetry(cougars_mb)
#Warnings: look into duplicate times and order for Elwha and Hoko


# check class
class(cougars)

# check number of datasets
length(cougars)

# names of cougars
names(cougars)

# summary of cougar data-- really nice output***
summary(cougars)


# plot cougars with spatially-separated rainbow of colors
COL <- color(cougars,by='individual')
plot(cougars,col=COL)
