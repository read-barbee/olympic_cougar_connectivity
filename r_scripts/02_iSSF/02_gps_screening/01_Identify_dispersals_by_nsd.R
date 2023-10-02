#### Identify dispersals by nsd ####

# Author: Read Barbee

# Date:2023-10-02

#incomplete


################################ Libraries #################################
library(tidyverse)
library(terra)
library(sf)
library(amt)
library(janitor)

#########################################################################
##
## 1. Import and format location data
##
##########################################################################
# Mountain lion location data (April 2023)

#import screened locations and split into residents and dispersers
locs <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_10-02-2023.csv")


#########################################################################
##
## 7. Identify dispersal events (only need to do this once and save to csv)
##
##########################################################################

# Filter to individual of interest
filt <- locs %>% 
  filter(animal_id=="Junior") 

#view locs before making track
mapview::mapview(filt %>%  sf::st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070))

#make track
track <- make_track(filt, lon_utm, lat_utm, date_time_utc, check_duplicates = TRUE, all_cols = TRUE, crs = 5070) 



#calculate net squared displacement from first location over time
track <- amt::add_nsd(track)

#plot nsd over time to identify dispersal date 
nsd_plot <- ggplot(data = track, 
                   aes(x = t_, y = nsd_)) +
  geom_point()


plotly::ggplotly(nsd_plot)

#recheck Hana, Lady, Lolli
#Lolli's track not trimmed correctly. Also doesn't seem like much of a dispersal
#Harder to determine dispersal date for females because they don't disperse as far

