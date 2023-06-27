#### Addressing Habitat Bias and GPS Error with AniMotum####

# Author: Read Barbee

# Date:2023-06-23 

# Purpose: Try out the new aniMotum package to fit continuous time movement model to biased and error-heavy GPS data


###############################################################################


#### Libraries ####
library(tidyverse)
library(aniMotum)
library(sf)
library(pathroutr)


##################################################################
##
## 1. Import and format location data 
##
##################################################################

#all location data
data <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_6-23-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = with_tz(date_time_local, tzone="US/Pacific"))

#subset location data to al as test. Round to nearest hour (only removes one minute)
al_dat <- data %>% 
  filter(animal_id =="Al") %>% 
  mutate(date_time_utc = round_date(date_time_utc, unit = "hour"),
         date_time_local = round_date(date_time_local, unit = "hour"),
         lc = "G") %>% 
  na.omit()

#convert Al's dataframe to sf object. Reproject to WGS 84 Pseudomercator to match polygon file. Recast as Multipoint
al_sf <- al_dat %>% st_as_sf(coords= c("lon_wgs84", "lat_wgs84"), crs = 4326) %>% 
  select(animal_id, date_time_utc, lon_utm, lat_utm) %>% 
  rename(id = animal_id, date = date_time_utc, x=lon_utm, y=lat_utm) %>% 
  mutate(type = "original") %>% 
  st_transform(crs = 3857) %>% 
  st_cast("MULTIPOINT")

##################################################################
##
## 2. Import and format barrier polygons
##
##################################################################

#Import water body polygons
water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp") 

water_polys_filtered <- water_polys %>% filter(WB_PERIOD_ =="PER" ) #|WB_PERIOD_ =="INT"

#Filter water body polygons to retain those within convex hull of all of Al's points
water_polys_cropped <- sf::st_buffer(al_sf, dist = 10000) %>%
  sf::st_union() %>%
  sf::st_convex_hull() %>%
  sf::st_intersection(water_polys_filtered) %>%
  st_collection_extract('POLYGON') %>%
  st_union() %>% 
  st_sf()

#view the water polygon features retained
mapview::mapview(water_polys_cropped)

##################################################################
##
## 3. Fit state-space model to location data
##
##################################################################

#Fit state-space model with velocity fiilter and automatic error filtering
#can also specify angles and max distance for outlier spikes
fit<- fit_ssm(al_dat, 
                vmax = 20, 
                model = "crw", 
                time.step = 2, 
                control = ssm_control(verbose = 0),
                id = "animal_id", 
                date = "date_time_utc", 
                lc = "lc", 
                coord = c("lon_wgs84", "lat_wgs84"), 
                tz = "UTC")

#extract fitted (irregular time) and predicted (regular time) locations from model
fitted <-  grab(fit, what = "fitted")
predicted <- grab(fit, what = "predicted")


#convert fitted and predicted locations to sf objects
fitted_sf <- fitted %>% st_as_sf(coords= c("lon", "lat"), crs = 4326) %>% 
  select(id:y) %>% 
  mutate(type="fitted")

predicted_sf <- predicted %>% st_as_sf(coords= c("lon", "lat"), crs = 4326) %>% 
  select(id:y) %>% 
  mutate(type = "predicted")

predicted_sf2 <- predicted_sf %>% 
  st_transform(crs = 3857) %>% 
  select(-type) %>% 
  st_cast("MULTIPOINT") %>% 
  summarise(do_union=FALSE) %>% 
  st_geometry()

##################################################################
##
## 3. Re-route path so predicted points don't end up in the middle of water bodies
##
##################################################################

#create buffer around barrier objects as visgraph for rerouting function (connects all verticies of barrier polygon with Delaunay triangle mesh and removes any edges that cross the barrier). Essentially it creates a roadmap of traversible terrain.
visgraph <- pathroutr::prt_visgraph(water_polys_cropped, buffer = 15)

#Create table of all consecutive track points that intersect with a barrier polygon and calculate the shortest path through the visibility network between the non-intersecting points.
# segs_tbl <- pathroutr::get_barrier_segments(predicted_sf2, water_polys_cropped) %>% 
#   prt_shortpath(visgraph, blend=TRUE)
# This is simplified by using the prt_reroute() funciton below:

#Reroute the path based on the visibility network
rerouted <- pathroutr::prt_reroute(predicted_sf2, water_polys_cropped, visgraph, blend = TRUE) %>% 
  pathroutr::prt_update_points(predicted_sf2)


##################################################################
##
## 4.Compare the predicted locations with the observed locations
##
##################################################################


#compare original and fitted locations
og_fitted <- rbind(al_sf, fitted_sf)
fit_pred <- rbind(predicted_sf, fitted_sf )
og_pred <- rbind(al_sf, predicted_sf)

rr <- rbind(al_sf, rerouted %>% select(-fid))

mapview::mapview(rr, zcol="type")


mapview::mapview(og_fitted, zcol="type")
mapview::mapview(fit_pred, zcol="type")
mapview::mapview(og_pred, zcol="type")


#Error check visualizations

# ggplot() + 
#   ggspatial::annotation_spatial(water_polys_cropped, fill = "cornsilk3", size = 0) +
#   ggspatial::layer_spatial(predicted_sf2) +
#   theme_void()
# 
# ggplot() + 
#   ggspatial::annotation_spatial(water_polys_cropped, fill = "cornsilk3", size = 0) +
#   ggspatial::layer_spatial(segs_tbl$geometry, color = "deepskyblue3") +
#   theme_void()
# 
# ggplot() + 
#   ggspatial::annotation_spatial(data = water_polys_cropped, 
#                                 fill = "cornsilk3", 
#                                 size = 0) +
#   ggspatial::layer_spatial(data = al_sf) +
#   theme_void()


###############################################################################  