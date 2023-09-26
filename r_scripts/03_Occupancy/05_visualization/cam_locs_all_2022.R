#### Inspecting full 2022 Camera Grid ####

# Author: Read Barbee

# Date:2023-02-10 

# Purpose: map of camera locations for all tribes and olympic national park


###############################################################################
#### Library / Functions / Data ####

#### Library ####
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggmap)
library(basemaps)
library(janitor)

#### Functions ####

#### Data ####

#PNPTC 2022
pnptc <- st_read("data/Camera Data/2022/2022_PNPTC/2022_PlacedCameras_PNPTC.shp")

#Elwha 2022
lekt <- st_read("data/Camera Data/2022/2022_LEKT/2022.LEKT.CamGrid.shp")

#Makah 2022 corrected for errors
makah <- st_read("data/Camera Data/2022/2022_Makah/Makah2022_CamGrid_Corrected.shp")

#Skokomish 2022
skok <- st_read("data/Camera Data/2022/2022_Skok/SkokStations2022.shp")

#Quinault 2022
q_res <- st_read("data/Camera Data/2022/2022_Quin_Res/Quin.Res.Grid.shp")
q_wyn <- st_read("data/Camera Data/2022/2022_Quin_Wyno/2022.Quin.Wynochee.Grid.shp")

#ONP Fisher Study 2013-2016
onp_raw <- read_csv("data/Camera Data/Olympic_National_Park/ONP_fisher_grid_2013-2016_raw.csv")

#179 hexes sampled 2013-2016, 3 cameras per hex. 537 unique station locations (from progress report)



############################################################################### 


### Clean data for each grid, add column for grid, and convert to crs 3857 (world mercator)

lekt <- lekt %>%
  rename(station=Station_nu,
         lat=Camera_lat,
         lon=Camera_lon,
         general_loc=General_lo) %>%
  mutate(grid="lekt") %>% 
  select(grid, everything(),-general_loc) %>% 
st_transform(3857)


makah <- makah %>% 
  mutate(grid="makah") %>% 
  select(grid, station, lat, lon, geometry) %>% 
  st_transform(3857)

pnptc <- pnptc %>% 
  select(Site, 
         Lat,
         Long,
         geometry) %>% 
  rename(station=Site,
         lat=Lat,
         lon=Long) %>% 
  mutate(grid="pnptc", .before=everything()) %>% 
  st_transform(3857)

skok <- skok %>% 
  rename(station=Site_numbe,
         lat=Latitude,
         lon=Longitude) %>% 
  mutate(grid="skok", .before=everything()) %>% 
  st_transform(3857)

q_res <- q_res %>% 
  rename(station=name) %>% 
  mutate(grid="q_res", .before=everything()) %>% 
  st_transform(3857)
  
q_wyn <- q_wyn %>% 
  select(ident,
         Latitude,
         Longitude,
         geometry) %>% 
  rename(station=ident,
         lat=Latitude,
         lon=Longitude) %>% 
  mutate(grid="q_wyn", .before=everything()) %>% 
  st_zm() %>% 
  st_transform(3857)


#clean ONP data and filter to unique site ids to add to map- 537 unique station locations (179 hexes) over all years, ~9 sampling events/hex
onp_locs <- onp_raw %>% 
  clean_names() %>% 
  filter(study_incid=="O") %>% 
  distinct(utm_e, utm_n, .keep_all=TRUE)



onp_hexes <- onp_raw %>% 
  clean_names() %>% 
  filter(study_incid=="O" & !is.na(hex_id)) %>% 
  distinct(hex_id, .keep_all=TRUE) %>% 
  select(hex_id, utm_e, utm_n) %>%
  mutate(grid="onp") %>% 
  rename(station=hex_id,
         lat=utm_e,
         lon=utm_n) %>% 
    select(grid, everything())

onp_hex_locs <- st_as_sf(onp_hexes, coords = c("lat", "lon"), remove=FALSE, crs=st_crs(32610))

onp_hex_locs_reproj <- st_transform(onp_hex_locs, crs=st_crs(3857))

plot(onp_hex_locs_reproj)

  

#Combine all grids into single dataframe
all_cams <- rbind(lekt, makah, pnptc, skok, q_res, q_wyn, onp_hex_locs_reproj)


#Set mapbox token
mapbox_token <- "pk.eyJ1IjoicmVhZGJhcmJlZSIsImEiOiJjbGUzNzFwNHIwMWlhM3ZtY290YWVrY2RqIn0.glD2Lh70Kyyd-rhaXwc6_Q"


#Check available maptypes for basemap
get_maptypes()

#Set basemap defaults
set_defaults(map_service = "mapbox", map_type = "satellite", map_token=mapbox_token)

#set bounding box for basemap with buffer around extent of all_cams
bbox <- st_bbox(st_buffer(all_cams, dist=20000))

#plot camera points on basemap
# cam_plot1 <- ggplot() + 
#   basemap_gglayer(ext=bbox) +
#   scale_fill_identity() + 
#   coord_sf() +
#   geom_spatvector(data= all_cams, aes(color=as.factor(grid))) + 
#   scale_color_hue(name="Camera Grid", labels=c("Lower Elwha", "Makah", "PNPTC", "Quinault-Reservation","Quinault-Wynoochee", "Skokomish")) +
#   labs(title="Olympic Cougar Project Camera Sites 2022")

#ggsave("data/Camera Data/2022/camera_locations_2022_v1.png", cam_plot1)



set_defaults(map_service = "mapbox", map_type = "light", map_token=mapbox_token)

#set bounding box for basemap with buffer around extent of all_cams

#plot camera points on basemap
cam_plot2 <- ggplot() + 
  basemap_gglayer(ext=bbox) +
  scale_fill_identity() + 
  coord_sf() +
  geom_spatvector(data= all_cams, aes(color=as.factor(grid))) + 
  scale_color_hue(name="Camera Grid", labels=c("Lower Elwha", "Makah", "Olympic National Park", "PNPTC", "Quinault-Reservation","Quinault-Wynoochee", "Skokomish")) +
  labs(title="Olympic Cougar Project Camera Sites 2013-2022 (n = 697)")

#ggsave("data/Camera Data/camera_locations_all_v1.png", cam_plot2, dpi="print")




