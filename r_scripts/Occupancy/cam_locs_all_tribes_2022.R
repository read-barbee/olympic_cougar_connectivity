#### Inspecting full 2022 Camera Grid ####

# Author: Read Barbee

# Date:2023-02-10 

# Purpose:


###############################################################################
#### Library / Functions / Data ####

#### Library ####
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggmap)
library(basemaps)

#### Functions ####

#### Data ####

pnptc <- st_read("data/Camera Data/2022/2022_PNPTC/2022_PlacedCameras_PNPTC.shp")
lekt <- st_read("data/Camera Data/2022/2022_LEKT/2022.LEKT.CamGrid.shp")
makah <- st_read("data/Camera Data/2022/2022_Makah/Makah2022_CamGrid_Corrected.shp")
skok <- st_read("data/Camera Data/2022/2022_Skok/SkokStations2022.shp")

#not reading in correctly
q_res <- st_read("data/Camera Data/2022/2022_Quin_Res/Quin.Res.Grid.shp")
q_wyn <- st_read("data/Camera Data/2022/2022_Quin_Wyno/2022.Quin.Wynochee.Grid.shp")


#### Write corrected makah shape file (corrected erroneous lat/long values) #####
# 
# makah <- st_read("data/Camera Data/2022/2022_Makah/Makah2022_CamGrid_Final.shp")
# 
# makah_corrected <- as.data.frame(makah) %>%
#   mutate(Longitude = ifelse(Longitude > 0, Longitude*-1, Longitude)) %>%
#   mutate(Latitude = ifelse(round(Latitude)==43,
#                            Latitude + 5,
#                            Latitude)) %>%
#   select(-geometry) %>%
#   st_as_sf(coords = c("Longitude", "Latitude"),
#            crs=4326,
#            remove=FALSE) %>%
#   rename(station=Station,
#          general_loc=General_Lo,
#          lat=Latitude,
#          lon=Longitude)
# 
# st_write(makah_corrected, "data/Camera Data/2022/2022_Makah/Makah2022_CamGrid_Corrected.shp")

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


#Combine all grids into single dataframe
all_cams <- rbind(lekt, makah, pnptc, skok, q_res, q_wyn)


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
  scale_color_hue(name="Camera Grid", labels=c("Lower Elwha", "Makah", "PNPTC", "Quinault-Reservation","Quinault-Wynoochee", "Skokomish")) +
  labs(title="Olympic Cougar Project Camera Sites 2022 (n = 518)")

# ggsave("data/Camera Data/2022/camera_locations_2022_v2.png", cam_plot2, dpi="print")




