#### Measure distances between camera stations of the same name between years ####

# Author: Read Barbee

# Date:2023-09-25

# Purpose:


################################ libraries #################################
library(tidyverse)
library(unmarked)
library(ubms)
library(camtrapR)
library(beepr)
library(doParallel)
library(DataExplorer)
library(stocc)
library(sf)




occ_dat <- read_csv("data/Camera_Data/master/ocp_occ_dat_9-18-23.csv") %>% 
  mutate(across(station_id_year:year, as.factor)) %>% 
  mutate(aspect_rad = (pi*aspect)/180, .after=aspect) %>%
  mutate(northing = cos(aspect_rad),
         easting = sin(aspect_rad), .after=aspect_rad) %>% 
  rename(aspect_deg = aspect) %>% 
  select(-c(aspect_deg, aspect_rad, land_cover_usfs, land_use_usfs))

proj_coords <- occ_dat %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(crs=5070) %>% 
  st_coordinates() %>% as.data.frame

occ_dat <- occ_dat %>% 
  mutate(x = proj_coords[,1],
         y = proj_coords[,2], .after = lat)

#scale covariates
occ_dat_scaled <- occ_dat %>% 
  mutate(across(tree_cover_hansen:dist_water, scale)) %>% 
  mutate(across(tree_cover_hansen:dist_water, as.numeric))


#remove station rows with missing covariate values
complete_cases <- occ_dat_scaled %>% select(tree_cover_hansen:dist_water) %>% complete.cases()

occ_dat_complete <- occ_dat_scaled %>% 
  mutate(comp = complete_cases) %>% 
  filter(comp==TRUE)


#convert to unmarked dataframe
umf <- unmarkedFrameOccu(y = occ_dat_complete %>% select(d_1:d_365),
                         siteCovs = occ_dat_complete %>% select(station_id, tree_cover_hansen:dist_water))

umf_cell <- unmarkedFrameOccu(y = occ_dat_complete %>% select(d_1:d_365),
                              siteCovs = occ_dat_complete %>% select(cell_id, tree_cover_hansen:dist_water))



mean_station_locs <- occ_dat %>% 
  group_by(cell_id) %>% 
  summarize(x = mean(x), y = mean(y)) %>% 
  st_as_sf(coords = c("x", "y"), crs = 5070)


station_locs <- occ_dat %>%  st_as_sf(coords = c("x", "y"), crs = 5070)


#PROBLEM: LEKT Station 73 and others are in wildly different locations between years

mapview::mapview(mean_station_locs)
mapview::mapview(station_locs %>% filter(station_id == "LEKT_Station73"))


#calculate maximum distances between stations with the same ids between years
occ_dat %>% group_by(station_id)

station_ids <- occ_dat %>% distinct(station_id) %>% pull()

dists <- vector()
for(i in 1:length(station_ids)){
  rows <- occ_dat %>% 
    filter(station_id == station_ids[i]) %>% 
    st_as_sf(coords=c("x", "y"), crs = 5070) 
  
  dist_pairs <- vector()
  for(j in 1:(nrow(rows)-1)){
    row1 <- rows[j,]
    row2 <- rows[j+1,]
    
    dist_pairs[j] <- st_distance(row1, row2)
  }
  
  first_last <- st_distance(rows[1,], rows[nrow(rows),])
  
  dists[i] <- max(c(dist_pairs, first_last))
  
}

max_dists <- tibble(station_id = station_ids,
                    max_distance = dists)

max_dists %>% filter(max_distance > 30) %>% 
  filter(str_detect(station_id, "LEKT"))
