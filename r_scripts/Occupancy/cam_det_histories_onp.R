#### Olympic National Park Detection History Formatting ####

# Author: Read Barbee

# Date:2023-09-11 

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(janitor)

#########################################################################
##
## 1. Import and format detection data 
##
##########################################################################

dat <- read_csv("data/Camera_Data/Olympic_National_Park/ONP_fisher_grid_2013-2016_raw_formatted.csv") %>% 
  clean_names() %>% 
  select(-c(rep_id, study_incid, waypoint_type)) %>% 
  mutate(cougar = case_when(cougar >1 ~ 1,
                             .default = cougar))


dat_wide<- dat %>% 
  mutate(station = paste0(hex_id, "_", station_num), .after=station_id) %>% 
  #group_by(station_id) %>%
  pivot_wider(names_from = visit_num, values_from =visit_date:cougar) %>% 
  select(station_id, station, year, cougar_1:cougar_5, everything())



#########################################################################
##
## 2. Extract covariate values for each site
##
##########################################################################

# Load the terra package
library(terra)
library(sf)

# Load your raster data
cov_stack <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

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


# Load your sf object with points for station locations. native crs is NAD83 / UTM zone 10N
sf_points <- dat_wide %>% st_as_sf(coords=c("utm_e", "utm_n"), crs = 26910, remove=FALSE) %>% 
  st_transform(crs=5070)

# Extract cell numbers for each station location
covs <-terra::extract(cov_stack, sf_points)

occ_dat <- dat_wide %>% 
  bind_cols(covs) %>%  
  select(-ID) %>% 
  mutate(station_id = str_remove(station_id, "_r0")) %>% 
  mutate(station_id = str_remove(station_id, "_r1")) %>% 
  mutate(station_id = str_remove(station_id, "_r2")) %>% 
  mutate(across(visit_date_1:visit_date_5, mdy)) %>% 
  mutate(across(visit_date_1:visit_date_5, yday))


#Create stacked unmarked frame
library(unmarked)

umf <- unmarkedFrameOccu(y = occ_dat[,4:8],
                         siteCovs = occ_dat %>% select(station_id, tree_cover_hansen:dist_water),
                         obsCovs = list(date = occ_dat %>% select(visit_date_1:visit_date_5),
                                        int = occ_dat %>% select(interval_1:interval_5),
                                        effort = occ_dat %>% select(camera_days_good_1:camera_days_good_5),
                                        bait_tree = occ_dat %>% select(bait_days_good_1:bait_days_good_5),
                                        bait_snare = occ_dat %>% select(snare_days_good_1:snare_days_good_5)))






library(ubms)

fit_stack <- stan_occu(~ scale(date) + scale(int) + scale(effort) + scale(bait_tree) + scale(bait_snare)  ~scale(tree_cover_hansen) + (1|station_id), 
                       data=umf, chains=3, iter=100)
fit_stack






