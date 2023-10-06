#### Combine camera detection histories all tribes all years####

# Author: Read Barbee

# Date:2023-09-11 

#Last updated: 2023-10-05

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(janitor)
library(sf)

#########################################################################
##
## 1. Import and format detection data from each partner and year
##
##########################################################################

#ONP 2013-2016
onp_fisher <- read_csv("data/Camera_Data/Olympic_National_Park/onp_fisher_2013_2016_det_hist.csv") 

#get wgs84 coords for park data
onp_wgs84 <- onp_fisher %>% 
  st_as_sf(coords=c("utm_e", "utm_n"), crs = 26910) %>% 
  st_transform(crs = 4326) %>% 
  st_coordinates()
  
onp_fisher <- onp_fisher %>% 
  unite("station", c(hex_id, station_num)) %>% 
  unite( "station_id", c(grid_id, station_id), sep = "_", remove = FALSE) %>% 
  mutate(lon = onp_wgs84[,1],
         lat = onp_wgs84[,2], .before = utm_e) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, bait_days_good, snare_days_good, effort_correction_factor, everything(), -utm_e, -utm_n)


#2019

lekt_2019 <- read_csv("data/Camera_Data/2019/LEKT_2019/lekt_2019_det_hist.csv") %>% 
  select(-cameras) %>% 
  rename(station = station_id,
         lat = latitude,
         lon = longitude) %>% 
  relocate(station, .before = year) %>%
  relocate(lon, .before = lat) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE)


#2020

lekt_2020 <- read_csv("data/Camera_Data/2020/LEKT_2020/lekt_2020_det_hist.csv") %>% 
  select(-c(camera_id:pull_date)) %>% 
  rename(station = station_id,
         lat = latitude,
         lon = longitude) %>% 
  relocate(station, .before = year) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE)


pnptc_2020 <- read_csv("data/Camera_Data/2020/PNPTC_2020/pnptc_2020_det_hist.csv") %>% 
  rename(station = station_id) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())
  


#2021

lekt_2021 <- read_csv("data/Camera_Data/2021/LEKT_2021/lekt_2021_det_hist.csv") %>% 
  select(-c(camera_id, cameras)) %>% 
  rename(station = station_id,
         lat = latitude,
         lon = longitude) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())

pnptc_2021 <- read_csv("data/Camera_Data/2021/PNPTC_2021/pnptc_2021_det_hist.csv") %>% 
  rename(station = station_id,
         lat = latitude,
         lon = longitude) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())

makah_2021 <- read_csv("data/Camera_Data/2021/Makah_2021/makh_2021_det_hist.csv") %>% 
  rename(station = station_id) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())

skok_2021 <- read_csv("data/Camera_Data/2021/Skok_2021/skok_2021_det_hist.csv") %>% 
  select(-c(set_date, pull_date)) %>% 
  rename(station = station_id) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())


#2022

lekt_2022 <- read_csv("data/Camera_Data/2022/LEKT_2022/lekt_2022_det_hist.csv") %>% 
  select(-c(camera_id, cameras)) %>% 
  rename(station = station_id,
         lat = latitude,
         lon = longitude) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())


pnptc_2022 <- read_csv("data/Camera_Data/2022/PNPTC_2022/pnptc_2022_det_hist.csv") %>% 
  rename(station = station_id,
         lat = latitude,
         lon = longitude) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())


makah_2022 <- read_csv("data/Camera_Data/2022/Makah_2022/makh_2022_det_hist.csv") %>% 
  rename(station = station_id) %>% 
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything(), -camera_id)


quin_res_2022 <- read_csv("data/Camera_Data/2022/Quin_Res_2022/quin_res_2022_det_hist.csv") %>% 
  rename(station = station_id,
         lat = latitude,
         lon = longitude) %>%
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(-c(camera_id, cameras)) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())

quin_wyn_2022 <- read_csv("data/Camera_Data/2022/Quin_Wyno_2022/quin_wyno_2022_det_hist.csv") %>% 
  rename(station = station_id,
         lat = latitude,
         lon = longitude) %>%
  unite( "station_id", c(grid_id, station, year), sep = "_", remove = FALSE) %>% 
  select(-c(camera_id, cameras)) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, everything())

#combine all data frames into a list
data_list <- list(lekt_2019, 
               lekt_2020, 
               pnptc_2020, 
               lekt_2021, 
               pnptc_2021, 
               makah_2021, 
               skok_2021, 
               lekt_2022,
               pnptc_2022,
               makah_2022,
               quin_res_2022,
               quin_wyn_2022)

#add columns to match ocp data list to onp format
data_list2 <- data_list %>% map(function(dat){dat %>% mutate(
  bait_days_good = 0,
  snare_days_good = 0,
  effort_correction_factor = 0,
  .after = cell_id)
  })

#add onp data to the list

data_list_full <- c(data_list2, list(onp_fisher))

#########################################################################
##
## 2. Combine frames and calculate start and end dates for each year
##
##########################################################################

#pivot all data frames in list to long form
dat_pivot <- map(data_list_full, function(dat) dat %>% 
                   pivot_longer(cols = -c(station_id:effort_correction_factor), names_to = "date", values_to = "value") %>% 
                   mutate(date = ymd(date),
                          obs = yday(date), #convert to julian_day
                          #obs = paste0(month(date), "_" ,day(date)),
                          .before = value))

dat_pivot <- bind_rows(dat_pivot)

#summarize start and end dates for each year across all grids
survey_periods <- dat_pivot %>% 
  bind_rows() %>% 
  group_by(year) %>% 
  summarize(min_date = min(date), 
            max_date = max(date))


#combine rows and pivot wider
dat_wide <- dat_pivot %>%
  arrange(date) %>% 
  select(-date) %>% 
  pivot_wider(names_from = obs, 
              values_from = value,
              #names_sort = TRUE,
              names_prefix = "d_") %>% View()


#########################################################################
##
## 3. Extract covariate values for each site
##
##########################################################################

# Load the terra package
library(terra)
library(sf)

# Load your raster data
cov_stack <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_asp.tif")


# Load your sf object with points for station locations
sf_points <- dat_wide %>% st_as_sf(coords=c("lon", "lat"), crs = 4326, remove=FALSE) %>% 
  st_transform(crs = 5070)

# Extract cell numbers for each station location
covs <-terra::extract(cov_stack, sf_points)

#combine covariates with detection histories for each site
occ_dat <- dat_wide %>% 
  bind_cols(covs) %>%  
  select(-ID) %>% 
  rename(station_id_year = station_id) %>% 
  unite("station_id", c(grid_id, station), remove = FALSE) %>% 
  mutate(station_id = as.factor(station_id),
         year = as.factor(year)) %>% 
  relocate(c(northing, easting), .after = aspect)


#write_csv(occ_dat, "data/Camera_Data/master/ocp_occ_dat_10-02-23.csv")



# #Create stacked unmarked frame
# library(unmarked)
# 
# umf <- unmarkedFrameOccu(y = occ_dat[,9:351],
#                   siteCovs = occ_dat %>% select(station_id, tree_cover_hansen:dist_water))
# 
# 
# library(ubms)
# 
# fit_stack <- stan_occu(~1~scale(tree_cover_hansen) + (1|station_id), 
#                        data=umf, chains=3, iter=100)
# fit_stack

