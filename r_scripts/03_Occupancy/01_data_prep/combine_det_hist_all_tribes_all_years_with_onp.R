#### Combine camera detection histories all tribes all years####

# Author: Read Barbee

# Date:2023-09-11 

#Last updated: 2023-10-06

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(janitor)
library(sf)

#########################################################################
##
## 1. Data import and initial formatting
##
##########################################################################

#Olympic National Park 2013-2016

#import ONP 2013-2016 detection history
onp_fisher <- read_csv("data/Camera_Data/Olympic_National_Park/onp_fisher_2013_2016_det_hist.csv") 

#import onp survey dates for each year for later formatting
onp_survey_periods <- read_csv("data/Camera_Data/Olympic_National_Park/onp_survey_periods_2013_2016.csv")

#get wgs84 coords for onp data
onp_wgs84 <- onp_fisher %>% 
  st_as_sf(coords=c("utm_e", "utm_n"), crs = 26910) %>% 
  st_transform(crs = 4326) %>% 
  st_coordinates()

#format columns and add wgs84 coordinates
onp_fisher <- onp_fisher %>% 
  unite("station", c(hex_id, station_num)) %>% 
  unite( "station_id", c(grid_id, station_id), sep = "_", remove = FALSE) %>% 
  mutate(lon = onp_wgs84[,1],
         lat = onp_wgs84[,2], .before = utm_e) %>% 
  select(station_id, grid_id, station, year, lon, lat, cell_id, bait_days_good, snare_days_good, effort_correction_factor, everything(), -utm_e, -utm_n)

#split onp data into separate frames by survey year
onp_split <- split(onp_fisher, as.factor(onp_fisher$year))
names(onp_split) <- c("onp_2013", "onp_2014", "onp_2015", "onp_2016")


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

#########################################################################
##
## 2. Batch format OCP data to match ONP format
##
##########################################################################

#combine all OCP data frames into a list
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

names(data_list) <- c(list("lekt_2019", 
                          "lekt_2020", 
                          "pnptc_2020", 
                          "lekt_2021", 
                          "pnptc_2021", 
                          "makah_2021", 
                          "skok_2021", 
                          "lekt_2022",
                          "pnptc_2022",
                          "makah_2022",
                          "quin_res_2022",
                          "quin_wyn_2022"))

#add columns to match ocp data list to onp format
data_list2 <- data_list %>% map(function(dat){dat %>% mutate(
  bait_days_good = 0,
  snare_days_good = 0,
  effort_correction_factor = 0,
  .after = cell_id)
  })


#########################################################################
##
## 3. Trim ONP frames for each year to survey period
##
##########################################################################

#function to remove date columns outside of a particular year for each onp list element - Working
trim_onp_dates <- function(dat){
  
  curr_year <- as.character(dat$year[1])
  
  survey_p <- onp_survey_periods %>% filter(year == curr_year)
  
  lower <- survey_p$min_date
  upper <- survey_p$max_date
  
  trimmed_dat <- dat %>% 
    pivot_longer(cols = -c(station_id:effort_correction_factor), names_to = "date", values_to = "value") %>% 
    mutate(date = ymd(date)) %>% 
    filter(date >= lower) %>% 
    filter(date <= upper) %>% 
    mutate(date = as.character(date)) %>% 
    pivot_wider(names_from = date, values_from = value)
  
  return(trimmed_dat)
}

#apply trimming function
onp_trimmed <- onp_split %>% map(trim_onp_dates)

#add onp survey frames to ocp list
data_list_full <- c(data_list2, onp_trimmed)


#########################################################################
##
## 4. Batch format all data: combine frames/calculate start and end dates for each year
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
              names_sort = TRUE,
              names_prefix = "d_")


#########################################################################
##
## 5. Extract covariate values for each site
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


#write_csv(occ_dat, "data/Camera_Data/master/ocp_onp_occ_dat_10-06-23.csv")



#########################################################################
##
## 6. Test: make sure final output can be converted to unmarked frame
##
##########################################################################

#tried to pivot detection covariates wider, but they need to be their own daily matrices :(
# occ_dat_v2 <- occ_dat %>% pivot_longer(cols = -c(station_id_year:effort_correction_factor), names_to = "date", values_to = "det") %>% 
#   pivot_wider(names_from = date, values_from = c(det, bait_days_good:effort_correction_factor))

# #Create stacked unmarked frame
# library(unmarked)
# 
# umf <- unmarkedFrameOccu(y = occ_dat[,12:398],
#                   siteCovs = occ_dat %>% select(cell_id, tree_cover_hansen:dist_water, bait_days_good:effort_correction_factor))
# 
# 
# library(ubms)
# 
# fit_stack <- stan_occu(~1~scale(tree_cover_hansen) + (1|station_id),
#                        data=umf, chains=3, iter=100)
# fit_stack


################################ Graveyard #################################


# trim_onp_dates <- function(dat){
#   
#   year <- dat$year[1]
#   
#   lower <- ymd(paste0(year, "-01-01")) - days(1)
#   upper <- ymd(paste0(year, "-12-31")) + days(1)
#   
#   trimmed_dat <- dat %>% 
#     pivot_longer(cols = -c(station_id:effort_correction_factor), names_to = "date", values_to = "value") %>% 
#     mutate(date = ymd(date)) %>% 
#     filter(date > lower) %>% 
#     filter(date < upper) %>% 
#     mutate(date = as.character(date)) %>% 
#     pivot_wider(names_from = date, values_from = value)
#   
#   return(trimmed_dat)
# }

#test <- onp_fisher %>% 
# filter(station_id =="ONP_21_1_2016") %>% 
#   pivot_longer(cols = -c(station_id:effort_correction_factor), names_to = "date", values_to = "value") %>%
#   mutate(date = ymd(date),
#          obs = yday(date)) %>% 
#   arrange(date) %>% 
#   select(-date) %>% 
#   pivot_wider(names_from = obs, 
#               values_from = value,
#               names_sort = TRUE,
#               names_prefix = "d_") %>% View()

#              values_fn = ~sum(.x, na.rm = TRUE)

