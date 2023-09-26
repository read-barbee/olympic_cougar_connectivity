#### Camera Detection Histories 2020####

# Author: Read Barbee

# Date:2023-08-07 

# Purpose:

#multi-camera station situations
#1. camera replaced. merge rows
#2. camera moved. treat as two separate stations
#3. two cameras deployed simultaneously. Merge rows. Value= sum of effort

##Process
#1. consolidate activity history to station level
#2. pivot long for vector replacement
#3. replace 0s in activity sheet with NAs and 1s with 0s
#4. insert 1s by station and date for cougar detections
#5. Pivot wide


################################ Libraries #################################
library(tidyverse)
library(janitor)


#########################################################################
##
## 1. Import and format activity and deteciton data
##
##########################################################################
load("data/Camera_Data/2020/LEKT_2020/dtbs.2020.Rdata")

#camera activity-- aggregating to camera level isn't necessary because no stations were moved or had multiple cameras
lekt_2020_act <-camact %>%
  clean_names() %>% 
  select(station:camera_id, date, active) %>% 
  mutate(date = ymd(as.character(date))) 

#species detections
lekt_2020_det <- read_csv("data/Camera_Data/2020/LEKT_2020/S3071_20200513_20201210_2022.12.12_22.30_metadata_tbl.V2.csv") %>% 
  clean_names() %>% 
  mutate(date = date(mdy_hm(date_time_original))) %>% 
  select(station, cameras, x, y, date, species) %>% 
  filter(species=="Puma")


#########################################################################
##
## 2. Join cougar detection counts to activity history by station name and date
##
##########################################################################

dets <- lekt_2020_det %>% group_by(station, date) %>% count() %>% rename(value = n)

det_hist_full <- lekt_2020_act %>% 
  mutate(active = case_when(active==0 ~ NA,
                           .default = active)) %>% 
  left_join(dets, join_by(station, date)) %>% 
  mutate(cougar_detections_count = case_when(is.na(value)==F ~ value,
                             active > 0 & is.na(value)==T ~ 0)) %>% 
  mutate(cougar_detections_binary = case_when(cougar_detections_count >=1 ~ 1,
                                              .default = cougar_detections_count)) %>% 
  select(-value) %>% 
  rename(effort = active,
         latitude = x,
         longitude = y) %>% 
  select(station, camera_id, latitude:date, cougar_detections_count, cougar_detections_binary, effort) %>% 
  arrange(station, date)

#########################################################################
##
## 3. Create binary and count detection histories
##
##########################################################################

det_hist_binary <- det_hist_full %>% 
  select(-c(cougar_detections_count, effort)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_binary) %>% 
  mutate(grid_id = "LEKT",
         year = "2020", .before=station)


det_hist_counts <- det_hist_full %>% 
  select(-c(cougar_detections_binary, effort)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_count) %>% 
  mutate(grid_id = "LEKT",
         year = "2020", .before=station)

effort_matrix <- det_hist_full %>% 
  select(-c(cougar_detections_binary, cougar_detections_count)) %>% 
  pivot_wider(names_from = date, values_from = effort) %>% 
  mutate(grid_id = "LEKT",
         year = "2020", .before=station)



# write_csv(det_hist_binary, "data/Camera_Data/2020/LEKT_2020/lekt_2020_det_hist.csv")
# write_csv(det_hist_counts, "data/Camera_Data/2020/LEKT_2020/lekt_2020_det_hist_counts.csv")
# write_csv(det_hist_full, "data/Camera_Data/2020/LEKT_2020/lekt_2020_det_hist_long.csv")








