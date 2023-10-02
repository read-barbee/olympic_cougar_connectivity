#### Camera Detection Histories Quinault Wynoochee 2022 ####

# Author: Read Barbee

# Date:2023-10-02

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(janitor)

################################ Helper functions #################################

#merge rows function:
#if the latitude and longitude of row 1 = row 2:
#1. combine the camera numbers
#2. set the set up and pull dates to the minimum and maximum dates of the two rows
#3. sum the detection histories

#if the latitude and longitude of row 1 != row 2:
# keep rows separate and rename to station_a and station_b

#these two functions merge activity histories for cameras within the same 30x30mm raster pixel
merge_rows <- function(row1, row2){
  if(row1$cell_id==row2$cell_id){
    row_new <- row1[1:6]
    row_new$CameraID <- paste0(row1$CameraID, "_", row2$CameraID)
    hist1 <- row1[7:length(row1)]
    hist2 <- row2[7:length(row2)]
    hist_merge <- vector()
    for(i in 1:length(hist1)){
      if(is.na(hist1[i])==T & is.na(hist2[i])==F){
        hist_merge[i] <- hist2[i]
      } else if(is.na(hist1[i])==F & is.na(hist2[i])==T){
        hist_merge[i] <- hist1[i]
      } else if(is.na(hist1[i])==F & is.na(hist2[i])==F){
        hist_merge[i] <- (hist1[i] + hist2[i])
      }
      else if(is.na(hist1[i])==T & is.na(hist2[i])==T){
        hist_merge[i] <- NA
      }
    }
    hist_merge <- unlist(hist_merge)
    hist_merge <- tibble(date = names(row1[7:length(row1)]), value = hist_merge) %>% pivot_wider(names_from = date, values_from = value)
    row_new <- bind_cols(row_new, hist_merge) 
  }
  
  else if (row1$cell_id!=row2$cell_id){
    row1$station_id <- paste0(row1$station_id, "_a")
    row2$station_id <- paste0(row2$station_id, "_b")
    
    row_new <- bind_rows(row1, row2)
  }
  
  return(row_new)
}


station_level <- function(rows) {
  if (nrow(rows) == 1) {
    return(rows[1, ])
  }
  else{
  new_row <- rows[1, ]
  for (i in 2:nrow(rows)) {
    new_row <- merge_rows(new_row, rows[i, ])
  }
  
  return(new_row)
}}


#Function to merge images containing cougars that are fewer than 15 min apart
merge_events <- function(df) {
  df %>%
    arrange(station_id,timestamp) %>%
    group_by(station_id) %>%
    mutate(time_diff = c(0, diff(timestamp))) %>%
    mutate(grp = cumsum(time_diff >= 900)) %>%
    group_by(station_id,grp) %>%
    summarize(
      min_time = min(timestamp),
      max_time = max(timestamp),
      image_count = n()
    ) %>%
    mutate(date = date(min_time), .after=grp) %>% 
    ungroup() 
}



#########################################################################
##
## 1. Import and format activity and deteciton data
##
##########################################################################
#load("data/Camera_Data/2020/LEKT_2020/dtbs.2020.Rdata")

#camera activity-- aggregating to camera level isn't necessary because no stations were moved or had multiple cameras
cam_act <-read_csv("data/Camera_Data/2022/Quin_Wyno_2022/OlympicWyno_2022_ActivitySheet_FINAL.csv") %>%
  select(-c(Study, `General Location`, `Viewshed Area`, Elevation)) %>% 
  rename(longitude = Longitude,
         latitude = Latitude,
         station_id=Station,
         cameras = Cameras)

#species detections
cougar_det <- read_csv("data/Camera_Data/2022/Quin_Wyno_2022/OlympicWyno_2022_SpeciesDetections_FINAL.csv", col_types = cols(datetime = col_character())) %>% 
  clean_names() %>% 
  rename(timestamp = datetime,
         station_id = station) %>% 
  #filter(is.na(station_id) == F ) %>% 
  mutate(date_time = timestamp) %>% 
  separate_wider_delim(cols = date_time, delim = " ", names = c("date", "time")) %>% 
  mutate(timestamp = ymd_hms(timestamp),
         date = ymd(date),
         time = hms(time)) %>%
  #select(station_id, camera_id, date, time, timestamp, species) %>% 
  filter(species=="Cougar")

#########################################################################
##
## 2. Extract raster cell number for each camera station
##
##########################################################################
# Load the terra package
library(terra)
library(sf)

# Load your raster data
raster <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

temp <- raster[[1]]

# Load your sf object with points for station locations
sf_points <- cam_act %>% st_as_sf(coords=c("longitude", "latitude"), crs = 4326, remove=FALSE) %>% 
  st_transform(crs = 5070)

# Extract cell numbers for each station location
xy <- sf_points %>% st_coordinates()


#append raster cell id to each station location
cam_act <- cam_act %>% mutate(cell_id = cellFromXY(temp, xy), .after = latitude)

get_dupes(cam_act, cell_id)
get_dupes(cam_act, station_id)

########################################################################
##
## 3. Combine camera activity histories by raster cell
##
##########################################################################

#get list of distinct station names to iterate over
station_names <- cam_act %>% distinct(station_id) %>% pull()

#Execute aggregation function for each station
stations <- list()
for(i in 1:length(station_names)){
  rows <- cam_act %>% filter(station_id == station_names[i]) %>% arrange(CameraID)
  stations[[i]] <- station_level(rows)
  
  print(i)
}

#bind list elements from for loop into single dataframe
stations <- bind_rows(stations)

#check for duplicates or errors in aggregation
get_dupes(stations, station_id)

#pivot longer and clean up
stations_long <- stations %>%
  pivot_longer(cols=c(-c(station_id:cameras)), names_to= "date", values_to = "cam_status") %>%
  clean_names() %>%
  mutate(date = mdy(date))


#########################################################################
##
## 4. Generate count of independent cougar detections for each day at each camera station
##
##########################################################################

cougar_counts <- cougar_det %>% 
  merge_events() %>% 
  group_by(station_id, date) %>% 
  count() %>% 
  ungroup()


#########################################################################
##
## 5. Join cougar detection counts to activity history by station name and date
##
##########################################################################

#make long form detection history
#cam_status: NA = not deployed; 0 = deployed and inactive, 1 = deployed and active
#effort: 0 = undeployed/inactive, 1 = deployed and active, 2 = multiple cameras deployed and active
det_hist_full <- stations_long %>% 
  mutate(effort = case_when(is.na(cam_status)==T ~ 0,
                            .default = cam_status)) %>% 
  left_join(cougar_counts, join_by(station_id, date)) %>% 
  mutate(cougar_detections_count = case_when(is.na(n)==F ~ n,
                                             effort > 0 & is.na(n)==T ~ 0,
                                             is.na(cam_status)==T ~ NA)) %>%
  mutate(cougar_detections_binary = case_when(cougar_detections_count>=1 ~ 1,
                                              .default = cougar_detections_count)) %>% 
  select(-n) %>%  
  relocate(effort, .after = cam_status)



#########################################################################
##
## 6. Create binary and count detection histories
##
##########################################################################

det_hist_binary <- det_hist_full %>% 
  select(-c(cougar_detections_count, effort, cam_status)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_binary) %>% 
  mutate(grid_id = "QUIN_WYNO",
         year = "2022", .before=station_id)


det_hist_counts <- det_hist_full %>% 
  select(-c(cougar_detections_binary, effort, cam_status)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_count) %>% 
  mutate(grid_id = "QUIN_WYNO",
         year = "2022", .before=station_id)

effort_matrix <- det_hist_full %>% 
  select(-c(cougar_detections_binary, cougar_detections_count, cam_status)) %>% 
  pivot_wider(names_from = date, values_from = effort) %>% 
  mutate(grid_id = "QUIN_WYNO",
         year = "2022", .before=station_id)




# write_csv(det_hist_binary, "data/Camera_Data/2022/Quin_Wyno_2022/quin_wyno_2022_det_hist.csv")
# write_csv(det_hist_counts, "data/Camera_Data/2022/Quin_Wyno_2022/quin_wyno_2022_det_hist_counts.csv")
# write_csv(det_hist_full, "data/Camera_Data/2022/Quin_Wyno_2022/quin_wyno_2022_det_hist_long.csv")
