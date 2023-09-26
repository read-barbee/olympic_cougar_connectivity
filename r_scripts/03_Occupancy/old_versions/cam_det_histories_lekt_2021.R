#### Camera Detection Histories LEKT 2021 ####

# Author: Read Barbee

# Date:2023-08-11 

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

merge_rows <- function(row1, row2){
  if(row1$latitude==row2$latitude & row1$longitude==row2$longitude){
    row_new <- row1[1:10]
    row_new$cam_num <- paste0(row1$cam_num, "_", row2$cam_num)
    row_new$set_up_date <- min(c(row1$set_up_date, row2$set_up_date))
    row_new$pull_date <- max(c(row1$set_up_date, row2$set_up_date))
    hist1 <- row1[11:length(row1)]
    hist2 <- row2[11:length(row2)]
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
    hist_merge <- tibble(date = names(row1[11:length(row1)]), value = hist_merge) %>% pivot_wider(names_from = date, values_from = value)
    row_new <- bind_cols(row_new, hist_merge) 
  }
  
  else if (row1$latitude!=row2$latitude | row1$longitude!=row2$longitude){
    row1$station <- paste0(row1$station, "_a")
    row2$station <- paste0(row2$station, "_b")
    
    row_new <- bind_rows(row1, row2)
  }
  
  return(row_new)
}
station_level <- function(rows){
  if(nrow(rows)==1){
    new_row <- rows[1,]
  }
  else if(nrow(rows)==2){
    new_row <- merge_rows(rows[1,], rows[2,])
  }
  else if(nrow(rows)==3){
    out <- merge_rows(rows[1,], rows[2,])
    new_row <- merge_rows(out, rows[3,])
  }
  else if(nrow(rows)==4){
    out <- merge_rows(rows[1,], rows[2,])
    out2 <- merge_rows(out, rows[3,])
    new_row <- merge_rows(out2, rows[4,])
  }
  return(new_row) 
}


#########################################################################
##
## 1. Import and format activity and deteciton data
##
##########################################################################
#load("data/Camera_Data/2020/LEKT_2020/dtbs.2020.Rdata")

#camera activity-- aggregating to camera level isn't necessary because no stations were moved or had multiple cameras
lekt_2021_act <-read_csv("data/Camera_Data/2021/LEKT_2021/2021.Camera.Activity.Sheet.csv") %>%
  select(-c(study, elevation)) %>% 
  rename(longitude = x,
         latitude = y)

#species detections
lekt_2021_det <- read_csv("data/Camera_Data/2021/LEKT_2021/S3071_20210419_20211209_2022.09.12_20.04_metadata_tbl.csv") %>% 
  clean_names() %>% 
  select(station, cameras, x, y, date, species) %>% 
  rename(longitude = x,
         latitude = y) %>% 
  mutate(date = ymd(date)) %>% 
  filter(species=="Puma")

#########################################################################
##
## 2. Aggregate activity histories from the camera level to the station level
##
##########################################################################

#get list of distinct station names to iterate over
station_names <- lekt_2021_act %>% distinct(station) %>% pull()

#Execute aggregation function for each station
stations <- list()
for(i in 1:length(station_names)){
  rows <- lekt_2021_act %>% filter(station == station_names[i]) %>% arrange(cameras)
  stations[[i]] <- station_level(rows)
  
  print(i)
}

#bind list elements from for loop into single dataframe
stations <- bind_rows(stations)

stations2 <- stations %>% pivot_longer(cols =c(-c(station:cameras)),  names_to= "date", values_to = "value") %>% mutate(date = mdy(date))


#########################################################################
##
## 3. Join cougar detection counts to activity history by station name and date
##
##########################################################################

dets <- lekt_2021_det %>% group_by(station, date) %>% count() %>% rename(value = n) %>% 
  mutate(station = case_when(station=="Station96" ~ "Station96_b",
                             .default = station))



det_hist_full <- stations2 %>% 
  mutate(effort = value, 
         value = case_when(value==0 ~ NA,
                           .default = value)) %>% 
  left_join(dets, join_by(station, date)) %>% 
  mutate(value.x = case_when(is.na(value.y)==F ~ value.y,
                             value.x > 0 & is.na(value.y)==T ~ 0)) %>% 
  mutate(cougar_detections_binary = case_when(value.x >=1 ~ 1,
                                              .default = value.x)) %>% 
  select(-value.y) %>%  
  rename(cougar_detections_count = value.x) %>% 
  select(station:latitude, camera_id:date, cougar_detections_count, cougar_detections_binary, effort)


#########################################################################
##
## 3. Create binary and count detection histories
##
##########################################################################

det_hist_binary <- det_hist_full %>% 
  select(-c(cougar_detections_count, effort)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_binary) %>% 
  mutate(grid_id = "LEKT",
         year = "2021", .before=station)


det_hist_counts <- det_hist_full %>% 
  select(-c(cougar_detections_binary, effort)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_count) %>% 
  mutate(grid_id = "LEKT",
         year = "2021", .before=station)

effort_matrix <- det_hist_full %>% 
  select(-c(cougar_detections_binary, cougar_detections_count)) %>% 
  pivot_wider(names_from = date, values_from = effort) %>% 
  mutate(grid_id = "LEKT",
         year = "2021", .before=station)



# write_csv(det_hist_binary, "data/Camera_Data/2021/LEKT_2021/lekt_2021_det_hist.csv")
# write_csv(det_hist_counts, "data/Camera_Data/2021/LEKT_2021/lekt_2021_det_hist_counts.csv")
# write_csv(det_hist_full, "data/Camera_Data/2021/LEKT_2021/lekt_2021_det_hist_long.csv")
