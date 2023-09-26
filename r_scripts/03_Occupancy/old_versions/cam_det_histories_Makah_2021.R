#### Makah Detection History Formatting 2021 ####

# Author: Read Barbee

# Date:2023-08-31 

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
  row1$X <- round(row1$X, digits=5)
  row1$Y <- round(row1$Y, digits=5)
  row2$X <- round(row2$X, digits=5)
  row2$Y <- round(row2$Y, digits=5)
  if(row1$Y==row2$Y & row1$X==row2$X){
    row_new <- row1[1:7]
    row_new$CameraID <- paste0(row1$CameraID, "_", row2$CameraID)
    hist1 <- row1[8:length(row1)]
    hist2 <- row2[8:length(row2)]
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
    hist_merge <- tibble(date = names(row1[8:length(row1)]), value = hist_merge) %>% pivot_wider(names_from = date, values_from = value)
    row_new <- bind_cols(row_new, hist_merge) 
  }
  
  else if (row1$Y!=row2$Y | row1$X!=row2$X){
    row1$Station <- paste0(row1$Station, "_a")
    row2$Station <- paste0(row2$Station, "_b")
    
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
  else if(nrow(rows)==5){
    out <- merge_rows(rows[1,], rows[2,])
    out2 <- merge_rows(out, rows[3,])
    out3 <- merge_rows(out2, rows[4,])
    new_row <- merge_rows(out3, rows[5,])
  }
  else if(nrow(rows)==6){
    out <- merge_rows(rows[1,], rows[2,])
    out2 <- merge_rows(out, rows[3,])
    out3 <- merge_rows(out2, rows[4,])
    out4 <- merge_rows(out3, rows[5,])
    new_row <- merge_rows(out4, rows[6,])
  }
  return(new_row) 
}


#########################################################################
##
## 1. Create cougar detection file from individual camera check csv files
##
##########################################################################

# Set the path to the folder containing CSV files
folder_path <- "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Camera_Data/2021/Makah_2021/Makah_CameraGrid2021_ForRead/makah_checks_2021"

# Get a list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

# Create a function to read and process each CSV file
process_csv <- function(file_path) {
  # Extract the filename without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Add a new column with the filename
  data_with_filename <- data %>%
    clean_names() %>% 
    select(species, capture_date_local) %>% 
    mutate(filename = file_name, .before=species,
           capture_date_local = ymd(capture_date_local)) %>% 
    separate_wider_delim(cols = filename, delim = ".", names = c("tribe", "year", "camera", "check"))
  
  return(data_with_filename)
}

# Read and process each CSV file
all_data <- lapply(csv_files, process_csv)

# Combine all processed data into a single data frame
combined_data <- bind_rows(all_data)


#filter just for cougar detections and relevant columns
cougar_det <- combined_data %>% 
  filter(species=="Cougar") %>% 
  select(camera, species, capture_date_local) %>% 
  rename(camera_id = camera, 
         date = capture_date_local)

#create dataframe key linking stations to camera names
cam_to_stat <- cam_act %>% select(Station, CameraID) %>% clean_names()

#aggregate individual cougar detections by station
cougar_counts <- cougar_det %>% 
  group_by(camera_id, date) %>% 
  count() %>% 
  ungroup() %>% 
  left_join(cam_to_stat, by=join_by(camera_id)) %>% 
  select(station, date, n)


#########################################################################
##
## 2. Import and format camera activity file
##
##########################################################################


cam_act <- read_csv("data/Camera_Data/2021/Makah_2021/Makah_CameraGrid2021_ForRead/ActivitySheet_FINAL_2.1.22_Number2.csv")


#########################################################################
##
## 3. Aggregate activity histories from the camera level to the station level
##
##########################################################################

#get list of distinct station names to iterate over
station_names <- cam_act %>% distinct(Station) %>% pull()

#Execute aggregation function for each station
stations <- list()
for(i in 1:length(station_names)){
  rows <- cam_act %>% filter(Station == station_names[i]) %>% arrange(Cameras)
  stations[[i]] <- station_level(rows)
  
  print(i)
}

#bind list elements from for loop into single dataframe
stations <- bind_rows(stations)

#pivot longer and clean up
stations2 <- stations %>%
  pivot_longer(cols=c(-c(Study:Elevation)), names_to= "date", values_to = "value") %>%
  clean_names() %>%
  select(-c(study, elevation)) %>%
  mutate(date = dmy(date))


#########################################################################
##
## 4. Join cougar detection counts to activity history by station name and date
##
##########################################################################

det_hist_full <- stations2 %>% 
  mutate(effort = value, 
         value = case_when(value==0 ~ NA,
                           .default = value)) %>% 
  left_join(cougar_counts, join_by(station, date)) %>% 
  mutate(value = case_when(is.na(n)==F ~ n,
                             value > 0 & is.na(n)==T ~ 0)) %>%
  mutate(cougar_detections_binary = case_when(value>=1 ~ 1,
                                              .default = value)) %>% 
  select(-n) %>%  
  rename(cougar_detections_count = value) %>% 
  select(-cameras) %>% 
  relocate(effort, .after = cougar_detections_binary)


#########################################################################
##
## 5. Create binary and count detection histories
##
##########################################################################

det_hist_binary <- det_hist_full %>% 
  select(-c(cougar_detections_count, effort)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_binary) %>% 
  mutate(grid_id = "MAKH",
         year = "2021", .before=station)


det_hist_counts <- det_hist_full %>% 
  select(-c(cougar_detections_binary, effort)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_count) %>% 
  mutate(grid_id = "MAKH",
         year = "2021", .before=station)

effort_matrix <- det_hist_full %>% 
  select(-c(cougar_detections_binary, cougar_detections_count)) %>% 
  pivot_wider(names_from = date, values_from = effort) %>% 
  mutate(grid_id = "MAKH",
         year = "2021", .before=station)


# 
# write_csv(det_hist_binary, "data/Camera_Data/2021/Makah_2021/makh_2021_det_hist.csv")
# write_csv(det_hist_counts, "data/Camera_Data/2021/Makah_2021/makh_2021_det_hist_counts.csv")
# write_csv(det_hist_full, "data/Camera_Data/2021/Makah_2021/makh_2021_det_hist_long.csv")


