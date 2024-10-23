#### Makah combine years all species ####

# Author: Read Barbee

# Date:2024-10-13 

#Procedure

#1. For each activity matrix from each year:

#a. Consolidate activity matrices for each year by distinct location. If rows     have the same station name but different locations, name them station#_a       and station#_b, etc.

#b. Rename and reorder columns to fit standard format. Make sure each has a     column identifying the unique location-year combination (i.e. LEKT01_2020).    Pivot into long format so there is a column for each acive data and a column    for the activity status on that date (0, 1, NA). 


#2. Bind rows of long-form activity matrices into single data frame

#3. Pivot wider to get a column for each date across all years

#4. Check the wide-form activity history for errors

#5. Format detection files for each year

#6. Merge detection files into single dataframe

#7. Rename detection deployments to correspond with new activity matrix deployments

#8. Make sure camtrapR can make a detection history given the new activity and detection dataframes

#9. Export long and wide form activity matrices and formatted detection file

################################ Helper functions #################################

## sum_or_na ##
#for grouped rows in wideform activity histories: take the sum of rows by column or return NA if all values are NA
sum_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA)  # Return NA if all values are NA
  } else {
    return(sum(x, na.rm = TRUE))  # Return the sum if at least one value is not NA
  }
}


## merge_by_location ##
#merges activity histories for cameras at the same station and location in the same deployment year:

#For each station_id in a given year:
#if the latitude and longitude of rows match:
#1. sum the activity histories using the sum_or_na function.
#2. set the set up and pull dates to the minimum and maximum dates of the      matching rows
#if the latitude and longitude of rows DON'T match:
#1. Keep the rows separate and assign the station names suffixes _a, _b, etc.


merge_by_location <- function(data){
  out <- data %>% 
    group_by(station_id, longitude, latitude) %>%
    arrange(set_date) %>% 
    summarize(set_date = min(set_date),
              pull_date = max(pull_date),
              viewshed = first(viewshed),
              rows_merged = n(),
              across(-c(!matches("\\d")), sum_or_na), .groups = "drop") %>%
    ungroup() %>% 
    group_by(station_id) %>%                  # Group by station_id
    arrange(set_date) %>%                      # Arrange by date
    mutate(same_location = n_distinct(longitude, latitude) == 1,   # Check if all latitudes are the same
           suffix = ifelse(!same_location,                    # Add suffix if latitudes differ
                           paste0("_", letters[1:n()]), ""),
           station_id := paste0(station_id, suffix)) %>%
    ungroup() %>%
    select(-suffix, -same_location) %>% 
    arrange(station_id)
  
  return(out)
}


## start_end_from_matrix ##
#calculate start and end dates from dates between NA values for each row
start_end_from_matrix <- function(data){
  out <- data %>%
    arrange(station_id, act_date) %>%
    group_by(station_id) %>%
    summarise(
      set_date = min(act_date[!is.na(status) & lag(is.na(status), default = TRUE)], na.rm = TRUE),
      pull_date = coalesce(
        max(act_date[is.na(status) & lag(!is.na(status), default = FALSE)], na.rm = TRUE),
        max(act_date, na.rm = TRUE)  # Return the last date if no transition
      )
    ) %>%
    ungroup()
  
  return(out)
}

## fill_missing_dates ##
#add columns of NAs for each day between the min set date and max pull date across all years (mode = "survey") or all dates for relevant years (mode = "all") that don't exist in the current wide-form activity matrix

fill_missing_dates <- function(data, mode){
  
  date_cols <- data %>% 
    select(-c(survey_id:dep_year)) %>% names()
  
  if(mode == "survey"){
    start_date <- min(ymd(date_cols))
    end_date <- max(ymd(date_cols))
    
  }else if (mode == "all" | is.null(mode)){
    start_year <- year(min(ymd(date_cols)))
    end_year <- year(max(ymd(date_cols)))
    
    start_date <- ymd(paste0(start_year, "/01/01"))
    end_date <- ymd(paste0(end_year, "/12/31"))
  }
  
  date_seq <- seq(start_date, end_date, by = "day") %>% as.character()
  
  for(i in date_seq){
    if(!(as.character(i) %in% date_cols))
      data <- data %>% 
        mutate(!!as.character(i) := NA)
  }
  
  out <- data %>%
    pivot_longer(-c(survey_id:dep_year), names_to = "name", values_to = "value") %>%
    mutate(name = ymd(name)) %>%
    arrange(name) %>%
    #mutate(name = as.character(name)) %>%
    pivot_wider(names_from = name, values_from = value)
  
  return(out)
}

## fill_jdays ##
#fill missing juilan days in each year for stacking the detection history
fill_jdays <- function(data){
  
  date_cols <- data %>%
    select(-deployment_id) %>% names()
  
  for(i in 1:365){
    if(!(as.character(i) %in% date_cols))
      data <- data %>%
        mutate(!!as.character(i) := NA)
  }
  
  out <- data %>%
    pivot_longer(-deployment_id, names_to = "name", values_to = "value") %>%
    mutate(name = as.numeric(name)) %>%
    arrange(name) %>%
    pivot_wider(names_from = name, values_from = value)
  
  return(out)
}


# Function to read and process each MAKAH CSV file
process_csv <- function(file_path) {
  # Extract the filename without extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Add a new column with the filename
  data_with_filename <- data %>%
    clean_names() %>% 
    select(species, capture_date_local, capture_time_local) %>% 
    rename(date = capture_date_local,
           time = capture_time_local) %>% 
    mutate(filename = file_name, .before=species,
           timestamp = as.character(ymd_hms(paste(date, time))),
           date = as.character(ymd(date)),
           time = as.character(time)) %>% 
    separate_wider_delim(cols = filename, delim = ".", names = c("grid", "dep_year", "camera", "check"))
  
  return(data_with_filename)
}

################################ Libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)

###############################################################################

#1. Import viewshed files

###############################################################################

#None

###############################################################################

#2. Format individual activity sheets

###############################################################################

################################ 2021 #################################
makah_act_2021 <- read_csv("data/Camera_Data/all_species/MAKAH/2021/ActivitySheet_FINAL_2.1.22_Number2_edited.csv")


#Create a key to translate between old and new deployment names
makah_key_2021 <- makah_act_2021 %>% 
  clean_names() %>%
  select(station, cameras, x, y) %>% 
  rename(station_id = station,
         camera_id = cameras,
         longitude = x,
         latitude = y) %>% 
  mutate(dep_year = "2021", .before = station_id)


#check the number of distinct locations (for number of output rows)
target_rows_2021 <- makah_act_2021 %>% 
  distinct(X, Y) %>% 
  nrow()


makah_act_2021_f <- makah_act_2021 %>% 
  clean_names() %>% 
  select(-c(study, elevation, camera_id)) %>% 
  rename(station_id = station,
         longitude = x,
         latitude = y,
         camera_id = cameras) %>% 
  mutate(station_id =  str_replace(station_id,  "(\\D)(\\d)\\b", "\\10\\2")) %>% 
  mutate(viewshed = NA,
         station_id = str_replace(station_id, coll("Station"), coll("MAKAH"))) %>% 
  relocate(viewshed, .after = camera_id) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(act_date = dmy(act_date)) %>% 
  relocate(camera_id, .after = station_id) 

start_end_dates_2021 <- makah_act_2021_f %>%
  start_end_from_matrix() 


makah_act_2021_wide <- makah_act_2021_f %>% 
  left_join(start_end_dates_2021, by = "station_id") %>%
  relocate(set_date, pull_date, .after = viewshed) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  merge_by_location()


#check if the correct number of rows were returned
nrow(makah_act_2021_wide) == target_rows_2021



makah_act_2021_long <- makah_act_2021_wide %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(dep_year = "2021", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>%
  mutate(survey_id = "MAKAH_2021", .before = station_id) %>% 
  mutate(act_date = as.character(act_date),
         set_date = as.character(set_date),
         pull_date = as.character(pull_date))

  


###############################################################################

#3. Combine activity sheets in long and wide formats

###############################################################################

makah_act_all <- makah_act_2021_long


loc_keys_all <- makah_key_2021


makah_act_all_long <- makah_act_all %>% 
  unite("deployment_id", survey_id, station_id,  sep="_", remove = F) %>% 
  relocate(survey_id, .before = deployment_id) %>% 
  select(-act_year)


#check for duplicates in deployment id and activity date
makah_act_all_long  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date)


#fill missing dates (periods of NA between surveys that weren't included in annual activity matrices)
makah_act_all_long_date_fill <- makah_act_all_long %>% 
  pivot_wider(names_from = act_date, values_from = status) %>%
  fill_missing_dates(mode = "survey") %>% 
  pivot_longer(-c(survey_id:dep_year), names_to = "act_date", values_to = "status")



#check for duplicates in deployment id and activity date
makah_act_all_long_date_fill  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date) #%>% distinct(deployment_id)

#lmit total number of active cameras per day at a location to 1 in case of overlap errors in matrix
makah_act_all_final_long <- makah_act_all_long_date_fill %>% 
  mutate(status = case_when(status > 1 ~ 1,
                            .default = status)) 

# #camera_level_test
makah_act_all_final_long  %>%
  #filter(dep_year == "2019") %>% 
  arrange(dep_year, set_date) %>% 
  select(deployment_id, act_date, status) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  column_to_rownames("deployment_id") %>%
  as.matrix() %>%
  camtrapR:::camopPlot(., lattice = TRUE)


#pivot from long form to wide form
makah_act_all_final_wide <- makah_act_all_final_long %>%
  pivot_wider(names_from = act_date, values_from = status)

#check for duplicates in deployment_id
makah_act_all_final_wide %>% get_dupes(deployment_id)



#########################################################################
##
## 1. Combine detections for all species across years
##
##########################################################################
# Set the path to the folder containing CSV files
folder_path <- "data/Camera_Data/all_species/MAKAH/2021/makah_checks_2021"

# Get a list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

# Read and process each CSV file
dets_2021 <- map(csv_files, process_csv)

# Combine all processed data into a single data frame
makah_dets_all <- bind_rows(dets_2021) %>% 
  select(-c(grid, check)) %>% 
  rename(station_id = camera) %>% 
  mutate(station_id = str_replace(station_id, coll("CAM"), coll("Station"))) %>% 
  mutate(station_id = str_replace(station_id, coll("Cam"), coll("Station"))) %>% 
  mutate(camera_id = case_when(str_detect(station_id, "_2") ~ "Camera2",
                               str_detect(station_id, "_3") ~ "Camera3",
                               str_detect(station_id, "_4") ~ "Camera4",
                               str_detect(station_id, "_5") ~ "Camera5",
                               str_detect(station_id, "_6") ~ "Camera6",
                               .default = "Camera1"), .after = station_id)
  


###############################################################################

#5. Change detection deployment names to match activity history

###############################################################################

#Make a deployment key to link old detection metadata to new deployment structure
dep_key <- makah_act_all_final_wide %>% 
  select(deployment_id, longitude, latitude, dep_year) %>% 
  left_join(loc_keys_all, by = c("dep_year","longitude", "latitude"))


#Rename the deployment and station_ids in the detection file 
makah_dets_all_final <- makah_dets_all %>%
  mutate(station_id = str_remove(station_id, "_.*")) %>% 
  left_join(dep_key, by = c("dep_year", "station_id", "camera_id")) %>% 
  select(-c(station_id, camera_id)) %>% 
  mutate(date = as.character(date),
         time = as.character(time),
         act_year = year(timestamp)) %>% 
  mutate(timestamp = paste0(date, " ", time)) %>% 
  mutate(timestamp = as.character(timestamp)) %>% 
  select(deployment_id, longitude, latitude, dep_year, act_year, species, timestamp, date, time) %>% 
  filter(act_year!= "2020")
  

makah_dets_all_final %>% 
  mutate(unique_id = 1:nrow(.)) %>% #filter(unique_id == 37010)
  mutate(timestamp = ymd_hms(timestamp)) %>% filter(is.na(timestamp)) 

#check for mismatches between deployment ids in detection and activty history frames
det_stations <- makah_dets_all_final %>% distinct(deployment_id) %>% pull(deployment_id)

act_stations <- makah_act_all_final_long %>% distinct(deployment_id) %>% pull(deployment_id)

setdiff(act_stations, det_stations) #stations in act not in det
setdiff(det_stations, act_stations) #stations in det not in act

###############################################################################

#6. Test the camtrapR detectionHistory function using the new frames

###############################################################################


makah_act_all_final_wide_mat <- makah_act_all_final_wide %>% 
  arrange(set_date) %>% 
  select(-c(survey_id, station_id:dep_year)) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

test <- detectionHistory(makah_dets_all_final,
                         "Deer",
                         makah_act_all_final_wide_mat,
                         output = "binary",
                         stationCol = "deployment_id",
                         speciesCol = "species",
                         recordDateTimeCol = "timestamp",
                         recordDateTimeFormat = "ymd HMS",
                         occasionLength = 1,
                         includeEffort = F,
                         day1 = "survey",
                         datesAsOccasionNames = F)

det_hist <- test$detection_history

colnames(det_hist) <- colnames(makah_act_all_final_wide_mat)

camtrapR:::camopPlot(det_hist, lattice = TRUE)

###############################################################################

#7. Try to stack the detection history-- (for when more years are added)

###############################################################################


# test2 <- det_hist %>%
#   as.data.frame() %>% 
#   rownames_to_column("deployment_id") %>%
#   pivot_longer(-deployment_id, names_to = "act_date", values_to = "detection") %>% 
#   mutate(act_date = ymd(act_date),
#          jday = yday(act_date)) %>% 
#   mutate(act_year = as.factor(year(act_date))) 
# 
# 
# test3 <- test2 %>% 
#   split(.$act_year) %>% 
#   map(., function(x){
#     x %>% select(-c(act_year, act_date)) %>% 
#       pivot_wider(names_from = jday, values_from = detection)
#   }) %>% 
#   map(., fill_jdays)
# 
# for(i in 1:length(test3)){
#   year_name <- names(test3)[i]
#   
#   test3[[i]] <- test3[[i]] %>% 
#     filter(str_detect(deployment_id, pattern = year_name))
# }
# 
# test4 <- test3 %>% 
#   bind_rows() #%>% 
# #distinct(deployment_id, .keep_all = T) 
# 
# 
# #filter(!all(is.na(c_across(-deployment_id))))
# 
# 
# date_seq <- seq(mdy("01/01/2021"), mdy("12/31/2021"), by = "day") %>% as.character()
# 
# plot_mat <- test4 %>% 
#   column_to_rownames("deployment_id") %>% 
#   as.matrix()
# 
# colnames(plot_mat) <- date_seq
# 
# camtrapR:::camopPlot(plot_mat, lattice = TRUE)

###############################################################################

#8. Export

###############################################################################

# write_csv(makah_act_all_final_long, "data/Camera_Data/all_species/makah_cam_act_2021_long.csv")
# 
# write_csv(makah_act_all_final_wide, "data/Camera_Data/all_species/makah_cam_act_2021_wide.csv")
# 
# write_csv(makah_dets_all_final, "data/Camera_Data/all_species/makah_detections_all_species_2021.csv")
# 
# write_csv(dep_key, "data/Camera_Data/all_species/makah_deployment_key_2021.csv")




