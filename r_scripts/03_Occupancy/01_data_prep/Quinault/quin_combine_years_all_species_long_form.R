#### Camera Detection Histories Quinault Wynoochee 2022 ####

# Author: Read Barbee

# Date:2023-10-02

# Purpose:
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



################################ Libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)


#########################################################################
##
## 1. Import and format activity and deteciton data
##
##########################################################################

################################ 2022 - Res #################################
quin_act_2022_res <- read_csv("data/Camera_Data/all_species/QUIN/2022/Res/OlympicQuin_2022_ActivitySheet_FINAL.csv")

#Create a key to translate between old and new deployment names
loc_key_2022_res <- quin_act_2022_res %>% 
  clean_names() %>% 
  select(station, camera_id, longitude, latitude) %>% 
  rename(station_id = station) %>% 
  # mutate(camera_id = case_when(str_detect(camera_id, coll("_2")) ~ "Camera2",
  #                              str_detect(camera_id, coll("_3")) ~ "Camera3",
  #                                         .default = "Camera1")) %>% 
  mutate(dep_year = "2022", .before = station_id)

#check the number of distinct locations (for number of output rows)
target_rows_2022_res <- quin_act_2022_res %>% 
  clean_names() %>% 
  distinct(longitude, latitude) %>% 
  nrow()


#format column names and merge rows by location
quin_act_2022_res_f <- quin_act_2022_res %>% 
  clean_names() %>% 
  select(-c(study, general_location, cameras, elevation)) %>% 
  rename(station_id = station,
         viewshed = viewshed_area) %>% 
  relocate(longitude, .before = latitude) %>% 
pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>%
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(act_date = mdy(act_date)) 


start_end_dates_2022_res <- quin_act_2022_res_f %>%
  start_end_from_matrix() 

quin_act_2022_res_wide <- quin_act_2022_res_f %>% 
  left_join(start_end_dates_2022_res, by = "station_id") %>%
  relocate(set_date, pull_date, .after = viewshed) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  merge_by_location()

#check if the correct number of rows were returned
nrow(quin_act_2022_res_wide) == target_rows_2022_res


quin_act_2022_long_res <- quin_act_2022_res_wide  %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(dep_year = "2022", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>%
  mutate(survey_id = "QIR_2022", .before = station_id) %>% 
  mutate(act_date = as.character(act_date),
         set_date = as.character(set_date),
         pull_date = as.character(pull_date))
  
  
################################ 2022 - Wyn #################################
quin_act_2022_wyn <- read_csv("data/Camera_Data/all_species/QUIN/2022/Wyno/OlympicWyno_2022_ActivitySheet_FINAL.csv")

#Create a key to translate between old and new deployment names
loc_key_2022_wyn <- quin_act_2022_wyn %>% 
  clean_names() %>% 
  select(station, camera_id, longitude, latitude) %>% 
  rename(station_id = station) %>%
  # mutate(camera_id = case_when(str_detect(camera_id, coll("_2")) ~ "Camera2",
  #                              str_detect(camera_id, coll("_3")) ~ "Camera3",
  #                              .default = "Camera1")) %>% 
  mutate(dep_year = "2022", .before = station_id)

#check the number of distinct locations (for number of output rows)
target_rows_2022_wyn <- quin_act_2022_wyn %>% 
  clean_names() %>% 
  distinct(longitude, latitude) %>% 
  nrow()


#format column names and merge rows by location
quin_act_2022_wyn_f <- quin_act_2022_wyn %>% 
  clean_names() %>% 
  select(-c(study, general_location, cameras, elevation)) %>% 
  rename(station_id = station,
         viewshed = viewshed_area) %>% 
  relocate(longitude, .before = latitude) %>% 
  rowwise() %>% 
  filter(!all(is.na(c_across(-c(station_id:viewshed))))) %>% #remove stations with all NAS in activity
  ungroup() %>%
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>%
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(act_date = mdy(act_date))


start_end_dates_2022_wyn <- quin_act_2022_wyn_f %>%
  start_end_from_matrix() 

quin_act_2022_wyn_wide <- quin_act_2022_wyn_f %>% 
  left_join(start_end_dates_2022_res, by = "station_id") %>%
  relocate(set_date, pull_date, .after = viewshed) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  merge_by_location()

#check if the correct number of rows were returned
nrow(quin_act_2022_wyn_wide) == target_rows_2022_wyn - 2


quin_act_2022_long_wyn <- quin_act_2022_wyn_wide  %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(dep_year = "2022", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>%
  mutate(survey_id = "WYN_2022", .before = station_id) %>% 
  mutate(act_date = as.character(act_date),
         set_date = as.character(set_date),
         pull_date = as.character(pull_date))



###############################################################################

#3. Combine activity sheets in long and wide formats

###############################################################################

quin_act_all <- bind_rows(quin_act_2022_long_res,
                          quin_act_2022_long_wyn)


loc_keys_all <- bind_rows(loc_key_2022_res,
                          loc_key_2022_wyn)


quin_act_all_long <- quin_act_all %>% 
  unite("deployment_id", survey_id, station_id,  sep="_", remove = F) %>% 
  relocate(survey_id, .before = deployment_id) %>% 
  select(-act_year)


#check for duplicates in deployment id and activity date
quin_act_all_long  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date)


#fill missing dates (periods of NA between surveys that weren't included in annual activity matrices)
quin_act_all_long_date_fill <- quin_act_all_long %>% 
  pivot_wider(names_from = act_date, values_from = status) %>%
  fill_missing_dates(mode = "survey") %>%
  pivot_longer(-c(survey_id:dep_year), names_to = "act_date", values_to = "status")



#check for duplicates in deployment id and activity date
quin_act_all_long_date_fill  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date) #%>% distinct(deployment_id)

#lmit total number of active cameras per day at a location to 1 in case of overlap errors in matrix
quin_act_all_final_long <- quin_act_all_long_date_fill %>% 
  mutate(status = case_when(status > 1 ~ 1,
                            .default = status))

# #camera_level_test
quin_act_all_final_long  %>%
  #filter(dep_year == "2019") %>% 
  arrange(dep_year, station_id, act_date) %>% 
  select(deployment_id, act_date, status) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  column_to_rownames("deployment_id") %>%
  as.matrix() %>%
  camtrapR:::camopPlot(., lattice = TRUE)


#pivot from long form to wide form
quin_act_all_final_wide <- quin_act_all_final_long %>%
  pivot_wider(names_from = act_date, values_from = status)

#check for duplicates in deployment_id
quin_act_all_final_wide %>% get_dupes(deployment_id)


#Make a deployment key to link old detection metadata to new deployment structure
dep_key <- quin_act_all_final_wide %>%
  select(deployment_id, longitude, latitude, dep_year) %>% 
  left_join(loc_keys_all, by = c("dep_year","longitude", "latitude")) 


#########################################################################
##
## 4. Combine detections for all species across years
##
##########################################################################

################################ 2022 res #################################

quin_det_2022_res <- read_csv("data/Camera_Data/all_species/QUIN/2022/Res/OlympicQuin_2022_SpeciesDetections_FINAL.csv") %>% 
  clean_names() %>% 
  rename(timestamp = datetime,
         station_id = station) %>% 
  mutate(timestamp = paste0(timestamp, ":00")) %>% 
  mutate(timestamp = as.character(mdy_hms(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2022", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  mutate(timestamp = case_when(is.na(time) ~ paste0(timestamp, " 00:00:00"),
                               .default = timestamp),
         time = case_when(is.na(time) ~ "00:00:00",
                               .default = time)) %>% 
  #mutate(timestamp = as.character(ymd_hms(timestamp))) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)

ymd_hms(quin_det_2022_res$timestamp)

# quin_det_2022_res %>% 
#   mutate(timestamp = ymd_hms(timestamp)) %>% 
#   filter(is.na(timestamp))

################################ 2022 wyno #################################
quin_det_2022_wyn <- read_csv("data/Camera_Data/all_species/QUIN/2022/Wyno/OlympicWyno_2022_SpeciesDetections_FINAL.csv") %>% 
  clean_names() %>% 
  rename(timestamp = datetime,
         station_id = station,
         camera_id = camera) %>% 
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2022", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)


ymd_hms(quin_det_2022_wyn$timestamp)

################################ Combine grids #################################
quin_dets_all <- bind_rows(quin_det_2022_res,
                           quin_det_2022_wyn) 


###############################################################################

#5. Change detection deployment names to match activity history

###############################################################################

dep_key2 <- dep_key %>% 
  mutate(camera_id = case_when(deployment_id == "WYN_2022_WYN55_a" ~ paste0(camera_id, "a"),
                               deployment_id == "WYN_2022_WYN55_b" ~ paste0(camera_id, "b"),
                               .default = camera_id),
         station_id = case_when(deployment_id == "WYN_2022_WYN55_a" ~ paste0(station_id, "a"),
                               deployment_id == "WYN_2022_WYN55_b" ~ paste0(station_id, "b"),
                               .default = station_id))

#Rename the deployment and station_ids in the detection file 
quin_dets_all_final <- quin_dets_all %>%
  mutate(camera_id = str_replace(camera_id, coll("Cam"), coll("CAM"))) %>% 
  left_join(dep_key2, by = c("dep_year", "station_id", "camera_id")) %>% 
  select(-c(station_id, camera_id)) %>% 
  select(deployment_id, longitude, latitude, everything())



#check for mismatches between deployment ids in detection and activty history frames
det_stations <- quin_dets_all_final %>% distinct(deployment_id) %>% pull(deployment_id)

act_stations <- quin_act_all_final_long %>% distinct(deployment_id) %>% pull(deployment_id)

setdiff(act_stations, det_stations) #stations in act not in det
setdiff(det_stations, act_stations) #stations in det not in act



###############################################################################

#6. Test the camtrapR detectionHistory function using the new frames

###############################################################################

#remove any NA values in the species column
quin_dets_all_final2 <- quin_dets_all_final %>% 
  filter(!is.na(species))

quin_act_all_final_wide_mat <- quin_act_all_final_wide %>% 
  select(-c(survey_id, station_id:dep_year)) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

test <- detectionHistory(quin_dets_all_final2,
                         species = "Cougar",
                         quin_act_all_final_wide_mat,
                         output = "binary",
                         stationCol = "deployment_id",
                         speciesCol = "species",
                         recordDateTimeCol = "timestamp",
                         recordDateTimeFormat = "ymd HMS",
                         occasionLength = 1,
                         includeEffort = F,
                         day1 = "survey",
                         datesAsOccasionNames = T
                         )

det_hist <- test$detection_history

#colnames(det_hist) <- colnames(quin_act_all_final_wide_mat)

camtrapR:::camopPlot(det_hist, lattice = TRUE)


###############################################################################

#7. Try to stack the detection history- for when there is more than one year

###############################################################################

# 
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
#distinct(deployment_id, .keep_all = T) 


#filter(!all(is.na(c_across(-deployment_id))))

# 
# date_seq <- seq(mdy("01/01/2020"), mdy("12/31/2020"), by = "day") %>% as.character()
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
# 
# write_csv(quin_act_all_final_long, "data/Camera_Data/all_species/quin_cam_act_2022_long.csv")
# 
# write_csv(quin_act_all_final_wide, "data/Camera_Data/all_species/quin_cam_act_2022_wide.csv")
# 
# write_csv(quin_dets_all_final2, "data/Camera_Data/all_species/quin_detections_all_species_2022.csv")
# 
# write_csv(dep_key2, "data/Camera_Data/all_species/quin_deployment_key_2022.csv")
# 





