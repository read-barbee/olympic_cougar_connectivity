#### LEKT Camera Activity Long Form All Years ####

# Author: Read Barbee

# Date:2024-10-16 

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

################################ libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)

###############################################################################

#1. Import viewshed files

###############################################################################

#2019 
vs_2019 <- read_csv("data/Camera_Data/all_species/LEKT/2019/LEKT_2019_viewsheds_with_camtypes.csv") %>% 
  clean_names() %>% 
  select(station, sum_area) %>% 
  rename(station_id = station,
         viewshed = sum_area)

#2020
vs_2020 <- read_csv("data/Camera_Data/all_species/LEKT/2020/lekt_2020_viewsheds_edited.csv") %>% 
  rename(station_id = `2020.Sta.Num`, 
         viewshed = Total_Area,
         exclude = Exclude) %>%
  select(station_id, viewshed, exclude) %>% 
  filter(exclude == FALSE) %>% 
  mutate(station_id = str_remove(station_id, "_.*")) %>% 
  select(-exclude)
  
vs_2021 <- read_csv("data/Camera_Data/all_species/LEKT/2021/lekt_2021_viewsheds_edited.csv") %>% 
  rename(station_id = `2021.Sta.Num`, 
         viewshed = Total_Area) %>%
  select(station_id, viewshed) %>%
  mutate(camera_id = case_when(str_detect(station_id, coll("b")) ~ "Camera2",
                               .default = "Camera1"), .after = station_id) %>% 
  mutate(station_id = str_remove(station_id, "a")) %>% 
  mutate(station_id = str_remove(station_id, "b")) %>% 
  mutate(station_id = paste0("Station", station_id))


vs_2022 <- read_csv("data/Camera_Data/all_species/LEKT/2022/lekt_2022_viewsheds_edited.csv") %>% 
  rename(station_id = cam,
         viewshed = area)

vs_2023 <- read_csv("data/Camera_Data/all_species/LEKT/2023/lekt_2023_viewsheds_edited.csv") %>% 
  rename(station_id = `2023.Sta.Num`,
         viewshed = Total_Area)

###############################################################################

#2. Format individual activity sheets

###############################################################################

################################ 2019 #################################
lekt_act_2019 <- read_csv("data/Camera_Data/all_species/LEKT/2019/2019.Activity.Sheet_corrected.csv")

#Create a key to translate between old and new deployment names
loc_key_2019 <- lekt_act_2019 %>% 
  clean_names() %>% 
  select(station, cameras, longitude, latitude) %>% 
  rename(station_id = station,
         camera_id = cameras) %>% 
  mutate(year = "2019", .before = station_id)

#check the number of distinct locations (for number of output rows)
target_rows_2019 <- lekt_act_2019 %>% 
  distinct(longitude, latitude) %>% 
  nrow()

#format column names and merge rows by location
lekt_act_2019_f <- lekt_act_2019 %>%
  select(-c(study, relative_location, camera_id, cam_num)) %>% 
  rename(station_id = station,
         camera_id = cameras,
         set_date = set_up_date) %>% 
  left_join(vs_2019, by = "station_id") %>%
  relocate(viewshed, .after = camera_id) %>%
  mutate(station_id =  str_replace(station_id,  "(\\D)(\\d)\\b", "\\10\\2")) %>% 
  mutate(station_id = str_replace(station_id, coll("Station"), coll("LEKT"))) %>% 
  relocate(longitude, .before = latitude) %>% 
  mutate(set_date = mdy(set_date),
         pull_date = mdy(pull_date)) %>%
merge_by_location() 
  
  
#check if the correct number of rows were returned
nrow(lekt_act_2019_f) == target_rows_2019
  

#pivot the merged and formatted activity matrix to long form for joining with other years. Add relevant columns for joining
lekt_act_2019_long <- lekt_act_2019_f %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = as.character(mdy(act_date))) %>% 
  mutate(dep_year = "2019", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>% 
  mutate(survey_id = "LEKT_2019", .before = station_id) %>% 
  mutate(set_date= as.character(set_date),
         pull_date = as.character(pull_date))


################################ 2020 #################################

lekt_act_2020 <-read_csv("data/Camera_Data/all_species/LEKT/2020/2020.Camera.Activity.Sheet.IDS.copy.csv")

loc_key_2020 <- lekt_act_2020 %>% 
  select(station, cameras, x, y) %>% 
  rename(station_id = station,
         camera_id = cameras,
         longitude = x,
         latitude = y) %>% 
  mutate(year = "2020", .before = station_id)

#check the number of distinct locations (for number of output rows)
target_rows_2020 <- lekt_act_2020 %>%
  distinct(x, y) %>% 
  nrow()

lekt_act_2020_f <- lekt_act_2020 %>% 
  select(-c(study, camera_id, time_cam_on, elevation)) %>%
  rename(longitude = x,
         latitude = y,
         station_id = station,
         camera_id = cameras,
         set_date = set_up_date) %>% 
  left_join(vs_2020, by = "station_id") %>%
  relocate(viewshed, .after = camera_id) %>% 
  mutate(station_id =  str_replace(station_id,  "(\\D)(\\d)\\b", "\\10\\2")) %>% 
  mutate(station_id = str_replace(station_id, coll("Station"), coll("LEKT"))) %>%
  relocate(viewshed, .after = camera_id) %>%
  relocate(longitude, .before = latitude) %>% 
  mutate(pull_date = case_when(pull_date == "Logged" ~ NA,
                               .default = pull_date)) %>% 
  mutate(set_date = mdy(set_date),
         pull_date = mdy(pull_date)) %>%
  merge_by_location()
  
  
#check if the correct number of rows were returned
nrow(lekt_act_2020_f) == target_rows_2020


lekt_act_2020_long <- lekt_act_2020_f %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>%
  mutate(act_date = as.character(dmy(act_date))) %>%
  mutate(dep_year = "2020", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>%
  mutate(survey_id = "LEKT_2020", .before = station_id) %>% 
  mutate(set_date = as.character(set_date),
         pull_date = as.character(pull_date))
  
  
################################ 2021 #################################

lekt_act_2021 <-read_csv("data/Camera_Data/all_species/LEKT/2021/2021.Camera.Activity.Sheet.csv")

loc_key_2021 <- lekt_act_2021 %>% 
  select(station, cameras, x, y) %>% 
  rename(station_id = station,
         camera_id = cameras,
         longitude = x,
         latitude = y) %>% 
  mutate(year = "2021", .before = station_id)

#check the number of distinct locations (for number of output rows)
target_rows_2021 <- lekt_act_2021 %>%
  distinct(x, y) %>% 
  nrow()


lekt_act_2021_f <- lekt_act_2021 %>% 
  select(-c(study, camera_id, elevation)) %>% 
  rename(longitude = x,
         latitude = y,
         station_id=station,
         camera_id = cameras) %>%
  left_join(vs_2021, by = c("station_id", "camera_id")) %>%
  relocate(viewshed, .after = camera_id) %>% 
  mutate(station_id =  str_replace(station_id,  "(\\D)(\\d)\\b", "\\10\\2")) %>% 
  mutate(station_id = str_replace(station_id, coll("Station"), coll("LEKT"))) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>%
  mutate(act_date = mdy(act_date))


start_end_dates_2021 <- lekt_act_2021_f %>% 
start_end_from_matrix() 
  

lekt_act_2021_wide <- lekt_act_2021_f %>% 
  left_join(start_end_dates_2021, by = "station_id") %>% 
  relocate(set_date, pull_date, .after = viewshed) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  merge_by_location()


  #check if the correct number of rows were returned
  nrow(lekt_act_2021_wide) == target_rows_2021

lekt_act_2021_long <- lekt_act_2021_wide %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(dep_year = "2021", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>%
  mutate(survey_id = "LEKT_2021", .before = station_id) %>% 
  mutate(act_date = as.character(act_date),
         set_date = as.character(set_date),
         pull_date = as.character(pull_date)) 


################################ 2022 #################################
lekt_act_2022 <-read_csv("data/Camera_Data/all_species/LEKT/2022/2022.Activity.Sheet.NEW.csv")

loc_key_2022 <- lekt_act_2022 %>%
  clean_names() %>% 
  select(station, cameras, x, y) %>% 
  rename(station_id = station,
         camera_id = cameras,
         longitude = x,
         latitude = y) %>% 
  mutate(year = "2022", .before = station_id)

#check the number of distinct locations (for number of output rows)
target_rows_2022 <- lekt_act_2022 %>% 
  distinct(X, Y) %>% 
  nrow()

lekt_act_2022_f <- lekt_act_2022 %>%
  clean_names() %>% 
  select(-c(study, camera_id, elevation)) %>% 
  rename(longitude = x,
         latitude = y,
         station_id=station,
         camera_id = cameras) %>%
  left_join(vs_2022, by = "station_id") %>% 
  relocate(viewshed, .after = camera_id) %>% 
  mutate(station_id =  str_replace(station_id,  "(\\D)(\\d)\\b", "\\10\\2")) %>% 
  mutate(station_id = str_replace(station_id, coll("Station"), coll("LEKT"))) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(act_date = mdy(act_date)) %>% 
  relocate(camera_id, .after = station_id) 
    
  start_end_dates_2022 <- lekt_act_2022_f %>%
    start_end_from_matrix() 
  
  
  lekt_act_2022_wide <- lekt_act_2022_f %>% 
    left_join(start_end_dates_2022, by = "station_id") %>%
    relocate(set_date, pull_date, .after = viewshed) %>% 
    pivot_wider(names_from = act_date, values_from = status) %>% 
    merge_by_location()
  
  
  #check if the correct number of rows were returned
  nrow(lekt_act_2022_wide) == target_rows_2022
  
  lekt_act_2022_long <- lekt_act_2022_wide %>% 
    select(-rows_merged) %>% 
    pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
    mutate(dep_year = "2022", 
           act_year = year(ymd(act_date)),
           .before = act_date) %>%
    mutate(survey_id = "LEKT_2022", .before = station_id) %>% 
    mutate(act_date = as.character(act_date),
           set_date = as.character(set_date),
           pull_date = as.character(pull_date))
  
################################ 2023 #################################
lekt_act_2023 <- read_csv("data/Camera_Data/all_species/LEKT/2023/2023.Activity Sheet Only_edited.csv") 
  
loc_key_2023 <- lekt_act_2023 %>% 
  clean_names() %>% 
    select(cam, x, y) %>% 
    rename(station_id = cam,
           longitude = x,
           latitude = y) %>% 
  mutate(camera_id = case_when(str_detect(station_id, "_2") ~ "Camera2",
                               .default = "Camera1"), .after = station_id) %>% 
  mutate(station_id = paste0("Station", station_id)) %>%
  mutate(station_id = str_remove(station_id, "_.*")) %>% 
    mutate(year = "2023", .before = station_id)

#check the number of distinct locations (for number of output rows)
target_rows_2023 <- lekt_act_2023 %>% 
  distinct(X, Y) %>% 
  nrow()
  
lekt_act_2023_wide <- lekt_act_2023 %>% 
  clean_names() %>% 
  select(-c(viewshed, start_time:check3, end_time, end)) %>%
  rename(longitude = x,
         latitude = y,
         station_id=cam,
         set_date = start_date,
         pull_date = end_date) %>%
  left_join(vs_2023, by = "station_id") %>% 
  relocate(viewshed, .after = pull_date) %>% 
  mutate(station_id = paste0("LEKT", station_id)) %>%
  mutate(station_id = str_remove(station_id, "_.*")) %>%
  mutate(station_id =  str_replace(station_id,  "(\\D)(\\d)\\b", "\\10\\2")) %>% 
  mutate(set_date = mdy(set_date),
         pull_date = mdy(pull_date)) %>% 
  merge_by_location()
  
  #check if the correct number of rows were returned
  nrow(lekt_act_2023_wide) == target_rows_2023


lekt_act_2023_long <- lekt_act_2023_wide %>%   
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(dep_year = "2023", 
         act_year = year(mdy(act_date)),
         .before = act_date) %>%
  mutate(survey_id = "LEKT_2023", .before = station_id) %>% 
  mutate(act_date = as.character(mdy(act_date)),
         set_date = as.character(set_date),
         pull_date = as.character(pull_date)) 


###############################################################################

#3. Combine activity sheets in long and wide formats

###############################################################################

lekt_act_all <- bind_rows(lekt_act_2019_long,
                          lekt_act_2020_long,
                          lekt_act_2021_long,
                          lekt_act_2022_long,
                          lekt_act_2023_long)


loc_keys_all <- bind_rows(loc_key_2019,
                          loc_key_2020,
                          loc_key_2021,
                          loc_key_2022,
                          loc_key_2023) %>% 
  rename(dep_year = year)


lekt_act_all_long <- lekt_act_all %>% 
  unite("deployment_id", survey_id, station_id,  sep="_", remove = F) %>%  #camera_id,
  relocate(survey_id, .before = deployment_id) %>% 
  select(-act_year)


#check for duplicates in deployment id and activity date
lekt_act_all_long  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date)


#calculate the start and end date for each deployment from the activity matrix.
#returns warnings if there is not final NA (i.e. 0s or 1s extend all the way to the end of the activity matrix for that year). In these cases, use the final date in the activity matrix for that deployment.

start_end_dates <- lekt_act_all_long %>%
  arrange(deployment_id, act_date) %>%
  group_by(deployment_id) %>%
  summarise(
    set_date = min(act_date[!is.na(status) & lag(is.na(status), default = TRUE)], na.rm = TRUE),
    pull_date = coalesce(
      max(act_date[is.na(status) & lag(!is.na(status), default = FALSE)], na.rm = TRUE),
      max(act_date, na.rm = TRUE)  # Return the last date if no transition
    )
  ) %>%
  ungroup()

#join and reformat columns
lekt_act_all_long2 <- lekt_act_all_long %>% 
  select(-c(set_date, pull_date)) %>%
  left_join(start_end_dates, by = "deployment_id") %>% 
  relocate(c(set_date, pull_date), .after = viewshed)

#fill missing dates (periods of NA between surveys that weren't included in annual activity matrices)
lekt_act_all_long_date_fill <- lekt_act_all_long2 %>% 
  pivot_wider(names_from = act_date, values_from = status) %>%
  fill_missing_dates(mode = "survey") %>% 
  pivot_longer(-c(survey_id:dep_year), names_to = "act_date", values_to = "status")
  


#check for duplicates in deployment id and activity date
lekt_act_all_long_date_fill  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date) #%>% distinct(deployment_id)

#lmit total number of active cameras per day at a location to 1 in case of overlap errors in matrix
lekt_act_all_final_long <- lekt_act_all_long_date_fill %>% 
  mutate(status = case_when(status > 1 ~ 1,
                            .default = status))

# #camera_level_test
lekt_act_all_final_long  %>%
  #filter(dep_year == "2019") %>% 
  arrange(dep_year, set_date) %>% 
  select(deployment_id, act_date, status) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  column_to_rownames("deployment_id") %>%
  as.matrix() %>%
  camtrapR:::camopPlot(., lattice = TRUE)


#pivot from long form to wide form
lekt_act_all_final_wide <- lekt_act_all_final_long %>%
  pivot_wider(names_from = act_date, values_from = status)
  
#check for duplicates in deployment_id
lekt_act_all_final_wide %>% get_dupes(deployment_id)


#Make a deployment key to link old detection metadata to new deployment structure
dep_key <- lekt_act_all_final_wide %>%
  select(deployment_id, longitude, latitude, dep_year) %>% 
  left_join(loc_keys_all, by = c("dep_year","longitude", "latitude")) 

#########################################################################
##
## 4. Combine detections for all species across years
##
##########################################################################

################################ 2019 #################################
lekt_det_2019 <- read_csv("data/Camera_Data/all_species/LEKT/2019/SNA_20190703_20200214_2022.09.19_21.23_metadata_tbl.csv") %>% 
  clean_names() %>% 
  select(station,cameras, date_time_original, species) %>% 
  rename(timestamp = date_time_original,
         station_id = station,
         camera_id = cameras) %>% 
  mutate(timestamp = as.character(mdy_hm(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2019", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)


################################ 2020 #################################
lekt_det_2020 <- read_csv("data/Camera_Data/all_species/LEKT/2020/S3071_20200513_20201210_2022.12.12_22.30_metadata_tbl.V2.csv") %>% 
  clean_names() %>% 
  select(station, cameras, date_time_original, species) %>% 
  rename(timestamp = date_time_original,
         station_id = station,
         camera_id = cameras) %>% 
  filter(is.na(station_id) == F ) %>% 
  mutate(timestamp = as.character(mdy_hm(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2020", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)


################################ 2021 #################################
lekt_det_2021 <- read_csv("data/Camera_Data/2021/LEKT_2021/S3071_20210419_20211209_2022.09.12_20.04_metadata_tbl.csv", col_types = cols(DateTimeOriginal = col_character())) %>% 
  clean_names() %>% 
  select(station, cameras, date_time_original, species) %>% 
  rename(timestamp = date_time_original,
         station_id = station,
         camera_id = cameras) %>% 
  filter(!is.na(station_id)) %>% 
  mutate(timestamp = as.character(ymd_hms(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2021", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)

################################ 2022 #################################
lekt_det_2022 <- read_csv("data/Camera_Data/all_species/LEKT/2022/2022 OCP Camera Data.csv") %>% 
  clean_names() %>% 
  select(station_num, timestamp, common_name) %>% 
  rename(station_id = station_num,
         species = common_name) %>%
  mutate(camera_id = "Camera1", .after = station_id) %>% 
  filter(!is.na(station_id)) %>% 
  mutate(station_id = str_remove_all(station_id, "Check.*")) %>% 
  mutate(timestamp = as.character(mdy_hm(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2022", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)


################################ 2023 #################################
lekt_det_2023 <- read_csv("data/Camera_Data/all_species/LEKT/2023/2023.Data_key_spp_time_date_shifts_complete.csv",
                          col_types = cols(timestamp = col_character())) %>% 
  rename(station_id = Station,
         species = common_name) %>% 
  select(station_id, timestamp, species) %>% 
  filter(!is.na(station_id)) %>% 
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F) %>% 
  mutate(station_id = paste0("Station", station_id),
         camera_id = "Camera1") %>% 
  mutate(dep_year = "2023", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)


#fix issue with matching detections to stations 33a and b due to lack of camera ids in original activity file.
start_2023_33_a <- lekt_act_all_final_wide %>% 
  filter(deployment_id == "LEKT_2023_LEKT33_a") %>% pull(set_date)

end_2023_33_a <- lekt_act_all_final_wide %>% 
  filter(deployment_id == "LEKT_2023_LEKT33_a") %>% pull(pull_date)

start_2023_33_b<- lekt_act_all_final_wide %>% 
  filter(deployment_id == "LEKT_2023_LEKT33_b") %>% pull(set_date)

end_2023_33_b<- lekt_act_all_final_wide %>% 
  filter(deployment_id == "LEKT_2023_LEKT33_b") %>% pull(pull_date)


lekt_det_2023_fix <- lekt_det_2023 %>% 
  mutate(camera_id = case_when(station_id == "Station33" & 
                                 timestamp > start_2023_33_a &
                                 timestamp < end_2023_33_a  ~ "Camera1",
                               station_id == "Station33" & 
                                 timestamp > start_2023_33_b &
                                 timestamp < end_2023_33_b  ~ "Camera2",
                               .default = camera_id))


################################ Combine years #################################
lekt_dets_all <- bind_rows(lekt_det_2019,
                           lekt_det_2020,
                           lekt_det_2021,
                           lekt_det_2022,
                           lekt_det_2023_fix) 

###############################################################################

#5. Change detection deployment names to match activity history

###############################################################################

#Rename the deployment and station_ids in the detection file 
lekt_dets_all_final <- lekt_dets_all %>% 
  mutate(camera_id = case_when(station_id == "Station43" & 
                                 camera_id == "Camera2" &
                                 dep_year == "2019" ~ "Camera3",
                               .default = camera_id)) %>% 
  left_join(dep_key, by = c("dep_year", "station_id", "camera_id")) %>% 
  select(-c(station_id, camera_id)) %>% 
  select(deployment_id, longitude, latitude, everything()) %>% 
  mutate(time = case_when(is.na(time) ~ "00:00:00",
                          .default = time)) %>% 
  mutate(timestamp = paste0(date, " ", time))

#get rid of one human record before recorded start date at one deployment
lekt_dets_all_final <- lekt_dets_all_final %>% 
  filter(!(deployment_id == "LEKT_2019_LEKT33" &
             species == "Human" &
             timestamp < ymd("2019-07-18"))) %>% 
  mutate(act_year = as.character(act_year))
  
  
#lekt_dets_all_final %>% mutate(timestamp = ymd_hms(timestamp))

#check for mismatches between deployment ids in detection and activty history frames
det_stations <- lekt_dets_all_final %>% distinct(deployment_id) %>% pull(deployment_id)

act_stations <- lekt_act_all_final_long %>% distinct(deployment_id) %>% pull(deployment_id)

setdiff(act_stations, det_stations) #stations in act not in det
setdiff(det_stations, act_stations) #stations in det not in act

#Stations with no data: 
#LEKT_2021_LEKT65- stolen/not deployed
#LEKT_2022_LEKT92 - no entries in detection file
#LEKT_2023_LEKT33_a - stolen/not deployed


###############################################################################

#6. Test the camtrapR detectionHistory function using the new frames

###############################################################################


lekt_act_all_final_wide_mat <- lekt_act_all_final_wide %>% 
  select(-c(survey_id, station_id:dep_year)) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

test <- detectionHistory(lekt_dets_all_final,
                         "Puma",
                         lekt_act_all_final_wide_mat,
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

colnames(det_hist) <- colnames(lekt_act_all_final_wide_mat)

camtrapR:::camopPlot(det_hist, lattice = TRUE)


###############################################################################

#7. Try to stack the detection history

###############################################################################


test2 <- det_hist %>%
  as.data.frame() %>% 
  rownames_to_column("deployment_id") %>%
  pivot_longer(-deployment_id, names_to = "act_date", values_to = "detection") %>% 
  mutate(act_date = ymd(act_date),
         jday = yday(act_date)) %>% 
  mutate(act_year = as.factor(year(act_date))) 


test3 <- test2 %>% 
  split(.$act_year) %>% 
  map(., function(x){
    x %>% select(-c(act_year, act_date)) %>% 
      pivot_wider(names_from = jday, values_from = detection)
  }) %>% 
 map(., fill_jdays)
  
for(i in 1:length(test3)){
  year_name <- names(test3)[i]
  
  test3[[i]] <- test3[[i]] %>% 
    filter(str_detect(deployment_id, pattern = year_name))
}

test4 <- test3 %>% 
  bind_rows() #%>% 
  #distinct(deployment_id, .keep_all = T) 


  #filter(!all(is.na(c_across(-deployment_id))))


date_seq <- seq(mdy("01/01/2020"), mdy("12/31/2020"), by = "day") %>% as.character()

plot_mat <- test4 %>% 
column_to_rownames("deployment_id") %>% 
  as.matrix()

colnames(plot_mat) <- date_seq

camtrapR:::camopPlot(plot_mat, lattice = TRUE)

###############################################################################

#8. Export

###############################################################################

# write_csv(lekt_act_all_final_long, "data/Camera_Data/all_species/lekt_cam_act_2019_2023_long.csv")
# #
# write_csv(lekt_act_all_final_wide, "data/Camera_Data/all_species/lekt_cam_act_2019_2023_wide.csv")
# 
#write_csv(lekt_dets_all_final, "data/Camera_Data/all_species/lekt_detections_all_species_2019-2023.csv")
# 
# write_csv(dep_key, "data/Camera_Data/all_species/lekt_deployment_key_2019-2023.csv")

