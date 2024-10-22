#### Format SKOK Detection Histories All Species ####

# Author: Read Barbee

# Date:2024-09-13 

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


#make activity matrix from Wildlife Insights format deployment dates
dep_to_matrix <- function(data){
  
  min_date <- min(data$start_date)
  max_date <- max(data$end_date)
  
  date_seq <- seq(min_date, max_date, by = "day")
  
  date_matrix <- data %>%
    rowwise() %>%
    mutate(
      # Create a logical vector: TRUE if the date is between start and end, otherwise FALSE
      date_vector = list(as.numeric(date_seq >= start_date & date_seq <= end_date))
    ) %>%
    unnest_wider(date_vector, names_sep = "_") %>%
    ungroup()
  
  # Replace 0 with NA
  date_matrix[date_matrix == 0] <- NA
  
  #rename_columns
  date_matrix <- date_matrix %>% rename_with(.cols=-c(deployment_id:end_date), function(x){as.character(date_seq)})
  
  return(date_matrix)
}

################################ Libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)

###############################################################################

#2. Import viewshed files

###############################################################################

vs_2021 <- read_csv("data/Camera_Data/all_species/SKOK/2021/skok_viewsheds_2021_edited.csv") %>% 
  rename(set_date = start,
         pull_date = end)


vs_2022 <- read_csv("data/Camera_Data/all_species/SKOK/2022_2023/skok_viewsheds_2022.csv") %>% 
  rename(station_id = cam,
         viewshed = area) %>% 
  distinct(station_id, viewshed, .keep_all = F) %>% 
  mutate(dep_year = "2022", .after = station_id)


#########################################################################
##
## 3. Combine activity sheets for all years
##
##########################################################################

################################ 2021 #################################
skok_act_2021 <- read_csv("data/Camera_Data/all_species/SKOK/2021/Act_sheet_2021_final.csv")

#Create a key to translate between old and new deployment names
loc_key_2021 <- skok_act_2021 %>% 
  clean_names() %>% 
  select(station, cameras, x, y) %>%
  rename(station_id = station,
         camera_id = cameras,
         longitude = x,
         latitude = y) %>% 
  mutate(dep_year = "2021", .before = station_id)

#check the number of distinct locations (for number of output rows)
target_rows_2021 <- skok_act_2021 %>% 
  distinct(X, Y) %>% 
  nrow()

#initial formatting
skok_act_2021_f <- skok_act_2021 %>% 
  clean_names() %>%
  select(-c(study, camera_id, elevation)) %>% 
  rename(longitude = x,
         latitude = y,
         station_id=station,
         camera_id = cameras) %>% 
  left_join(vs_2021, by = "station_id") %>% 
  relocate(c(set_date, pull_date, viewshed), .after = camera_id) %>%
  mutate(station_id =  str_replace(station_id,  "(\\D)(\\d)\\b", "\\10\\2")) %>%
  merge_by_location() 
  
#check if the correct number of rows were returned
nrow(skok_act_2021_f) == target_rows_2021


skok_act_2021_long <- skok_act_2021_f %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(act_date = mdy(act_date),
         set_date = mdy_hm(set_date),
         pull_date = mdy_hm(pull_date)) %>% 
  group_by(station_id, longitude, latitude) %>% 
  mutate(status = case_when(act_date < date(set_date) ~ NA,
                            act_date > date(pull_date) ~ NA,
                            .default = status)) %>% 
  ungroup() %>% 
  relocate(longitude, .before = latitude) %>% 
  mutate(dep_year = "2021", .before = act_date) %>% 
  mutate(set_date = as.character(date(set_date)),
         pull_date = as.character(date(pull_date)),
         act_date = as.character(act_date)) %>% 
  mutate(survey_id = "SKOK_2021", .before = station_id) 



################################ 2022 #################################

#setdiff(vs_2022$station_id, skok_act_2022_2023$station_id)

skok_act_2022_2023 <- read_csv("data/Camera_Data/all_species/SKOK/2022_2023/deployments_fixed.csv") 

#parse deployment id column into station and year for easier use
skok_act_2022_2023 <- skok_act_2022_2023 %>% 
  separate_wider_delim(deployment_id, delim = "-", names = c("station_id", "dep_year"), cols_remove = F) %>% 
  mutate(dep_year = str_remove(dep_year, "[A-Za-z]")) %>% 
  relocate(deployment_id, .before = station_id)


#Create a key to translate between old and new deployment names
loc_key_2022_2023 <- skok_act_2022_2023 %>% 
  clean_names() %>% 
  select(deployment_id, longitude, latitude) %>% 
  mutate(dep_year = case_when(str_detect(deployment_id, coll("2022")) ~ "2022",
                              str_detect(deployment_id, coll("2023")) ~ "2023"),
         .before = longitude)

#separate out the 2022 entries
skok_act_2022 <- skok_act_2022_2023 %>% 
  filter(dep_year == "2022")

#check the number of distinct locations (for number of output rows)
target_rows_2022 <- skok_act_2022 %>% 
  distinct(longitude, latitude) %>% 
  nrow()


#format column names and merge rows by location
skok_act_2022_wide <- skok_act_2022 %>% 
  clean_names() %>% 
  select(deployment_id:end_date, -placename) %>% 
  mutate(start_date = mdy_hm(start_date),
         end_date = mdy_hm(end_date)) %>% 
  dep_to_matrix() %>% 
  rename(set_date = start_date,
         pull_date = end_date) %>% 
  mutate(station_id = str_replace_all(station_id, coll("CAM"), coll("Station"))) %>% 
  left_join(vs_2022, by = c("station_id", "dep_year")) %>%
  relocate(c(set_date, pull_date, viewshed), .after = latitude) %>% 
  merge_by_location()


#check the number of distinct locations (for number of output rows)
nrow(skok_act_2022_wide) == target_rows_2022


#pivot the merged and formatted activity matrix to long form for joining with other years. Add relevant columns for joining
skok_act_2022_long <- skok_act_2022_wide %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(dep_year = "2022",
         act_date = as.character(act_date),
         set_date = as.character(set_date),
         pull_date = as.character(pull_date),
         survey_id = paste0("SKOK_", dep_year)) %>% 
  select(survey_id, station_id, dep_year, longitude, latitude, set_date, pull_date, viewshed, act_date, status) 

################################ 2023 #################################

#separate out the 2022 entries
skok_act_2023 <- skok_act_2022_2023 %>% 
  filter(dep_year == "2023")

#check the number of distinct locations (for number of output rows)
target_rows_2023 <- skok_act_2023 %>% 
  distinct(longitude, latitude) %>%
  nrow()


#format column names and merge rows by location
skok_act_2023_wide <- skok_act_2023 %>% 
  clean_names() %>% 
  select(deployment_id:end_date, -placename) %>% 
  mutate(start_date = mdy_hm(start_date),
         end_date = mdy_hm(end_date)) %>% 
  dep_to_matrix() %>% 
  rename(set_date = start_date,
         pull_date = end_date) %>% 
  mutate(station_id = str_replace_all(station_id, coll("CAM"), coll("Station"))) %>% 
  mutate(viewshed = NA) %>% 
  relocate(c(set_date, pull_date, viewshed), .after = latitude) %>% 
  merge_by_location()


#check the number of distinct locations (for number of output rows)
nrow(skok_act_2023_wide) == target_rows_2023


#pivot the merged and formatted activity matrix to long form for joining with other years. Add relevant columns for joining
skok_act_2023_long <- skok_act_2023_wide %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(dep_year = "2023",
         act_date = as.character(act_date),
         set_date = as.character(set_date),
         pull_date = as.character(pull_date),
         survey_id = paste0("SKOK_", dep_year)) %>% 
  select(survey_id, station_id, dep_year, longitude, latitude, set_date, pull_date, viewshed, act_date, status) 


###############################################################################

#3. Combine activity sheets for all years in long and wide formats

###############################################################################
skok_act_all <- bind_rows(skok_act_2021_long,
                          skok_act_2022_long,
                          skok_act_2023_long) 

loc_keys_all <- bind_rows(loc_key_2021,
                          loc_key_2022_2023)



skok_act_all_long <- skok_act_all %>% 
  mutate(station_id = str_replace(station_id, coll("Station"), coll("SKOK"))) %>% 
  unite("deployment_id", survey_id, station_id, sep="_", remove = F) %>% 
  relocate(survey_id, .before = deployment_id) %>% 
  relocate(viewshed, .after = latitude) %>% 
  arrange(dep_year, station_id, set_date)

#check for duplicates in deployment id and activity date
skok_act_all_long  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date) 


#fill missing dates (periods of NA between surveys that weren't included in annual activity matrices)
skok_act_all_long_date_fill <- skok_act_all_long %>% 
  pivot_wider(names_from = act_date, values_from = status) %>%
  fill_missing_dates(mode = "survey") %>% 
  pivot_longer(-c(survey_id:dep_year), names_to = "act_date", values_to = "status")


#lmit total number of active cameras per day at a location to 1 in case of overlap errors in matrix
skok_act_all_final_long <- skok_act_all_long_date_fill %>% 
  mutate(status = case_when(status > 1 ~ 1,
                            .default = status))

# #camera_level_test
skok_act_all_final_long  %>%
  #filter(dep_year == "2023") %>% 
  arrange(dep_year, station_id, act_date) %>% 
  select(deployment_id, act_date, status) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  column_to_rownames("deployment_id") %>%
  as.matrix() %>%
  camtrapR:::camopPlot(., lattice = TRUE)


#pivot from long form to wide form
skok_act_all_final_wide <- skok_act_all_final_long %>%
  pivot_wider(names_from = act_date, values_from = status)

#check for duplicates in deployment_id
skok_act_all_final_wide %>% get_dupes(deployment_id)


#Make a deployment key to link old detection metadata to new deployment structure

dep_key_2021 <- skok_act_all_final_wide %>% 
  filter(dep_year == "2021") %>% 
  select(deployment_id, longitude, latitude, dep_year) %>% 
  left_join(loc_key_2021, by = c("dep_year","longitude", "latitude")) %>% 
  mutate(station_id =  str_replace(station_id,  "(\\D)(\\d)\\b", "\\10\\2")) 

dep_key_2022_2023 <- skok_act_all_final_wide %>%
  filter(dep_year != "2021") %>%
  select(deployment_id, longitude, latitude, dep_year) %>% 
  left_join(loc_key_2022_2023, by = c("dep_year","longitude", "latitude")) %>% 
  rename(deployment_id_new = deployment_id.x,
         deployment_id_old = deployment_id.y) %>% 
  select(deployment_id_new, deployment_id_old, longitude, latitude)

#########################################################################
##
## 1. Combine detections for all species across years. Change deployment names to match activity matrix
##
##########################################################################

################################ 2021 #################################
skok_det_2021 <- read_csv("data/Camera_Data/all_species/SKOK/2021/SKOK.2021.Metadata.Complete.modified.csv") %>% 
  clean_names() %>% 
  select(source_name_1, date_time_10, species) %>%
  separate_wider_delim(source_name_1, delim = ".", names = c("station_id", "check")) %>%
  mutate(station_id = str_replace_all(station_id, coll("CAM"), coll("Station")),
         camera_id = "Camera1") %>% 
  rename(timestamp = date_time_10) %>% 
  select(station_id, camera_id, timestamp, species) %>% 
  mutate(timestamp = mdy_hm(timestamp)) %>% 
  mutate(dep_year = "2021",
         act_year = year(timestamp),
         .after = camera_id) %>% 
  left_join(dep_key_2021, by = c("dep_year", "station_id", "camera_id")) %>% 
  mutate(timestamp = as.character(timestamp)) %>% 
  separate_wider_delim(timestamp, delim = " ", names = c("date", "time"), cols_remove = F, too_few = "align_start") %>%
  select(deployment_id, longitude, latitude, dep_year, act_year, species, timestamp, date, time) 


################################ 2021-2022 #################################


skok_det_2022_2023 <- read_csv("data/Camera_Data/all_species/SKOK/2022_2023/images_2006929.csv") %>% 
  select(deployment_id, timestamp, common_name) %>%  
  rename(deployment_id_old = deployment_id) %>% 
  mutate(deployment_id_old = case_when(deployment_id_old == "CAM38-20202" ~ "CAM38-2022",
                                       .default = deployment_id_old)) %>% 
  left_join(dep_key_2022_2023, by = "deployment_id_old") %>%
  rename(species = common_name,
         deployment_id = deployment_id_new) %>%
  mutate(dep_year = case_when(str_detect(deployment_id, "2022") ~ "2022",
                              str_detect(deployment_id, "2023") ~ "2023"),
         act_year = year(timestamp)) %>% 
  mutate(timestamp = as.character(timestamp)) %>% 
  separate_wider_delim(timestamp, delim = " ", names = c("date", "time"), cols_remove = F, too_few = "align_start") %>% 
  select(deployment_id, longitude, latitude, dep_year, act_year, species, timestamp, date, time)
  

################################ Combine years #################################
skok_det_all <- bind_rows(skok_det_2021,
                          skok_det_2022_2023) 


#fix issues with some timestamps
skok_det_all_final <- skok_det_all %>%
  mutate(timestamp = case_when(is.na(time) ~ paste0(timestamp, " 00:00:00"),
                               .default = timestamp)) %>% 
  mutate(time = case_when(is.na(time) ~ "00:00:00",
                               .default = time)) %>% 
  mutate(timestamp = ymd_hms(timestamp)) %>% 
  mutate(timestamp = case_when(deployment_id == "SKOK_2023_SKOK77" ~ timestamp + years(1),
         .default = timestamp)) %>% 
  mutate(timestamp = as.character(timestamp)) %>% 
  mutate(timestamp = case_when(is.na(timestamp) ~ paste0(date, " ", time),
                               .default = timestamp)) 

#workaround for parsing issue when timestamp is 00:00:00
# skok_det_all_final %>% mutate(timestamp = ymd_hms(as.POSIXct(timestamp, format = "%Y-%m-%d %T"))) %>% filter(is.na(timestamp)) %>% View()
  
# species_names <- skok_det_all_final %>% distinct(species) %>% pull(species)
# 
# test <- checkSpeciesNames(species_names, searchtype ="common")
#########################################################################
##
## 3. QA-QC
##
##########################################################################

det_stations <- skok_det_all %>% distinct(deployment_id) %>% pull(deployment_id)

act_stations <- skok_act_all_final_long %>% distinct(deployment_id) %>% pull(deployment_id)

setdiff(act_stations, det_stations) #stations in act not in det
setdiff(det_stations, act_stations) #stations in det not in act

#"SKOK_2021_SKOK07" -camera active but not in detection file
#"SKOK_2021_SKOK20" - camera active but not in detection file
#"SKOK_2021_SKOK34" - camera active but not in detection file
#"SKOK_2021_SKOK43" - camera active but not in detection file
#"SKOK_2021_SKOK49" - camera active but not in detection file
#"SKOK_2021_SKOK65" - camera active but not in detection file
#"SKOK_2022_SKOK39" - camera active but not in detection file. Images for 2023 deployment but not 2022
###############################################################################

#6. Test the camtrapR detectionHistory function using the new frames

###############################################################################


skok_act_all_final_wide_mat <- skok_act_all_final_wide %>% 
  select(-c(survey_id, station_id:dep_year)) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

test <- detectionHistory(skok_det_all_final,
                         "Puma",
                         skok_act_all_final_wide_mat,
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

colnames(det_hist) <- colnames(skok_act_all_final_wide_mat)

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
  bind_rows() 

date_seq <- seq(mdy("01/01/2023"), mdy("12/31/2023"), by = "day") %>% as.character()

plot_mat <- test4 %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

colnames(plot_mat) <- date_seq

camtrapR:::camopPlot(plot_mat, lattice = TRUE)

###############################################################################

#8. Export

###############################################################################

# write_csv(skok_act_all_final_long, "data/Camera_Data/skok_act_all_species_2021_2023_long.csv")
# 
# write_csv(skok_act_all_final_wide, "data/Camera_Data/skok_act_all_species_2021_2023_wide.csv")
# 
# 
# write_csv(skok_det_all_final, "data/Camera_Data/all_species/skok_detections_all_species_2019-2023.csv")
# 
# write_csv(dep_key_2021, "data/Camera_Data/all_species/skok_deployment_key_2021.csv")
# 
# write_csv(dep_key_2022_2023, "data/Camera_Data/all_species/skok_deployment_key_2022_2023.csv")










