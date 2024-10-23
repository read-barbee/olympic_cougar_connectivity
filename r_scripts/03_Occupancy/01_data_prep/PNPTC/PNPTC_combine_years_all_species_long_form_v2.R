#### PNPTC combine years all species ####

# Author: Read Barbee

# Date:2023-09-11 

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
    ungroup() %>% 
    rename_with(.cols=-c(!matches("\\d")), function(x){as.character(date_seq)}) 
  
  return(date_matrix)
}

################################ libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)

###############################################################################

#1. Import files with station coordinates

###############################################################################

locs_2020 <- read_csv("data/Camera_Data/all_species/PNPTC/2020/pnptc_2020_station_locations.csv")

locs_2021 <- read_csv("data/Camera_Data/all_species/PNPTC/2021/pnptc_2021_station_locations.csv")

locs_2022 <- read_csv("data/Camera_Data/all_species/PNPTC/2022/pnptc_2022_station_locations_edited.csv")

locs_2022_d <- locs_2022 %>% 
  rename(station_id = CameraID,
         camera_id = Camera,
         latitude = Latitude,
         longitude = Longitude) %>% 
  mutate(station_id= case_when(camera_id == "Camera2" ~ paste0(station_id, "_2"),
                             .default = station_id)) %>% 
  distinct(station_id, .keep_all = T)
  

###############################################################################

#2. Format individual activity sheets

###############################################################################

################################ 2020 #################################
pnptc_act_2020 <- read_csv("data/Camera_Data/all_species/PNPTC/2020/2020_DeployDates.csv")

pnptc_act_2020_locs <- pnptc_act_2020 %>% 
  mutate(start_date = as.character(mdy_hms(paste0(startdate, " ", starttime))),
         end_date = as.character(mdy_hms(paste0(enddate, " ", endtime)))) %>% 
  rename(station_id = cam,
         viewshed = area) %>% 
  select(station_id, start_date, end_date, viewshed) %>% 
  mutate(station_id = case_when(station_id == "FS53V2" ~  "FS53_b",
                                station_id == "FS53" ~ "FS53_a",
                                .default = station_id)) %>%
  left_join(locs_2020, by = "station_id") %>% 
  relocate(c(lon, lat), .before = viewshed) %>% 
  rename(longitude = lon,
         latitude = lat)

#Create a key to translate between old and new deployment names
loc_key_2020 <- pnptc_act_2020_locs %>% 
  clean_names() %>% 
  select(station_id, longitude, latitude) %>%
  mutate(station_id = case_when(station_id == "FS53_b" ~ "FS53V2",
                                station_id == "FS53_a" ~ "FS53",
                                .default = station_id)) %>%
  mutate(dep_year = "2020", .before = station_id) %>% 
  distinct(longitude, latitude, .keep_all = T)



#check the number of distinct locations (for number of output rows)
target_rows_2020 <- pnptc_act_2020_locs %>% 
  distinct(longitude, latitude) %>% 
  nrow()

#format column names and merge rows by location
pnptc_act_2020_f <- pnptc_act_2020_locs %>% 
  mutate(start_date = date(ymd_hms(start_date)),
         end_date = date(ymd_hms(end_date))) %>% 
  dep_to_matrix() %>%
  rename(set_date = start_date,
         pull_date = end_date) %>% 
  merge_by_location() %>% 
  pivot_longer(-c(!matches("\\d")), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = ymd(act_date)) %>% 
  group_by(station_id, latitude, longitude) %>% 
  mutate(status = case_when(act_date < date(set_date) ~ NA,
                            act_date > date(pull_date) ~ NA,
                            .default = status)) %>% 
  ungroup() %>%
  pivot_wider(names_from = act_date, values_from = status) %>% 
  mutate(station_id = paste0("PNPTC_", station_id))


# pnptc_act_2020_f %>%
#   arrange(set_date) %>%
# select(c(station_id, matches("\\d"))) %>% column_to_rownames("station_id") %>%   as.matrix() %>% camtrapR:::camopPlot(lattice=TRUE)


#check if the correct number of rows were returned
nrow(pnptc_act_2020_f) == target_rows_2020


#pivot the merged and formatted activity matrix to long form for joining with other years. Add relevant columns for joining
pnptc_act_2020_long <- pnptc_act_2020_f %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = as.character(ymd(act_date))) %>% 
  mutate(dep_year = "2020", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>% 
  mutate(survey_id = "PNPTC_2020", .before = station_id) %>% 
  mutate(set_date= as.character(set_date),
         pull_date = as.character(pull_date))

################################ 2021 #################################
#using the activity sheet sent previously since new deployment file is missing location and date info
pnptc_act_2021 <- read_csv("data/Camera_Data/all_species/PNPTC/2021/2021_CameraActivitySheet_edited.csv")


#Create a key to translate between old and new deployment names
loc_key_2021 <- pnptc_act_2021 %>% 
  clean_names() %>% 
  select(camera_id, longitude, latitude) %>% 
  rename(station_id = camera_id) %>% 
  mutate(dep_year = "2021", .before = station_id) %>% 
  distinct(longitude, latitude, .keep_all = T)

#check the number of distinct locations (for number of output rows)
target_rows_2021 <- pnptc_act_2021 %>% 
  clean_names() %>% 
  distinct(longitude, latitude) %>% 
  nrow()

#format column names and merge rows by location
pnptc_act_2021_f <- pnptc_act_2021 %>%
  clean_names() %>% 
  select(-c(country, study, elevation_m)) %>% 
  rename(station_id = camera_id,
         camera_id = camera) %>% 
  mutate(station_id = case_when(camera_id == "Camera2" ~ paste0(station_id, "_2"),
                                .default = station_id)) %>% 
  mutate(viewshed = NA,
         station_id = paste0("PNPTC_", station_id)) %>% 
  relocate(viewshed, .after = camera_id) %>%
  relocate(longitude, .before = latitude) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>%
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(act_date = mdy(act_date)) %>% 
  relocate(camera_id, .after = station_id) 
  
#warnings here if activity extends all the way to the end of the matrix 
start_end_dates_2021 <- pnptc_act_2021_f %>%
  start_end_from_matrix() %>% 
  mutate(pull_date = case_when(is.infinite(pull_date) ~ max(pnptc_act_2021_f$act_date),
                               .default = pull_date))


pnptc_act_2021_wide <- pnptc_act_2021_f %>% 
  left_join(start_end_dates_2021, by = "station_id") %>%
  relocate(set_date, pull_date, .after = viewshed) %>%  
  pivot_wider(names_from = act_date, values_from = status) %>% 
  mutate(station_id = str_remove(station_id, coll("_2"))) %>%
  merge_by_location()

#check if the correct number of rows were returned
nrow(pnptc_act_2021_wide) == target_rows_2021

#Internal NA values for FS02 actually encoded in the activity sheet dylan sent me
# pnptc_act_2021_wide %>%
#   arrange(set_date) %>%
# select(c(station_id, matches("\\d"))) %>% column_to_rownames("station_id") %>%   as.matrix() %>% camtrapR:::camopPlot(lattice=TRUE)


pnptc_act_2021_long <- pnptc_act_2021_wide %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(dep_year = "2021", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>%
  mutate(survey_id = "PNPTC_2021", .before = station_id) %>% 
  mutate(act_date = as.character(act_date),
         set_date = as.character(set_date),
         pull_date = as.character(pull_date)) 


  

################################ 2022 #################################
pnptc_act_2022 <- read_csv("data/Camera_Data/all_species/PNPTC/2022/2022_DeployDates.csv")


pnptc_act_2022_locs <- pnptc_act_2022 %>% 
  mutate(start_date = as.character(mdy_hms(paste0(startdate, " ", starttime))),
         end_date = as.character(mdy_hms(paste0(enddate, " ", endtime)))) %>% 
  rename(station_id = cam,
         viewshed = area) %>% 
  select(station_id, start_date, end_date, viewshed) %>% 
  left_join(locs_2022_d, by = "station_id") %>% 
  relocate(c(longitude, latitude), .before = viewshed) %>% 
  select(-camera_id) %>% 
  mutate(longitude = replace_na(longitude, -999), #replace missing location with -999 until I find the location
         latitude = replace_na(latitude, -999))

pnptc_act_2022_locs %>% filter(is.na(longitude)) 

#Create a key to translate between old and new deployment names
loc_key_2022 <- pnptc_act_2022_locs %>% 
  clean_names() %>% 
  select(station_id, longitude, latitude) %>% 
  mutate(dep_year = "2022", .before = station_id) %>% 
  distinct(longitude, latitude, .keep_all = T)



#check the number of distinct locations (for number of output rows)
target_rows_2022 <- pnptc_act_2022_locs %>% 
  distinct(longitude, latitude) %>% 
  nrow()

#format column names and merge rows by location
pnptc_act_2022_f <- pnptc_act_2022_locs %>% 
  mutate(start_date = date(ymd_hms(start_date)),
         end_date = date(ymd_hms(end_date))) %>% 
  dep_to_matrix() %>%
  rename(set_date = start_date,
         pull_date = end_date) %>% 
  merge_by_location() %>%
  pivot_longer(-c(!matches("\\d")), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = ymd(act_date)) %>% 
  group_by(station_id, latitude, longitude) %>% 
  mutate(status = case_when(act_date < date(set_date) ~ NA,
                            act_date > date(pull_date) ~ NA,
                            .default = status)) %>% 
  ungroup() %>%
  pivot_wider(names_from = act_date, values_from = status) %>%
  mutate(station_id = paste0("PNPTC_", station_id))

#weird that FS55 extends so far into 2023, but it has detection records to back it up
# pnptc_act_2022_f %>%
#   arrange(set_date) %>%
# select(c(station_id, matches("\\d"))) %>% column_to_rownames("station_id") %>%   as.matrix() %>% camtrapR:::camopPlot(lattice=TRUE)


#check if the correct number of rows were returned
nrow(pnptc_act_2022_f) == target_rows_2022


#pivot the merged and formatted activity matrix to long form for joining with other years. Add relevant columns for joining
pnptc_act_2022_long <- pnptc_act_2022_f %>% 
  select(-rows_merged) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = as.character(ymd(act_date))) %>% 
  mutate(dep_year = "2022", 
         act_year = year(ymd(act_date)),
         .before = act_date) %>% 
  mutate(survey_id = "PNPTC_2022", .before = station_id) %>% 
  mutate(set_date= as.character(set_date),
         pull_date = as.character(pull_date))

###############################################################################

#3. Combine activity sheets in long and wide formats

###############################################################################

pnptc_act_all <- bind_rows(pnptc_act_2020_long,
                           pnptc_act_2021_long,
                           pnptc_act_2022_long)


loc_keys_all <- bind_rows(loc_key_2020,
                          loc_key_2021,
                          loc_key_2022) 


pnptc_act_all_long <- pnptc_act_all %>% 
  unite("deployment_id", survey_id, station_id,  sep="_", remove = F) %>%  #camera_id,
  relocate(survey_id, .before = deployment_id) %>% 
  select(-act_year)


#check for duplicates in deployment id and activity date
pnptc_act_all_long  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date)


#fill missing dates (periods of NA between surveys that weren't included in annual activity matrices)
pnptc_act_all_long_date_fill <- pnptc_act_all_long %>% 
  pivot_wider(names_from = act_date, values_from = status) %>%
  fill_missing_dates(mode = "survey") %>% 
  pivot_longer(-c(survey_id:dep_year), names_to = "act_date", values_to = "status")



#check for duplicates in deployment id and activity date
pnptc_act_all_long_date_fill  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date) #%>% distinct(deployment_id)

#lmit total number of active cameras per day at a location to 1 in case of overlap errors in matrix
pnptc_act_all_final_long <- pnptc_act_all_long_date_fill %>% 
  mutate(status = case_when(status > 1 ~ 1,
                            .default = status)) %>% 
  filter(ymd(act_date) <= ymd("2022-11-22"))
  
  
  # mutate(act_date = case_when(ymd(act_date) > ymd("2022-11-22") ~ "2022-11-22",
  #                             .default = act_date))

# #camera_level_test
pnptc_act_all_final_long  %>%
  #filter(dep_year == "2021") %>% 
  arrange(dep_year, set_date) %>% 
  select(deployment_id, act_date, status) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  column_to_rownames("deployment_id") %>%
  as.matrix() %>%
  camtrapR:::camopPlot(., lattice = TRUE)


#pivot from long form to wide form
pnptc_act_all_final_wide <- pnptc_act_all_final_long %>%
  pivot_wider(names_from = act_date, values_from = status)

#check for duplicates in deployment_id
pnptc_act_all_final_wide %>% get_dupes(deployment_id)


#Make a deployment key to link old detection metadata to new deployment structure
dep_key <- pnptc_act_all_final_wide %>%
  select(deployment_id, longitude, latitude, dep_year) %>% 
  left_join(loc_keys_all, by = c("dep_year","longitude", "latitude")) 


#########################################################################
##
## 4. Combine detections for all species across years
##
##########################################################################

################################ 2020 #################################
pnptc_det_2020 <- read_csv("data/Camera_Data/all_species/PNPTC/2020/2020_AllPhotos_Clean.csv") 
  

pnptc_det_2020_f <- pnptc_det_2020 %>% 
  pivot_longer(-c(File:Time), names_to = "species", values_to = "status") %>% 
  filter(status ==1) %>% 
  unite("timestamp", Date, Time, sep = " ", remove = F) %>% 
  rename(station_id = File,
         date = Date,
         time = Time) %>% 
  mutate(timestamp = mdy_hms(timestamp),
         date = as.character(mdy(date))) %>% 
  mutate(dep_year = "2020",
         act_year = year(timestamp),
         timestamp = as.character(timestamp)) %>% 
  select(station_id, dep_year, act_year, species, timestamp, date, time)

################################ 2021 #################################
pnptc_det_2021 <- read_csv("data/Camera_Data/all_species/PNPTC/2021/2021_AllPhotos_NoSexOrAge.csv") 


pnptc_det_2021_f <- pnptc_det_2021 %>% 
  unite("timestamp", Date, Time, sep = " ", remove = F) %>% 
  rename(station_id = RelativePath,
         date = Date,
         time = Time,
         species = Species) %>% 
  mutate(timestamp = mdy_hms(timestamp),
         date = as.character(mdy(date))) %>% 
  mutate(dep_year = "2021",
         act_year = year(timestamp),
         timestamp = as.character(timestamp)) %>% 
  select(station_id, dep_year, act_year, species, timestamp, date, time)


################################ 2022 #################################
pnptc_det_2022 <- read_csv("data/Camera_Data/all_species/PNPTC/2022/AllData_2022_Clean.csv") 


pnptc_det_2022_f <- pnptc_det_2022 %>% 
  pivot_longer(-c(Cam:Time), names_to = "species", values_to = "status") %>% 
  filter(status ==1) %>% 
  unite("timestamp", Date, Time, sep = " ", remove = F) %>% 
  rename(station_id = Cam,
         date = Date,
         time = Time) %>% 
  mutate(timestamp = mdy_hms(timestamp),
         date = as.character(mdy(date))) %>% 
  mutate(dep_year = "2022",
         act_year = year(timestamp),
         timestamp = as.character(timestamp)) %>% 
  select(station_id, dep_year, act_year, species, timestamp, date, time)

################################ 2023 #################################
pnptc_det_2023 <- read_csv("data/Camera_Data/all_species/PNPTC/2023/2023_Images_Clean.csv") 


pnptc_det_2023_f <- pnptc_det_2023 %>% 
  pivot_longer(-c(Site:Time), names_to = "species", values_to = "status") %>% 
  filter(status ==1) %>% 
  unite("timestamp", Date, Time, sep = " ", remove = F) %>% 
  rename(station_id = Site,
         date = Date,
         time = Time) %>% 
  mutate(timestamp = mdy_hms(timestamp),
         date = as.character(mdy(date))) %>% 
  mutate(dep_year = "2022",
         act_year = year(timestamp),
         timestamp = as.character(timestamp)) %>% 
  select(station_id, dep_year, act_year, species, timestamp, date, time)


################################ Combine years #################################
pnptc_dets_all <- bind_rows(pnptc_det_2020_f,
                           pnptc_det_2021_f,
                           pnptc_det_2022_f) #not adding 2023 yet until I get the activity file


###############################################################################

#5. Change detection deployment names to match activity history

###############################################################################

dep_key2 <- dep_key %>% 
  mutate(station_id = case_when(deployment_id == "PNPTC_2021_PNPTC_FS66_b" ~ "FS66.V2",
         .default = station_id))
  

FS66_a_start_2021 <- "2021-06-04"
FS66_a_end_2021 <- "2021-11-18"

FS66_b_start_2021 <- ymd("2021-08-20")
FS66_b_end_2021 <- "2021-11-18"

#Rename the deployment and station_ids in the detection file. FS66 is weird
pnptc_dets_all_final <- pnptc_dets_all %>%
  filter(!str_detect(station_id, coll(".JPG"))) %>% 
  filter(station_id != "ElkCam03") %>%
  filter(!(station_id == "FS05" & dep_year == "2022")) %>%
  left_join(dep_key2, by = c("dep_year", "station_id")) %>%
  mutate(date = case_when(deployment_id == "PNPTC_2021_PNPTC_FS41" & 
                          act_year == "2020" ~ as.character(ymd(date) + years(1)),
                          .default = date)) %>% 
  mutate(timestamp = paste0(date, " ", time)) %>% 
  mutate(act_year = year(ymd_hms(timestamp))) %>% 
  select(deployment_id, longitude, latitude, dep_year, act_year, species, timestamp, date, time)

pnptc_dets_all_final <- pnptc_dets_all_final %>% 
  filter(!(deployment_id== "PNPTC_2022_PNPTC_FS55" & act_year == "2023"))

pnptc_dets_all_final %>% filter(deployment_id == "PNPTC_2022_PNPTC_FS55") %>% View()


pnptc_dets_all_final %>% mutate(timestamp = ymd_hms(timestamp)) %>% filter(is.na(timestamp))

#check for mismatches between deployment ids in detection and activty history frames
det_stations <- pnptc_dets_all_final %>% distinct(deployment_id) %>% pull(deployment_id)

act_stations <- pnptc_act_all_final_long %>% distinct(deployment_id) %>% pull(deployment_id)

setdiff(act_stations, det_stations) #stations in act not in det
setdiff(det_stations, act_stations) #stations in det not in act

#"PNPTC_2020_PNPTC_DNR93" is in the detection file but legitimately has no species detections in the file. 



###############################################################################

#6. Test the camtrapR detectionHistory function using the new frames

###############################################################################


pnptc_act_all_final_wide_mat <- pnptc_act_all_final_wide %>%
  arrange(set_date) %>% 
  select(-c(survey_id, station_id:dep_year)) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

test <- detectionHistory(pnptc_dets_all_final,
                         "Deer",
                         pnptc_act_all_final_wide_mat,
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

colnames(det_hist) <- colnames(pnptc_act_all_final_wide_mat)

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

# write_csv(pnptc_act_all_final_long, "data/Camera_Data/all_species/pnptc_cam_act_2020_2022_long.csv")
# 
# write_csv(pnptc_act_all_final_wide, "data/Camera_Data/all_species/pnptc_cam_act_2020_2022_wide.csv")
# 
# write_csv(pnptc_dets_all_final, "data/Camera_Data/all_species/pnptc_detections_all_species_2020-2022.csv")
# 
# write_csv(dep_key2, "data/Camera_Data/all_species/pnptc_deployment_key_2010-2022.csv")
# 






################################ Graveyard #################################
# 
# 
# select(deployment_id, station_id, longitude, latitude, everything()) %>% get_dupes(unique_id) %>% View()
# 
# pnptc_dets_all2 <- pnptc_dets_all %>% 
#   mutate(unique_id = 1:nrow(.), .before = station_id,
#          timestamp = ymd_hms(paste0(date, " ", time)))
# 
# 
# #fuzzy join by date: note--this removes about 7,000 records
# test <- pnptc_act_all_final_wide %>%
#   mutate(set_date = ymd(set_date),
#          pull_date = ymd(pull_date)) %>% 
#   select(deployment_id:pull_date) %>% 
#   fuzzyjoin::regex_left_join(pnptc_dets_all2, by = c("station_id" = "station_id")) %>%
#   filter(timestamp >= set_date & timestamp <= pull_date) %>% 
#   select(deployment_id, station_id.x, station_id.y, longitude, latitude, dep_year, act_year, species, timestamp, date, time) %>% 
#   rename(station_id_new = station_id.x,
#          station_id_old = station_id.y)


