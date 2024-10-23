#### Olympic National Park Detection History Formatting ####

# Author: Read Barbee

# Date:2023-09-11 
#Last updated: 2023-10-13

# Purpose:
#updated so that unknown bait, snare and camera days are treated as the first days of an interval

################################ Libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)
library(terra)
library(sf)


################################ Helper functions #################################



#test_dat <- dat %>% filter(station_id=="603_r1_s1_2013")

#function to combine visit periods into single row for each station_year
format_station <- function(station_dat, column) {
  name <- as.name(column)
  
  round_1 <- station_dat %>%
    mutate(problem = case_when(interval != !!name ~ TRUE,
                               interval == !!name ~ FALSE), .after = camera_days_good) %>% 
    mutate(int_start = visit_date - days(interval), .before = visit_date) %>% 
    #mutate(problem_days = interval - !!name, .after = camera_days_good) %>% 
    rename(int_end = visit_date) %>% 
    mutate(problem_start = case_when(problem == TRUE ~ int_start + days(!!name),
                                     problem == FALSE ~ NA),
           problem_end = case_when(problem == TRUE ~ int_end, #- days(interval - !!name)
                                   problem == FALSE ~ NA), .after = problem)
  
  if(nrow(station_dat) == 1){
    
    station_dat_formatted <- round_1 %>% 
      mutate(visit_num = 1) %>% 
      pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
      reframe(.by = station_id,
              year = first(year),
              utm_e = first(utm_e),
              utm_n = first(utm_n),
              hex_id = first(hex_id),
              station_num = first(station_num),
              set_date = min(int_start),
              pull_date = max(int_end),
              #problem_days = sum(problem_days),
              problem_periods = sum(problem),
              problem_start_1 = min(problem_start_1, na.rm = T),
              problem_end_1 = max(problem_end_1, na.rm = T),
              problem_start_2 = NA,
              problem_end_2 = NA,
              problem_start_3 = NA,
              problem_end_3 =NA,
              problem_start_4 = NA,
              problem_end_4 = NA,
              camera_days_good = sum(camera_days_good, na.rm = T),
              bait_days_good = sum(bait_days_good, na.rm = T),
              snare_days_good = sum(snare_days_good, na.rm = T),
              #effort_correction_factor = sum(effort_correction_factor, na.rm = T),
              photo_count = sum(photo_count, na.rm = T),
              cougar = sum(cougar)
      ) %>% 
      mutate_all(~ replace(., is.infinite(.), NA))
  }
  
  if(nrow(station_dat) == 2){
    
    station_dat_formatted <- round_1 %>% 
      mutate(visit_num = 1:nrow(station_dat)) %>% 
      pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
      reframe(.by = station_id,
              year = first(year),
              utm_e = first(utm_e),
              utm_n = first(utm_n),
              hex_id = first(hex_id),
              station_num = first(station_num),
              set_date = min(int_start),
              pull_date = max(int_end),
             # problem_days = sum(problem_days),
              problem_periods = sum(problem),
              problem_start_1 = min(problem_start_1, na.rm = T),
              problem_end_1 = max(problem_end_1, na.rm = T),
              problem_start_2 = min(problem_start_2, na.rm = T),
              problem_end_2 = max(problem_end_2, na.rm = T),
              problem_start_3 = NA,
              problem_end_3 =NA,
              problem_start_4 = NA,
              problem_end_4 = NA,
              camera_days_good = sum(camera_days_good, na.rm = T),
              bait_days_good = sum(bait_days_good, na.rm = T),
              snare_days_good = sum(snare_days_good, na.rm = T),
              #effort_correction_factor = sum(effort_correction_factor, na.rm = T),
              photo_count = sum(photo_count, na.rm = T),
              cougar = sum(cougar)
      ) %>% 
      mutate_all(~ replace(., is.infinite(.), NA))
  }
  
  if(nrow(station_dat) == 3){
    
    station_dat_formatted <- round_1 %>% 
      mutate(visit_num = 1:nrow(station_dat)) %>% 
      pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>%
      
      reframe(.by = station_id,
              year = first(year),
              utm_e = first(utm_e),
              utm_n = first(utm_n),
              hex_id = first(hex_id),
              station_num = first(station_num),
              set_date = min(int_start),
              pull_date = max(int_end),
              #problem_days = sum(problem_days),
              problem_periods = sum(problem),
              problem_start_1 = min(problem_start_1, na.rm = T),
              problem_end_1 = max(problem_end_1, na.rm = T),
              problem_start_2 = min(problem_start_2, na.rm = T),
              problem_end_2 = max(problem_end_2, na.rm = T),
              problem_start_3 = min(problem_start_3, na.rm = T),
              problem_end_3 = max(problem_end_3, na.rm = T),
              problem_start_4 = NA,
              problem_end_4 = NA,
              camera_days_good = sum(camera_days_good, na.rm = T),
              bait_days_good = sum(bait_days_good, na.rm = T),
              snare_days_good = sum(snare_days_good, na.rm = T),
              #effort_correction_factor = sum(effort_correction_factor, na.rm = T),
              photo_count = sum(photo_count, na.rm = T),
              cougar = sum(cougar)
      ) %>% 
      mutate_all(~ replace(., is.infinite(.), NA))
  }
  
  if(nrow(station_dat) == 4){
    
    station_dat_formatted <- round_1 %>% 
      mutate(visit_num = 1:nrow(station_dat)) %>% 
      pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
      reframe(.by = station_id,
              year = first(year),
              utm_e = first(utm_e),
              utm_n = first(utm_n),
              hex_id = first(hex_id),
              station_num = first(station_num),
              set_date = min(int_start),
              pull_date = max(int_end),
              #problem_days = sum(problem_days),
              problem_periods = sum(problem),
              problem_start_1 = min(problem_start_1, na.rm = T),
              problem_end_1 = max(problem_end_1, na.rm = T),
              problem_start_2 = min(problem_start_2, na.rm = T),
              problem_end_2 = max(problem_end_2, na.rm = T),
              problem_start_3 = min(problem_start_3, na.rm = T),
              problem_end_3 = max(problem_end_3, na.rm = T),
              problem_start_4 = min(problem_start_4, na.rm = T),
              problem_end_4 = max(problem_end_4, na.rm = T),
              camera_days_good = sum(camera_days_good, na.rm = T),
              bait_days_good = sum(bait_days_good, na.rm = T),
              snare_days_good = sum(snare_days_good, na.rm = T),
              #effort_correction_factor = sum(effort_correction_factor, na.rm = T),
              photo_count = sum(photo_count, na.rm = T),
              cougar = sum(cougar)
      ) %>% 
      mutate_all(~ replace(., is.infinite(.), NA))
  }
  
  
  return(station_dat_formatted)
}

#function to format yearly data by station_id
format_station_by_year <- function(year_dat, column){
  round_1 <- year_dat %>%  
    mutate(station_id = as.character(station_id)) %>% 
    mutate(station_id = as.factor(station_id))
  dat_split <- split(round_1, round_1$station_id)
  formatted <- dat_split %>% map(format_station, column)
  formatted <- bind_rows(formatted)
  return(formatted)
  }


prep_matrix <- function(dat, column){
  #split data into list elements by year to iterate over 
  dat_split <- split(dat, dat$year)
  
  #format each year by station_id
  dat_format <- dat_split %>% map(format_station_by_year, column)
  
  #bind the formatted list elements into a single dataframe
  dat_format <- bind_rows(dat_format)
  
  # Load your sf object with points for station locations
  sf_points <- dat_format %>% st_as_sf(coords=c("utm_e", "utm_n"), crs = 26910, remove=FALSE) %>% 
    st_transform(crs = 5070)
  
  # Extract cell numbers for each station location
  xy <- sf_points %>% st_coordinates()
  
  dat_format <- dat_format %>%
    mutate(cell_id = cellFromXY(temp, xy), .after = utm_n)
  
  return(dat_format)
}

prep_matrix_part2 <- function(part1){
  
 
  #format dataframe for captrapR::cameraOperation function 
  dat_format2 <- part1 %>% 
    rename(Problem1_from = problem_start_1,
           Problem1_to = problem_end_1,
           Problem2_from = problem_start_2,
           Problem2_to = problem_end_2,
           Problem3_from = problem_start_3,
           Problem3_to = problem_end_3,
           Problem4_from = problem_start_4,
           Problem4_to = problem_end_4,) %>% 
    select(station_id, cell_id, set_date, pull_date, Problem1_from:Problem4_to) %>% 
    mutate(across(c(set_date:Problem4_to), as.character))
  
  return(dat_format2)
}

#function to convert all values to be 0, 1, or NA
fix_camop <- function(col){
  col <- case_when(
    col < 0 ~ 0,
    col > 1 ~ 1,
    col > 0 & col < 1 ~ 1,
    .default = col
  )
}


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


#########################################################################
##
## 1. Import and format detection data 
##
##########################################################################

dat <- read_csv("data/Camera_Data/Olympic_National_Park/raw/ONP_fisher_grid_2013-2016_raw_formatted.csv") %>% 
  clean_names() %>% 
  select(-c(rep_id, study_incid, waypoint_type)) %>% 
  mutate(cougar = case_when(cougar >1 ~ 1,
                             .default = cougar),
         visit_date = mdy(visit_date)) %>%  #station_id = str_replace(station_id, coll("."), coll("_"))
  mutate(station_id = as.factor(station_id),
         station = paste0(hex_id, "_", station_num),
         camera_days_good = case_when(camera_days_good > interval ~ interval,
                                      .default = camera_days_good)) %>% 
  mutate(across(c(bait_days_good, snare_days_good), \(x) replace_na(x, 0))) %>% 
  mutate(across(camera_days_good:snare_days_good, \(x) case_when(x > interval ~ interval,
                                                                 .default = x))) %>% 
  select(station_id, hex_id, station_num, station, everything(), -rep)


# Load raster data
#raster <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_asp.tif")

#temp <- raster[[1]]



#########################################################################
##
## 2. Summarize start and end dates for each year and output for later
##
##########################################################################

survey_periods <- dat %>% 
  group_by(year) %>% 
  summarize(min_date = min(visit_date), 
            max_date = max(visit_date))

#write_csv(survey_periods, "data/Camera_Data/Olympic_National_Park/onp_survey_periods_2013_2016.csv")


###############################################################################

#1. Import activity matrix that I made in 2023

###############################################################################
cam_act <- read_csv("data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_cam_act_10-11-23.csv")

cam_act %>% get_dupes(station_id)


#convert to long form format for merging with OCP data

cam_act_f1 <- cam_act %>%
  pivot_longer(-c(station_id:snare_days_good), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = ymd(act_date),
         viewshed = NA) 

start_end_dates <- cam_act_f1 %>% start_end_from_matrix() %>% 
  mutate(pull_date = case_when(is.infinite(pull_date) ~NA,
                               .default = pull_date))

cam_act_f2 <- cam_act_f1 %>% 
  left_join(start_end_dates, by = "station_id") %>% 
  relocate(set_date, pull_date, .after = viewshed) %>% 
  mutate(act_year = year(act_date)) %>% 
  mutate(deployment_id = paste0("ONP_", hex_id, "_", station_num, "_", year), .before = station_id) %>% 
  mutate(station_id = paste0("ONP_", hex_id, "_", station_num)) %>% 
  rename(dep_year = year) %>% 
  mutate(survey_id = paste0("ONP_", as.character(dep_year)), .before = deployment_id) %>% 
  select(survey_id, deployment_id, station_id, utm_e, utm_n, viewshed, set_date, pull_date, dep_year, act_year, act_date, status) 
  

cam_act_sf <- cam_act_f2 %>% 
  st_as_sf(coords = c(x = "utm_e", y = "utm_n"), crs = 26910) %>% 
  st_transform(crs = 4326)

cam_act_sf %>% slice_sample(n=100) %>% mapview::mapview()

cam_act_lat_long <- st_coordinates(cam_act_sf) %>% as.data.frame()


cam_act_f3 <- cam_act_f2 %>% 
  mutate(longitude = cam_act_lat_long$X,
         latitude = cam_act_lat_long$Y) %>% 
 select(survey_id, deployment_id, station_id, longitude, latitude, viewshed, set_date, pull_date, dep_year, act_date, status)

# cam_act_f4 <- cam_act_f3 %>%
#   distinct(longitude, latitude, dep_year, .keep_all = T)

#fill missing dates (periods of NA between surveys that weren't included in annual activity matrices)
cam_act_long_final <- cam_act_f3 %>% 
  pivot_wider(names_from = act_date, values_from = status) %>%
  fill_missing_dates(mode = "survey") %>% 
  pivot_longer(-c(survey_id:dep_year), names_to = "act_date", values_to = "status")

cam_act_plot_obj <- cam_act_long_final %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  arrange(dep_year, set_date) %>% 
  select(-c(survey_id, station_id:dep_year)) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

#plot camera operation matrix
camtrapR:::camopPlot(camOp = cam_act_plot_obj, palette = "Heat", lattice = TRUE)


cam_act_wide <- cam_act_long_final %>% 
  pivot_wider(names_from = act_date, values_from = status) 

#########################################################################
##
## 3. Construct daily camera activity matrix from 2 week interval data
##
##########################################################################

#alternatively, import existing activity matrix
 # cam_act <- read_csv("data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_cam_act_10-11-23.csv")

# #prep data
# act_mat <- prep_matrix(dat, "camera_days_good")
# 
# act_mat2 <- prep_matrix_part2(act_mat)

#convert the dataframe into a daily actiity matrix
# camop <- camtrapR::cameraOperation(CTtable = act_mat2,
#                                    stationCol = "station_id",
#                                    setupCol = "set_date",
#                                    retrievalCol = "pull_date",
#                                    writecsv = FALSE,
#                                    hasProblems = TRUE,
#                                    dateFormat = "ymd"
# )
# 
# 
# #fix negative and non-integer values using the above function. 
# camop_fixed <- camop %>% as.data.frame() %>% 
#   mutate(across(everything(), fix_camop)) %>% 
#   as.matrix()


#plot camera operation matrix
#camtrapR:::camopPlot(camOp = cam_act_plot_obj, palette = "Heat", lattice = TRUE)
# 
# 
# #rejoin columns dropped by camtrapR funtction
# cols_to_join <- act_mat %>% select(station_id:station_num, camera_days_good, bait_days_good, snare_days_good)
# 
# cam_act_final <- camop_fixed  %>% as.data.frame() %>% rownames_to_column("station_id") %>% 
#   left_join(cols_to_join, by = join_by(station_id)) %>%
#   select(station_id, hex_id, station_num, year, utm_e, utm_n, cell_id, everything(), -c(camera_days_good, bait_days_good, snare_days_good)) 
#camera_days_good, bait_days_good, snare_days_good, 

#write_csv(cam_act_final, "data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_cam_act_10-11-23.csv")

# act_mat_long <- cam_act_final  %>% 
#   pivot_longer(cols= -c(station_id:cell_id), names_to = "date", values_to = "cam_status") %>% 
#   mutate(date = ymd(date)) %>% 
#   unite("station_id", c(hex_id, station_num, year), sep = "_", remove = TRUE)

#########################################################################
##
## 4.  Construct daily bait matrix from 2 week interval data
##
##########################################################################

#prep data
# bait_mat <- prep_matrix(dat, "bait_days_good")
# 
# bait_mat2 <- prep_matrix_part2(bait_mat)
# 
# #convert the dataframe into a daily actiity matrix
# camop_bait <- camtrapR::cameraOperation(CTtable = bait_mat2,
#                                    stationCol = "station_id",
#                                    setupCol = "set_date",
#                                    retrievalCol = "pull_date",
#                                    writecsv = FALSE,
#                                    hasProblems = TRUE,
#                                    dateFormat = "ymd"
# )
# 
# 
# #fix negative and non-integer values using the above function. 
# camop_bait_fixed <- camop_bait %>% as.data.frame() %>% 
#   mutate(across(everything(), fix_camop)) %>% 
#   as.matrix()
# 
# #rejoin columns dropped by camtrapR funtction
# cols_to_join_bait <- bait_mat %>% select(station_id:station_num, bait_days_good, snare_days_good)
# 
# bait_mat_final <- camop_bait_fixed  %>% as.data.frame() %>% rownames_to_column("station_id") %>% 
#   left_join(cols_to_join_bait, by = join_by(station_id)) %>%
#   select(station_id, hex_id, station_num, year, utm_e, utm_n, cell_id, everything(), -c(bait_days_good, snare_days_good))
# 
# 
# #write_csv(bait_mat_final, "data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_bait_mat_10-11-23.csv")
# 
# bait_mat_long <- bait_mat_final  %>% 
#   pivot_longer(cols= -c(station_id:cell_id), names_to = "date", values_to = "bait_status") %>% 
#   mutate(date = ymd(date)) %>% 
#   unite("station_id", c(hex_id, station_num, year), sep = "_", remove = TRUE) %>% 
#   select(station_id, date, bait_status)

#########################################################################
##
## 5.  Construct daily snare matrix from 2 week interval data
##
##########################################################################
# 
# #prep data
# snare_mat <- prep_matrix(dat, "snare_days_good")
# 
# snare_mat2 <- prep_matrix_part2(snare_mat)
# 
# #convert the dataframe into a daily actiity matrix
# camop_snare <- camtrapR::cameraOperation(CTtable = snare_mat2,
#                                         stationCol = "station_id",
#                                         setupCol = "set_date",
#                                         retrievalCol = "pull_date",
#                                         writecsv = FALSE,
#                                         hasProblems = TRUE,
#                                         dateFormat = "ymd"
# )
# 
# 
# #fix negative and non-integer values using the above function. 
# camop_snare_fixed <- camop_snare %>% as.data.frame() %>% 
#   mutate(across(everything(), fix_camop)) %>% 
#   as.matrix()
# 
# #rejoin columns dropped by camtrapR funtction
# cols_to_join_snare <- snare_mat %>% select(station_id:station_num, bait_days_good, snare_days_good)
# 
# snare_mat_final <- camop_snare_fixed  %>% as.data.frame() %>% rownames_to_column("station_id") %>% 
#   left_join(cols_to_join_snare, by = join_by(station_id)) %>%
#   select(station_id, hex_id, station_num, year, utm_e, utm_n, cell_id, everything(), -c(bait_days_good, snare_days_good))
# 
# 
# #write_csv(snare_mat_final, "data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_snare_mat_10-11-23.csv")
# 
# snare_mat_long <- snare_mat_final  %>% 
#   pivot_longer(cols= -c(station_id:cell_id), names_to = "date", values_to = "snare_status") %>% 
#   mutate(date = ymd(date)) %>% 
#   unite("station_id", c(hex_id, station_num, year), sep = "_", remove = TRUE) %>% 
#   select(station_id, date, snare_status)
# 


#########################################################################
##
## 6. Import and format species detections
##
##########################################################################

species_det <- read_csv("data/Camera_Data/Olympic_National_Park/raw/raw/onp_species_det_2013_2016.csv") %>% 
  clean_names()


coords_to_join <- cam_act_f3 %>% 
  distinct(deployment_id, .keep_all = T) %>% 
  mutate(dep_year = as.character(dep_year)) %>% 
  select(deployment_id, longitude, latitude)


species_det_f <- species_det %>% 
  mutate(across(hex_id:station_num, as.character)) %>% 
  mutate(station_num = str_remove(station_num, coll(".0"))) %>% 
  mutate(deployment_id = paste0("ONP_", hex_id, "_",station_num,"_", rep_year),
         station_id = paste0("ONP_", hex_id, "_", station_num),
         dep_year = rep_year) %>% 
  filter(!(is.na(start_date) & is.na(start_time))) %>%
  mutate(start_time = as.character(start_time)) %>% 
  mutate(start_time = case_when(is.na(start_time) ~ "00:00:00",
                                .default = start_time)) %>% 
  mutate(timestamp = mdy_hms(paste0(start_date, " ", start_time))) %>%
  mutate(act_year = year(timestamp)) 

species_det_f2 <- species_det_f %>% 
  left_join(coords_to_join, by = c("deployment_id")) %>% 
  select(deployment_id, station_id, longitude, latitude, dep_year, act_year, species, timestamp) 


species_det_f2 %>% filter(is.na(latitude))

species_det_final <- species_det_f2 %>% 
  filter(!(deployment_id=="ONP_171_1_2013" & timestamp > ymd("2013-11-05")))

###############################################################################

#6. Test the camtrapR detectionHistory function using the new frames

###############################################################################


onp_act_all_final_wide_mat <- cam_act_wide %>% 
  arrange(dep_year, set_date) %>% 
  select(-c(survey_id, station_id:dep_year)) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

test <- detectionHistory(species_det_final,
                         "Black_tailed_deer",
                         onp_act_all_final_wide_mat,
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

colnames(det_hist) <- colnames(onp_act_all_final_wide_mat)

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


date_seq <- seq(mdy("01/01/2016"), mdy("12/31/2016"), by = "day") %>% as.character()

plot_mat <- test4 %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix() 

colnames(plot_mat) <- date_seq

camtrapR:::camopPlot(plot_mat, lattice = TRUE)

###############################################################################

#8. Export

###############################################################################

# write_csv(cam_act_long_final, "data/Camera_Data/all_species/onp_cam_act_2013-2016_long.csv")
# #
# write_csv(cam_act_wide, "data/Camera_Data/all_species/onp_cam_act_2013-2016_wide.csv")
# 
# write_csv(species_det_final, "data/Camera_Data/all_species/onp_detections_all_species_2013-2016.csv")



