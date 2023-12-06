#### Olympic National Park Detection History Formatting ####

# Author: Read Barbee

# Date:2023-09-11 
#Last updated: 2023-10-11

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)
library(terra)
library(sf)


################################ Helper functions #################################

#Function to merge images containing cougars that are fewer than 15 min apart
merge_events <- function(df) {
  df %>%
    arrange(station_id, date) %>%
    group_by(station_id) %>%
    mutate(time_diff = c(0, diff(date))) %>%
    mutate(grp = cumsum(time_diff >= 900)) %>%
    group_by(station_id,grp) %>%
    summarize(
      min_time = min(date),
      max_time = max(date),
      image_count = n()
    ) %>%
    mutate(date = date(min_time), .after=grp) %>% 
    ungroup() 
}


#genralized functions
format_station <- function(station_dat, column){
  
  name <- as.name(column)
  
  round_1 <- station_dat %>% 
    mutate(problem = case_when(interval != !!name ~ TRUE,
                               interval == !!name ~ FALSE), .after = camera_days_good) %>% 
    mutate(int_start = visit_date - days(interval), .before = visit_date) %>% 
    mutate(problem_days = interval - !!name, .after = camera_days_good) %>% 
    rename(int_end = visit_date) %>% 
    mutate(problem_start = case_when(problem == TRUE ~ int_start,
                                     problem == FALSE ~ NA),
           problem_end = case_when(problem == TRUE ~ int_end,
                                   problem == FALSE ~ NA), .after = problem,
           effort_correction_factor = case_when(problem == TRUE ~ !!name,
                                                problem == FALSE ~ 0))
  
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
              problem_days = sum(problem_days),
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
              effort_correction_factor = sum(effort_correction_factor, na.rm = T),
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
              problem_days = sum(problem_days),
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
              effort_correction_factor = sum(effort_correction_factor, na.rm = T),
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
              problem_days = sum(problem_days),
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
              effort_correction_factor = sum(effort_correction_factor, na.rm = T),
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
              problem_days = sum(problem_days),
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
              effort_correction_factor = sum(effort_correction_factor, na.rm = T),
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
  select(station_id, hex_id, station_num, station, everything(), -rep)


# Load raster data
raster <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_asp.tif")

temp <- raster[[1]]



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


#########################################################################
##
## 3. Construct daily camera activity matrix from 2 week interval data
##
##########################################################################

#alternatively, import existing activity matrix
 #cam_act <- read_csv("data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_cam_act_10-11-23.csv")

#prep data
act_mat <- prep_matrix(dat, "camera_days_good")

act_mat2 <- prep_matrix_part2(act_mat)

#convert the dataframe into a daily actiity matrix
camop <- camtrapR::cameraOperation(CTtable = act_mat2,
                                   stationCol = "station_id",
                                   setupCol = "set_date",
                                   retrievalCol = "pull_date",
                                   writecsv = FALSE,
                                   hasProblems = TRUE,
                                   dateFormat = "ymd"
)


#fix negative and non-integer values using the above function. 
camop_fixed <- camop %>% as.data.frame() %>% 
  mutate(across(everything(), fix_camop)) %>% 
  as.matrix()


#plot camera operation matrix
#camtrapR:::camopPlot(camOp = camop_fixed, palette = "Heat", lattice = TRUE)


#rejoin columns dropped by camtrapR funtction
cols_to_join <- act_mat %>% select(station_id:station_num, effort_correction_factor, camera_days_good, bait_days_good, snare_days_good)

cam_act_final <- camop_fixed  %>% as.data.frame() %>% rownames_to_column("station_id") %>% 
  left_join(cols_to_join, by = join_by(station_id)) %>%
  select(station_id, hex_id, station_num, year, utm_e, utm_n, cell_id, effort_correction_factor, camera_days_good, bait_days_good, snare_days_good, everything())

#write_csv(cam_act_final, "data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_cam_act_10-11-23.csv")

act_mat_long <- cam_act_final  %>% 
  pivot_longer(cols= -c(station_id:snare_days_good), names_to = "date", values_to = "cam_status") %>% 
  mutate(date = ymd(date)) %>% 
  unite("station_id", c(hex_id, station_num, year), sep = "_", remove = TRUE)

#########################################################################
##
## 4.  Construct daily bait matrix from 2 week interval data
##
##########################################################################

#prep data
bait_mat <- prep_matrix(dat, "bait_days_good")

bait_mat2 <- prep_matrix_part2(bait_mat)

#convert the dataframe into a daily actiity matrix
camop_bait <- camtrapR::cameraOperation(CTtable = bait_mat2,
                                   stationCol = "station_id",
                                   setupCol = "set_date",
                                   retrievalCol = "pull_date",
                                   writecsv = FALSE,
                                   hasProblems = TRUE,
                                   dateFormat = "ymd"
)


#fix negative and non-integer values using the above function. 
camop_bait_fixed <- camop_bait %>% as.data.frame() %>% 
  mutate(across(everything(), fix_camop)) %>% 
  as.matrix()

#rejoin columns dropped by camtrapR funtction
cols_to_join_bait <- bait_mat %>% select(station_id:station_num, effort_correction_factor, bait_days_good, snare_days_good)

bait_mat_final <- camop_bait_fixed  %>% as.data.frame() %>% rownames_to_column("station_id") %>% 
  left_join(cols_to_join_bait, by = join_by(station_id)) %>%
  select(station_id, hex_id, station_num, year, utm_e, utm_n, cell_id, effort_correction_factor, bait_days_good, snare_days_good, everything()) %>% 
  rename(effort_correction_factor_bait = effort_correction_factor)


#write_csv(bait_mat_final, "data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_bait_mat_10-11-23.csv")

bait_mat_long <- bait_mat_final  %>% 
  pivot_longer(cols= -c(station_id:snare_days_good), names_to = "date", values_to = "bait_status") %>% 
  mutate(date = ymd(date)) %>% 
  unite("station_id", c(hex_id, station_num, year), sep = "_", remove = TRUE) %>% 
  select(station_id, date, bait_status, effort_correction_factor_bait)

#########################################################################
##
## 5.  Construct daily snare matrix from 2 week interval data
##
##########################################################################

#prep data
snare_mat <- prep_matrix(dat, "snare_days_good")

snare_mat2 <- prep_matrix_part2(snare_mat)

#convert the dataframe into a daily actiity matrix
camop_snare <- camtrapR::cameraOperation(CTtable = snare_mat2,
                                        stationCol = "station_id",
                                        setupCol = "set_date",
                                        retrievalCol = "pull_date",
                                        writecsv = FALSE,
                                        hasProblems = TRUE,
                                        dateFormat = "ymd"
)


#fix negative and non-integer values using the above function. 
camop_snare_fixed <- camop_snare %>% as.data.frame() %>% 
  mutate(across(everything(), fix_camop)) %>% 
  as.matrix()

#rejoin columns dropped by camtrapR funtction
cols_to_join_snare <- snare_mat %>% select(station_id:station_num, effort_correction_factor, bait_days_good, snare_days_good)

snare_mat_final <- camop_snare_fixed  %>% as.data.frame() %>% rownames_to_column("station_id") %>% 
  left_join(cols_to_join_snare, by = join_by(station_id)) %>%
  select(station_id, hex_id, station_num, year, utm_e, utm_n, cell_id, effort_correction_factor, bait_days_good, snare_days_good, everything()) %>% 
  rename(effort_correction_factor_snare = effort_correction_factor)


#write_csv(snare_mat_final, "data/Camera_Data/Olympic_National_Park/cam_act/onp_fisher_2013_2016_snare_mat_10-11-23.csv")

snare_mat_long <- snare_mat_final  %>% 
  pivot_longer(cols= -c(station_id:snare_days_good), names_to = "date", values_to = "snare_status") %>% 
  mutate(date = ymd(date)) %>% 
  unite("station_id", c(hex_id, station_num, year), sep = "_", remove = TRUE) %>% 
  select(station_id, date, snare_status, effort_correction_factor_snare)



#########################################################################
##
## 6. Import and format species detections
##
##########################################################################

species_det <- read_csv("data/Camera_Data/Olympic_National_Park/raw/raw/onp_species_det_2013_2016.csv") %>% 
  clean_names()

species_det2 <- species_det %>% 
  rename(year = rep_year) %>% 
  unite("station_id", c(hex_id, station_num, year), sep = "_", remove = FALSE) %>% 
  select(station_id, start_date, species) %>% 
  mutate(start_date = mdy(start_date)) %>% 
  rename(date = start_date)

cougar_det <- species_det2 %>% filter(species == "Cougar")

# 
# cam_act2 <- cam_act %>%
#   unite("station_id", c(hex_id, station_num, year), sep = "_", remove = FALSE)
# 
# 
# #pivot longer and clean up
# cam_act_long <- cam_act2 %>%
#   pivot_longer(cols=c(-c(station_id:snare_days_good)), names_to= "date", values_to = "cam_status") %>%
#   clean_names() %>%
#   mutate(date = ymd(date))

#########################################################################
##
## 7. Generate count of independent cougar detections for each day at each camera station
##
##########################################################################

cougar_counts <- cougar_det %>% 
  merge_events() %>% 
  group_by(station_id, date) %>% 
  count() %>% 
  ungroup()

#########################################################################
##
## 8. Join cougar detection counts to activity history by station name and date
##
##########################################################################

#make long form detection history
#cam_status: NA = not deployed; 0 = deployed and inactive, 1 = deployed and active
#effort: 0 = undeployed/inactive, 1 = deployed and active, 2 = multiple cameras deployed and active

#long form det history
det_hist_full <- act_mat_long %>% 
  left_join(cougar_counts, join_by(station_id, date)) %>% 
  left_join(bait_mat_long, join_by(station_id, date)) %>% 
  left_join(snare_mat_long, join_by(station_id, date)) %>% 
  mutate(cougar_detections_count = case_when(is.na(n)==F ~ n,
                                             cam_status > 0 & is.na(n)==T ~ 0,
                                             is.na(cam_status)==T ~ NA)) %>%
  mutate(cougar_detections_binary = case_when(cougar_detections_count>=1 ~ 1,
                                              .default = cougar_detections_count)) %>% 
  select(station_id:cell_id, camera_days_good, effort_correction_factor, bait_days_good, effort_correction_factor_bait, snare_days_good, effort_correction_factor_snare, date, cam_status, bait_status, snare_status, cougar_detections_count, cougar_detections_binary)


det_hist_corrected <- det_hist_full %>% 
  mutate(cam_status = case_when(cougar_detections_binary > 0 ~ 1,
                                .default = cam_status)) %>% 
  mutate(grid_id = "ONP", .before = station_id)
 

#########################################################################
##
## 9. Create binary and count detection histories and export
##
##########################################################################

det_hist_binary <- det_hist_full %>% 
  select(-c(cougar_detections_count, effort, cam_status)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_binary) %>% 
  mutate(grid_id = "ONP", .before=station_id)


det_hist_counts <- det_hist_full %>% 
  select(-c(cougar_detections_binary, effort, cam_status)) %>% 
  pivot_wider(names_from = date, values_from = cougar_detections_count) %>% 
  mutate(grid_id = "ONP", .before=station_id)

effort_matrix <- det_hist_full %>% 
  select(-c(cougar_detections_binary, cougar_detections_count, cam_status)) %>% 
  pivot_wider(names_from = date, values_from = effort) %>% 
  mutate(grid_id = "ONP", .before=station_id)


# survey_periods2 <- det_hist_binary %>% 
#   pivot_longer(cols = -c(grid_id:snare_days_good), names_to = "date", values_to = "value") %>% 
#   #bind_rows() %>% 
#   group_by(year) %>% 
#   summarize(min_date = min(date), 
#             max_date = max(date))



# write_csv(det_hist_binary, "data/Camera_Data/Olympic_National_Park/onp_fisher_2013_2016_det_hist.csv")
# write_csv(det_hist_counts, "data/Camera_Data/Olympic_National_Park/onp_fisher_2013_2016_det_hist_counts.csv")
#write_csv(det_hist_corrected, "data/Camera_Data/Olympic_National_Park/onp_fisher_2013_2016_det_hist_long_10-11-23.csv")













################################ Graveyard #################################

#########################################################################
##
## 6. Extract raster cell number for each camera station
##
##########################################################################


# Load your sf object with points for station locations
# sf_points <- dat_format %>% st_as_sf(coords=c("utm_e", "utm_n"), crs = 26910, remove=FALSE) %>% 
#   st_transform(crs = 5070)
# 
# # Extract cell numbers for each station location
# xy <- sf_points %>% st_coordinates()
# 
# 
# #append raster cell id to each station location
# dat_format <- dat_format %>% mutate(cell_id = cellFromXY(temp, xy), .after = utm_n)
# 
# #check for multiple cameras in same raster cell
# get_dupes(dat_format, station_id, cell_id)
# 


#Function to reframe data by station
# format_station <- function(station_dat){
#   
#   round_1 <- station_dat %>% 
#     mutate(problem = case_when(interval != camera_days_good ~ TRUE,
#                                interval == camera_days_good ~ FALSE), .after = camera_days_good) %>% 
#     mutate(int_start = visit_date - days(interval), .before = visit_date) %>% 
#     mutate(problem_days = interval - camera_days_good, .after = camera_days_good) %>% 
#     rename(int_end = visit_date) %>% 
#     mutate(problem_start = case_when(problem == TRUE ~ int_start,
#                                      problem == FALSE ~ NA),
#            problem_end = case_when(problem == TRUE ~ int_end,
#                                    problem == FALSE ~ NA), .after = problem,
#            effort_correction_factor = case_when(problem == TRUE ~ camera_days_good,
#                                                 problem == FALSE ~ 0))
#   
#   if(nrow(station_dat) == 1){
#     
#     station_dat_formatted <- round_1 %>% 
#       mutate(visit_num = 1) %>% 
#       pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
#       reframe(.by = station_id,
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               problem_start_1 = min(problem_start_1, na.rm = T),
#               problem_end_1 = max(problem_end_1, na.rm = T),
#               problem_start_2 = NA,
#               problem_end_2 = NA,
#               problem_start_3 = NA,
#               problem_end_3 =NA,
#               problem_start_4 = NA,
#               problem_end_4 = NA,
#               bait_days_good = sum(bait_days_good, na.rm = T),
#               snare_days_good = sum(snare_days_good, na.rm = T),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = T),
#               photo_count = sum(photo_count, na.rm = T),
#               cougar = sum(cougar)
#       ) %>% 
#       mutate_all(~ replace(., is.infinite(.), NA))
#   }
#   
#   if(nrow(station_dat) == 2){
#     
#     station_dat_formatted <- round_1 %>% 
#       mutate(visit_num = 1:nrow(station_dat)) %>% 
#       pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
#       reframe(.by = station_id,
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               problem_start_1 = min(problem_start_1, na.rm = T),
#               problem_end_1 = max(problem_end_1, na.rm = T),
#               problem_start_2 = min(problem_start_2, na.rm = T),
#               problem_end_2 = max(problem_end_2, na.rm = T),
#               problem_start_3 = NA,
#               problem_end_3 =NA,
#               problem_start_4 = NA,
#               problem_end_4 = NA,
#               bait_days_good = sum(bait_days_good, na.rm = T),
#               snare_days_good = sum(snare_days_good, na.rm = T),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = T),
#               photo_count = sum(photo_count, na.rm = T),
#               cougar = sum(cougar)
#       ) %>% 
#       mutate_all(~ replace(., is.infinite(.), NA))
#   }
#   
#   if(nrow(station_dat) == 3){
#     
#     station_dat_formatted <- round_1 %>% 
#       mutate(visit_num = 1:nrow(station_dat)) %>% 
#       pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>%
#       
#       reframe(.by = station_id,
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               problem_start_1 = min(problem_start_1, na.rm = T),
#               problem_end_1 = max(problem_end_1, na.rm = T),
#               problem_start_2 = min(problem_start_2, na.rm = T),
#               problem_end_2 = max(problem_end_2, na.rm = T),
#               problem_start_3 = min(problem_start_3, na.rm = T),
#               problem_end_3 = max(problem_end_3, na.rm = T),
#               problem_start_4 = NA,
#               problem_end_4 = NA,
#               bait_days_good = sum(bait_days_good, na.rm = T),
#               snare_days_good = sum(snare_days_good, na.rm = T),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = T),
#               photo_count = sum(photo_count, na.rm = T),
#               cougar = sum(cougar)
#       ) %>% 
#       mutate_all(~ replace(., is.infinite(.), NA))
#   }
#   
#   if(nrow(station_dat) == 4){
#     
#     station_dat_formatted <- round_1 %>% 
#       mutate(visit_num = 1:nrow(station_dat)) %>% 
#       pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
#       reframe(.by = station_id,
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               problem_start_1 = min(problem_start_1, na.rm = T),
#               problem_end_1 = max(problem_end_1, na.rm = T),
#               problem_start_2 = min(problem_start_2, na.rm = T),
#               problem_end_2 = max(problem_end_2, na.rm = T),
#               problem_start_3 = min(problem_start_3, na.rm = T),
#               problem_end_3 = max(problem_end_3, na.rm = T),
#               problem_start_4 = min(problem_start_4, na.rm = T),
#               problem_end_4 = max(problem_end_4, na.rm = T),
#               bait_days_good = sum(bait_days_good, na.rm = T),
#               snare_days_good = sum(snare_days_good, na.rm = T),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = T),
#               photo_count = sum(photo_count, na.rm = T),
#               cougar = sum(cougar)
#       ) %>% 
#       mutate_all(~ replace(., is.infinite(.), NA))
#   }
#   
#   
#   return(station_dat_formatted)
# }
# 
# 
# 
# 
# 
# #function to format yearly data by station_id
# format_station_by_year <- function(year_dat){
#   round_1 <- year_dat %>%  
#     mutate(station_id = as.character(station_id)) %>% 
#     mutate(station_id = as.factor(station_id))
#   dat_split <- split(round_1, round_1$station_id)
#   formatted <- dat_split %>% map(format_station)
#   formatted <- bind_rows(formatted)
#   return(formatted)}

#function to format yearly data by station_id
# format_bait_by_year <- function(year_dat){
#   round_1 <- year_dat %>%  
#     mutate(station_id = as.character(station_id)) %>% 
#     mutate(station_id = as.factor(station_id))
#   dat_split <- split(round_1, round_1$station_id)
#   formatted <- dat_split %>% map(format_station_bait)
#   formatted <- bind_rows(formatted)
#   return(formatted)}


# format_station_bait <- function(station_dat){
#   
#   round_1 <- station_dat %>% 
#     mutate(problem = case_when(interval != bait_days_good ~ TRUE,
#                                interval == bait_days_good ~ FALSE), .after = bait_days_good) %>% 
#     mutate(int_start = visit_date - days(interval), .before = visit_date) %>% 
#     mutate(problem_days = interval - bait_days_good, .after = camera_days_good) %>% 
#     rename(int_end = visit_date) %>% 
#     mutate(problem_start = case_when(problem == TRUE ~ int_start,
#                                      problem == FALSE ~ NA),
#            problem_end = case_when(problem == TRUE ~ int_end,
#                                    problem == FALSE ~ NA), .after = problem,
#            effort_correction_factor = case_when(problem == TRUE ~ bait_days_good,
#                                                 problem == FALSE ~ 0))
#   
#   if(nrow(station_dat) == 1){
#     
#     station_dat_formatted <- round_1 %>% 
#       mutate(visit_num = 1) %>% 
#       pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
#       reframe(.by = station_id,
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               problem_start_1 = min(problem_start_1, na.rm = T),
#               problem_end_1 = max(problem_end_1, na.rm = T),
#               problem_start_2 = NA,
#               problem_end_2 = NA,
#               problem_start_3 = NA,
#               problem_end_3 =NA,
#               problem_start_4 = NA,
#               problem_end_4 = NA,
#               bait_days_good = sum(bait_days_good, na.rm = T),
#               snare_days_good = sum(snare_days_good, na.rm = T),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = T),
#               photo_count = sum(photo_count, na.rm = T),
#               cougar = sum(cougar)
#       ) %>% 
#       mutate_all(~ replace(., is.infinite(.), NA))
#   }
#   
#   if(nrow(station_dat) == 2){
#     
#     station_dat_formatted <- round_1 %>% 
#       mutate(visit_num = 1:nrow(station_dat)) %>% 
#       pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
#       reframe(.by = station_id,
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               problem_start_1 = min(problem_start_1, na.rm = T),
#               problem_end_1 = max(problem_end_1, na.rm = T),
#               problem_start_2 = min(problem_start_2, na.rm = T),
#               problem_end_2 = max(problem_end_2, na.rm = T),
#               problem_start_3 = NA,
#               problem_end_3 =NA,
#               problem_start_4 = NA,
#               problem_end_4 = NA,
#               bait_days_good = sum(bait_days_good, na.rm = T),
#               snare_days_good = sum(snare_days_good, na.rm = T),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = T),
#               photo_count = sum(photo_count, na.rm = T),
#               cougar = sum(cougar)
#       ) %>% 
#       mutate_all(~ replace(., is.infinite(.), NA))
#   }
#   
#   if(nrow(station_dat) == 3){
#     
#     station_dat_formatted <- round_1 %>% 
#       mutate(visit_num = 1:nrow(station_dat)) %>% 
#       pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>%
#       
#       reframe(.by = station_id,
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               problem_start_1 = min(problem_start_1, na.rm = T),
#               problem_end_1 = max(problem_end_1, na.rm = T),
#               problem_start_2 = min(problem_start_2, na.rm = T),
#               problem_end_2 = max(problem_end_2, na.rm = T),
#               problem_start_3 = min(problem_start_3, na.rm = T),
#               problem_end_3 = max(problem_end_3, na.rm = T),
#               problem_start_4 = NA,
#               problem_end_4 = NA,
#               bait_days_good = sum(bait_days_good, na.rm = T),
#               snare_days_good = sum(snare_days_good, na.rm = T),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = T),
#               photo_count = sum(photo_count, na.rm = T),
#               cougar = sum(cougar)
#       ) %>% 
#       mutate_all(~ replace(., is.infinite(.), NA))
#   }
#   
#   if(nrow(station_dat) == 4){
#     
#     station_dat_formatted <- round_1 %>% 
#       mutate(visit_num = 1:nrow(station_dat)) %>% 
#       pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
#       reframe(.by = station_id,
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               problem_start_1 = min(problem_start_1, na.rm = T),
#               problem_end_1 = max(problem_end_1, na.rm = T),
#               problem_start_2 = min(problem_start_2, na.rm = T),
#               problem_end_2 = max(problem_end_2, na.rm = T),
#               problem_start_3 = min(problem_start_3, na.rm = T),
#               problem_end_3 = max(problem_end_3, na.rm = T),
#               problem_start_4 = min(problem_start_4, na.rm = T),
#               problem_end_4 = max(problem_end_4, na.rm = T),
#               bait_days_good = sum(bait_days_good, na.rm = T),
#               snare_days_good = sum(snare_days_good, na.rm = T),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = T),
#               photo_count = sum(photo_count, na.rm = T),
#               cougar = sum(cougar)
#       ) %>% 
#       mutate_all(~ replace(., is.infinite(.), NA))
#   }
#   
#   
#   return(station_dat_formatted)
# }

#possible simpler version of formatting function

# format_station <- function(station_dat) {
#   library(dplyr)
#   
#   # Calculate the problem_start and problem_end columns
#   formatted_data <- station_dat %>%
#     mutate(problem = interval != camera_days_good,
#            int_start = visit_date - days(interval),
#            problem_days = interval - camera_days_good) %>%
#     rename(int_end = visit_date) %>%
#     group_by(station_id) %>%
#     summarize(across(starts_with("problem_start"), ~ min(., na.rm = TRUE), .names = "min_{.col}"),
#               across(starts_with("problem_end"), ~ max(., na.rm = TRUE), .names = "max_{.col}"),
#               year = first(year),
#               utm_e = first(utm_e),
#               utm_n = first(utm_n),
#               hex_id = first(hex_id),
#               station_num = first(station_num),
#               set_date = min(int_start),
#               pull_date = max(int_end),
#               problem_days = sum(problem_days),
#               problem_periods = sum(problem),
#               bait_days_good = sum(bait_days_good, na.rm = TRUE),
#               snare_days_good = sum(snare_days_good, na.rm = TRUE),
#               effort_correction_factor = sum(effort_correction_factor, na.rm = TRUE),
#               photo_count = sum(photo_count, na.rm = TRUE),
#               cougar = sum(cougar))
#   
#   # Replace infinite values with NA
#   formatted_data <- formatted_data %>% mutate_all(~ replace(., is.infinite(.), NA))
#   
#   return(formatted_data)
# }

#########################################################################
##
## 2. Create biweekly detection history
##
##########################################################################

#this function takes a vector of visit dates and calculates the interval between each. (first is 0)
# visit_int_f <- function(vec){
#   visit_int <- vector()
#   visit_int[1] = 0
#   if(length(vec) > 1){
#   for(i in 2:length(vec)){
#   visit_int[i] <- ymd(vec[i]) - ymd(vec[i-1])}
#     }
#   return(visit_int)
# }
# 
# 
# #split data
# dat_split <- split(dat, dat$station_id)
# 
# 
# test_dat <- dat_split[[1]]
# 
# map_fun <- function(x) {
#   x %>% 
#     arrange(visit_date) %>% 
#     mutate(visit_num = 1:nrow(x), .after = visit_date) %>% 
#     #select(-c(visit_date)) %>% 
#     pivot_wider(names_from = visit_num, values_from = c(visit_date, interval:cougar)) #%>% 
#     # mutate(int_start_1 = visit_date_1 - interval_1,
#     #        int_end_1 = visit_date_2,
#     #        int_start_2 = visit_date_2,
#     #        int_end_2 = visit_date_3, .after = interval_3) %>% View()
# 
# }
# 
# dat_format <- dat_split %>% 
#   map(map_fun)
# 
# dat_format <- bind_rows(dat_format) %>% 
#   select(station_id, hex_id, station_num, station:utm_n, cougar_1, cougar_2, cougar_3, cougar_4, everything())

#write_csv(dat_format, "data/Camera_Data/Olympic_National_Park/onp_fisher_2013_2016_cam_act_biweekly.csv")
# dat_format2 %>% 
#   pivot_longer(cols = Problem1_from:Problem4_to, names_to = "problem", values_to = "value")


# test_ <- dat_format2[969:972,] %>% select(-c(Problem4_from, Problem4_to)) %>% 
#   mutate(across(c(set_date, pull_date), as_datetime))
#subset one station in one year to test formatting protocol


# test2 <- test %>% 
#   mutate(problem = case_when(interval != camera_days_good ~ TRUE,
#                              interval == camera_days_good ~ FALSE), .after = camera_days_good) %>% 
#   mutate(int_start = visit_date - days(interval), .before = visit_date) %>% 
#   mutate(problem_days = interval - camera_days_good, .after = camera_days_good) %>% 
#   rename(int_end = visit_date) %>% 
#   mutate(problem_start = case_when(problem == TRUE ~ int_start,
#                                    problem == FALSE ~ NA),
#          problem_end = case_when(problem == TRUE ~ int_end,
#                                    problem == FALSE ~ NA), .after = problem) %>% 
#   pivot_wider(names_from = visit_num, values_from = visit_date) 



#unnecessary attempt to merge problem date columns

# # test_co <- dat_format2 %>%
# #   select(Problem1_from:Problem4_to) %>% 
# #   mutate(start = pmin(ymd(start1), ymd(start2), na.rm = TRUE),
# #          end = pmax(ymd(end1), ymd(end2), na.rm = TRUE)) %>%
# #   select(start, end)
# 
# #if problem1_from is NA, 
# 
# co_test_dat <- dat_format2[711,]
# 
# if(co_test_dat$Problem1_to == co_test_dat$Problem2_from & co_test_dat$Problem2_to == co_test_dat$Problem3_from){TRUE} else{FALSE}
# 
# test_co <- co_test_dat %>% 
#   mutate(Problem1_to = case_when(Problem1_to == Problem2_from &  is.na(Problem3_from) ==T ~  Problem2_to,
#                                   Problem1_to == Problem2_from & Problem2_to == Problem3_from & is.na(Problem4_from)==T ~ Problem3_to,
#                                   Problem1_to == Problem2_from & Problem2_to == Problem3_from & Problem3_to == Problem4_from ~ Problem4_to,
#                                   .default = Problem1_to) )
# 
# for(i in 1:nrow(dat)){
#   if()
#   
#   
# }
# 
# test_co <- co_test_dat %>% 
#   pivot_longer(cols = contains("Problem"), names_to = "problem", values_to = "period") %>% 
#   filter(!is.na(period)) %>% 
#   pivot_wider(names_from = problem, values_from = period) %>% 
#   mutate(across(contains("Problem"), ymd)) %>% 
#   mutate(period_1 = interval(Problem1_from, Problem1_to),
#          period_2 = interval(Problem2_from, Problem2_to),
#          period_3 = interval(Problem3_from, Problem3_to)) %>% 
#   mutate(across(contains("period"), reduce, .f = union))
#   
#   reduce(c(test_co$period_1, test_co$period_2, test_co$period_3), .f = union)
#   
#   pivot_longer(cols = contains("period"), names_to = "period", values_to = "range")
# 
# 
# 
# 
# test3 <- interval(start = ymd(co_test_dat$Problem1_to), end = ymd(co_test_dat$Problem1_from))




