#### Olympic National Park Detection History Formatting ####

# Author: Read Barbee

# Date:2023-09-11 
#Last updated: 2023-10-02

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)

#########################################################################
##
## 1. Import and format detection data 
##
##########################################################################

dat <- read_csv("data/Camera_Data/Olympic_National_Park/ONP_fisher_grid_2013-2016_raw_formatted.csv") %>% 
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


#########################################################################
##
## 2. Create biweekly detection history
##
##########################################################################

#this function takes a vector of visit dates and calculates the interval between each. (first is 0)
visit_int_f <- function(vec){
  visit_int <- vector()
  visit_int[1] = 0
  if(length(vec) > 1){
  for(i in 2:length(vec)){
  visit_int[i] <- ymd(vec[i]) - ymd(vec[i-1])}
    }
  return(visit_int)
}


#split data
dat_split <- split(dat, dat$station_id)


test_dat <- dat_split[[1]]

map_fun <- function(x) {
  x %>% 
    arrange(visit_date) %>% 
    mutate(visit_num = 1:nrow(x), .after = visit_date) %>% 
    #select(-c(visit_date)) %>% 
    pivot_wider(names_from = visit_num, values_from = c(visit_date, interval:cougar)) #%>% 
    # mutate(int_start_1 = visit_date_1 - interval_1,
    #        int_end_1 = visit_date_2,
    #        int_start_2 = visit_date_2,
    #        int_end_2 = visit_date_3, .after = interval_3) %>% View()

}

dat_format <- dat_split %>% 
  map(map_fun)

dat_format <- bind_rows(dat_format) %>% 
  select(station_id, hex_id, station_num, station:utm_n, cougar_1, cougar_2, cougar_3, cougar_4, everything())


#########################################################################
##
## 2. Extract raster cell number for each camera station
##
##########################################################################
# Load the terra package
library(terra)
library(sf)

# Load your raster data
raster <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_asp.tif")

temp <- raster[[1]]

# Load your sf object with points for station locations
sf_points <- dat_format %>% st_as_sf(coords=c("utm_e", "utm_n"), crs = 26910, remove=FALSE) %>% 
  st_transform(crs = 5070)

# Extract cell numbers for each station location
xy <- sf_points %>% st_coordinates()


#append raster cell id to each station location
dat_format <- dat_format %>% mutate(cell_id = cellFromXY(temp, xy), .after = utm_n)

#check for multiple cameras in same raster cell
get_dupes(dat_format, station_id, cell_id)


#write_csv(dat_format, "data/Camera_Data/Olympic_National_Park/onp_fisher_2013_2016_cam_act_biweekly.csv")

#########################################################################
##
## 3. Construct daily camera activity matrix from 2 week interval data
##
##########################################################################

#Function to reframe data by station
format_station <- function(station_dat){

round_1 <- station_dat %>% 
  mutate(problem = case_when(interval != camera_days_good ~ TRUE,
                                    interval == camera_days_good ~ FALSE), .after = camera_days_good) %>% 
  mutate(int_start = visit_date - days(interval), .before = visit_date) %>% 
  mutate(problem_days = interval - camera_days_good, .after = camera_days_good) %>% 
  rename(int_end = visit_date) %>% 
  mutate(problem_start = case_when(problem == TRUE ~ int_start,
                                   problem == FALSE ~ NA),
         problem_end = case_when(problem == TRUE ~ int_end,
                                 problem == FALSE ~ NA), .after = problem,
         effort_correction_factor = case_when(problem == TRUE ~ camera_days_good,
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

#split data into list elements by station_id to iterate over 
dat_split <- split(dat, dat$station_id)

#apply the formatting function to each station
dat_format <- dat_split %>% map(format_station)

#bind the formatted list elements into a single dataframe
dat_format <- bind_rows(dat_format)

#format dataframe for captrapR::cameraOperation function 
dat_format2 <- dat_format %>% 
  rename(Problem1_from = problem_start_1,
         Problem1_to = problem_end_1,
         Problem2_from = problem_start_2,
         Problem2_to = problem_end_2,
         Problem3_from = problem_start_3,
         Problem3_to = problem_end_3,
         Problem4_from = problem_start_4,
         Problem4_to = problem_end_4,) %>% 
  select(station_id, set_date, pull_date, Problem1_from:Problem4_to) %>% 
  mutate(across(c(set_date:Problem4_to), as.character))


#convert the dataframe into a daily actiity matrix
camop <- camtrapR::cameraOperation(CTtable = dat_format2,
                         stationCol = "station_id",
                         setupCol = "set_date",
                         retrievalCol = "pull_date",
                         writecsv = FALSE,
                         hasProblems = TRUE,
                         dateFormat = "ymd"
)

#function to convert all values to be 0, 1, or NA
fix_camop <- function(col){
  col <- case_when(
    col < 0 ~ 0,
    col > 1 ~ 1,
    col > 0 & col < 1 ~ 1,
    .default = col
  )
}

#fix negative and non-integer values using the above function. 
camop_fixed <- camop %>% as.data.frame() %>% 
  mutate(across(everything(), fix_camop)) %>% 
  as.matrix()
  

#plot camera operation matrix
camtrapR:::camopPlot(camOp = camop_fixed, palette = "Heat", lattice = TRUE)


#########################################################################
##
## 3. Format activity matrix and export
##
##########################################################################

cols_to_join <- dat_format %>% select(station_id:station_num, effort_correction_factor, bait_days_good, snare_days_good)

cam_act_final <- camop_fixed  %>% as.data.frame() %>% rownames_to_column("station_id") %>% 
  left_join(cols_to_join, by = join_by(station_id)) %>%
  select(station_id, hex_id, station_num, year, utm_e, utm_n, effort_correction_factor, bait_days_good, snare_days_good, everything())


#write_csv(cam_act_final, "data/Camera_Data/Olympic_National_Park/onp_fisher_2013_2016_cam_act.csv")





################################ Graveyard #################################
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




