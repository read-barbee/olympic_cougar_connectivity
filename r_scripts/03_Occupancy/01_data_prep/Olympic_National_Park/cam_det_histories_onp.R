#### Olympic National Park Detection History Formatting ####

# Author: Read Barbee

# Date:2023-09-11 

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(janitor)

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
  mutate(station_id = as.factor(station_id))


#creates columns of binary cougar detection for each visit interval (~14 days each)
dat_wide<- dat %>% 
  mutate(station = paste0(hex_id, "_", station_num), .after=station_id) %>% 
  #group_by(station_id) %>%
  pivot_wider(names_from = visit_num, values_from =visit_date:cougar) %>% 
  select(station_id, station, year, cougar_1:cougar_5, everything())


#subset one station in one year to test formatting protocol
test <- dat %>% 
  filter(station_id=="109_r1_s1.1_2016")



#########################################################################
##
## 2. Attempt to construct daily activity matrix
##
##########################################################################

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
            photo_count = sum(photo_count, na.rm = T),
            cougar = sum(cougar)
    ) %>% 
    mutate_all(~ replace(., is.infinite(.), NA))
}


  return(station_dat_formatted)
}


dat_split <- split(dat, dat$station_id)


dat_format <- dat_split %>% map(format_station)

dat_format <- bind_rows(dat_format)

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


dat_format2 %>% 
  pivot_longer(cols = Problem1_from:Problem4_to, names_to = "problem", values_to = "value")


camop <- camtrapR::cameraOperation(CTtable = dat_format2,
                         stationCol = "station_id",
                         setupCol = "set_date",
                         retrievalCol = "pull_date",
                         writecsv = FALSE,
                         hasProblems = TRUE,
                         dateFormat = "ymd"
)



#########################################################################
##
## 2. Extract covariate values for each site
##
##########################################################################

# Load the terra package
library(terra)
library(sf)

# Load your raster data
cov_stack <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

names(cov_stack) <- c("tree_cover_hansen",
                      "gpp",
                      "infra_hii",
                      "landuse_hii",
                      "land_cover_usfs",
                      "land_use_usfs",
                      "npp",
                      "popdens_hii",
                      "power_hii",
                      "precip",
                      "rails_hii",
                      "roads_hii",
                      "elevation",
                      "slope",
                      "aspect",
                      "tri",
                      "tpi",
                      "perc_tree_cover",
                      "perc_nontree_veg",
                      "perc_nonveg",
                      "ndvi",
                      "evi",
                      "dist_water")


# Load your sf object with points for station locations. native crs is NAD83 / UTM zone 10N
sf_points <- dat_wide %>% st_as_sf(coords=c("utm_e", "utm_n"), crs = 26910, remove=FALSE) %>% 
  st_transform(crs=5070)

# Extract cell numbers for each station location
covs <-terra::extract(cov_stack, sf_points)

occ_dat <- dat_wide %>% 
  bind_cols(covs) %>%  
  select(-ID) %>% 
  mutate(station_id = str_remove(station_id, "_r0")) %>% 
  mutate(station_id = str_remove(station_id, "_r1")) %>% 
  mutate(station_id = str_remove(station_id, "_r2")) %>% 
  mutate(across(visit_date_1:visit_date_5, mdy)) %>% 
  mutate(across(visit_date_1:visit_date_5, yday))


#Create stacked unmarked frame
library(unmarked)

umf <- unmarkedFrameOccu(y = occ_dat[,4:8],
                         siteCovs = occ_dat %>% select(station_id, tree_cover_hansen:dist_water),
                         obsCovs = list(date = occ_dat %>% select(visit_date_1:visit_date_5),
                                        int = occ_dat %>% select(interval_1:interval_5),
                                        effort = occ_dat %>% select(camera_days_good_1:camera_days_good_5),
                                        bait_tree = occ_dat %>% select(bait_days_good_1:bait_days_good_5),
                                        bait_snare = occ_dat %>% select(snare_days_good_1:snare_days_good_5)))






library(ubms)

fit_stack <- stan_occu(~ scale(date) + scale(int) + scale(effort) + scale(bait_tree) + scale(bait_snare)  ~scale(tree_cover_hansen) + (1|station_id), 
                       data=umf, chains=3, iter=100)
fit_stack


################################ Graveyard #################################

#format by hex ID
# format_hex <- function(hex_dat){
# 
#   hex_dat_formatted <- hex_dat %>% 
#     mutate(problem = case_when(interval != camera_days_good ~ TRUE,
#                                interval == camera_days_good ~ FALSE), .after = camera_days_good) %>% 
#     mutate(int_start = visit_date - days(interval), .before = visit_date) %>% 
#     mutate(problem_days = interval - camera_days_good, .after = camera_days_good) %>% 
#     rename(int_end = visit_date) %>% 
#     mutate(problem_start = case_when(problem == TRUE ~ int_start,
#                                      problem == FALSE ~ NA),
#            problem_end = case_when(problem == TRUE ~ int_end,
#                                    problem == FALSE ~ NA), .after = problem) %>% 
#     pivot_wider(names_from = visit_num, values_from = c(problem_start, problem_end)) %>% 
#     reframe(.by = station_id,
#             year = first(year),
#             utm_e = first(utm_e),
#             utm_n = first(utm_n),
#             hex_id = first(hex_id),
#             station_num = first(station_num),
#             set_date = min(int_start),
#             pull_date = max(int_end),
#             problem_days = sum(problem_days),
#             problem_periods = sum(problem),
#             problem_start_1 = min(problem_start_1, na.rm = T),
#             problem_end_1 = max(problem_end_1, na.rm = T),
#             problem_start_2 = min(problem_start_2, na.rm = T),
#             problem_end_2 = max(problem_end_2, na.rm = T),
#             problem_start_3 = min(problem_start_3, na.rm = T),
#             problem_end_3 = max(problem_end_3, na.rm = T),
#             bait_days_good = sum(bait_days_good, na.rm = T),
#             snare_days_good = sum(snare_days_good, na.rm = T),
#             photo_count = sum(photo_count, na.rm = T),
#             cougar = sum(cougar)
#     ) %>% 
#     mutate_all(~ replace(., is.infinite(.), NA)) %>% View()
#   
#   return(hex_dat_formatted)
# }


