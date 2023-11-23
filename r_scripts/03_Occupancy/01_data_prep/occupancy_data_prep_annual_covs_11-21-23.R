#### Univariate occupancy analysis ####

# Author: Read Barbee

# Date:2023-10-12 

#Last updated: 2023-11-21

# Purpose:



library(tidyverse)
library(sf)
library(terra)
library(ubms)
library(beepr)
#library(doParallel)

################################ User-defined Parameters #################################

#  Sampling occasion specs
survey_period <- "day" 
number_of_periods <- 14 # how many days you want each survey to be

cov_folder_path <- "/Users/tb201494/Desktop/annual_cov_stacks_1km_buffer"

#if there are 362 distinct julian days represented in data, why are there only 331 in w_obs?
#there are the same number of unique intervals generated as unique dates 
#but any given cell_year doesn't have more than 331 days in its matrix

##test <- dat %>% group_by(cell_id, year) %>% distinct(date) %>% count()
#range(test$n)

#########################################################################
##
##  1. Import stacked occupancy data
##
##########################################################################
#  Source functions to make occupancy surveys of different lenghts
source("r_scripts/03_Occupancy/01_data_prep/Utility/sampling_interval_aggregation/original/make_int_fun.R")

#example data
#load("r_scripts/03_Occupancy/01_data_prep/Utility/sampling_interval_aggregation/original/dat_example.Rdata") # object called dat

make_interval2 <- function(x, dt_col, time_int = "month", increment = 1) {
  stopifnot(tibble::is_tibble(x))
  tmp <- dplyr::pull(x, {{dt_col}})
  if (time_int == "hour") {
    stopifnot(inherits(tmp, c("POSIXct", "POSIXlt")))
  }
  stopifnot(all(!is.na(tmp)))
  stopifnot(time_int %in% c("hour", "day", "week", "month"))
  if (dplyr::is_grouped_df(x)) {
    args_passed <- as.list(match.call())
    stop(
      "Data is grouped, try something like",
      "\n\t",
      paste0(
        "x %>% ",
        "dplyr::do(gbn_make_interval(",
        args_passed$dt_col, ", ",
        "'", args_passed$time_int, "', ",
        args_passed$increment, "))"
      )
    )
  }
  
  int <- lubridate::interval(min(tmp), tmp)
  
  new <- x %>%
    dplyr::mutate(
      rnd_dt = lubridate::round_date(tmp, unit = paste(increment, time_int)),
      interval = switch(
        time_int,
        "month" = int %/% lubridate::period(increment, units = "month"),
        "week" = int %/% lubridate::weeks(increment),
        "day" = int %/% lubridate::days(increment),
        "hour" = int %/% lubridate::hours(increment)
      ) + 1
    ) 
  return(new)
}

#activity/detection data
dat <- read_csv("data/Camera_Data/master/ocp_onp_occ_dat_long_11-21-23.csv")



#########################################################################
##
##  2. Import covariate stacks
##
##########################################################################


cov_files <- list.files(cov_folder_path, pattern = "\\.tif$", full.names = TRUE)

# Create an empty list to store the raster objects
cov_stacks <- list()

# Loop through each .tif file, import it, and assign the file name as the object name
for (tif_file in cov_files) {
  raster_name <- tools::file_path_sans_ext(basename(tif_file))
  raster <- rast(tif_file)
  cov_stacks[[raster_name]] <- raster
}

#rename layers missing "annual"
for (i in 1:length(cov_stacks)){
  stack_name <- names(cov_stacks)[i]
  year <- substr(stack_name, nchar(stack_name) - 3, nchar(stack_name))
  layer_names_old <- names(cov_stacks[[i]])
  layer_names_new <- layer_names_old
  layer_names_new[5] <- paste0("perc_tree_cov_annual_", year)
  layer_names_new[6] <- paste0("perc_nontree_veg_annual_", year)
  layer_names_new[7] <- paste0("perc_nonveg_annual_", year)
  
names(cov_stacks[[i]]) <- layer_names_new
}

#import static stack
static_stack <- rast("/Users/tb201494/Desktop/1km_buffer/static_stack_1km_buffer_single_tpi.tif")
mtpi <- rast("/Users/tb201494/Desktop/1km_buffer/static/mtpi_1km_buffer.tif")

static_stack$mtpi <- mtpi

static_stack <- static_stack[[c("elevation", "slope", "aspect", "aspect_northness", "aspect_eastness", "tpi", "mtpi", "tri", "distance_water")]]


#########################################################################
##
##  3. Format data for new survey interval
##
##########################################################################

# Note date column should be date format:
# $ date   : Date, format: "2009-05-02" "2009-08-10" "2009-05-14" ...

#check for missing days
missing_days <- which(!(1:365 %in% dat$j_day ))

missing_dates <- list()
for(i in 1:length(missing_days)){
missing_dates[[i]] <- ymd("2020-01-01") + days(missing_days[i])
}

#missing 4 days march 23-26

# dat %>% group_by(grid_id, year) %>% 
#   summarize(start = min(date),
#             end = max(date)) %>% View()


#Assign each row to an observation interval
obs_tmp1 <- make_interval2(dat, date, survey_period, number_of_periods) 

#add column for survey interval
obs_tmp2 <- obs_tmp1 %>% 
	dplyr::group_by(cell_id, year) %>% 
	dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
	dplyr::ungroup()


#calculate the maximum number of observations for each survey interval
obs_tmp3 <- obs_tmp2 %>%
	group_by(survey_interval, cell_id, year) %>% #station_id, 
  summarize(station_id = first(station_id),
            grid_id = first(grid_id),
            lon = first(lon),
            lat = first(lat),
            utm_e = first(utm_e),
            utm_n = first(utm_n),
            # annual_effort_correction_factor = first(effort_correction_factor),
            # annual_effort_correction_factor_bait = first(effort_correction_factor_bait),
            # annual_effort_correction_factor_snare = first(effort_correction_factor_snare),
            cam_days = sum(cam_status, na.rm = TRUE),
            bait_days = sum(bait_status, na.rm = TRUE),
            snare_days = sum(snare_status, na.rm = TRUE),
            cougar_detections_binary = case_when(all(is.na(cougar_detections_binary))==T ~NA,
                                                 sum(cougar_detections_binary, na.rm = TRUE) == 0 ~ 0,
                                                 sum(cougar_detections_binary, na.rm = TRUE) > 0 ~ 1))



w_obs <- obs_tmp3 %>%
	arrange(survey_interval, cell_id, year) %>%
	pivot_wider(values_from = c(cam_days, bait_days, snare_days, cougar_detections_binary), names_from = survey_interval) %>% 
 mutate(across(!contains("detections") , \(x) replace_na(x, 0))) %>% 
	as.data.frame()
	

#########################################################################
##
##  4. Extract annual covariate values by year
##
##########################################################################

#make spatial points data frame
sf_pts <- w_obs %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(crs = 5070)

#mapview::mapview(sf_pts)

#static covariates
sf_pts2 <- terra::extract(static_stack, sf_pts, bind = TRUE) 

#Function to extract annual covariates
extract_annual_covs <- function(df, cov_stack_list){
  #create year column
  df2 <- df %>%
    mutate(year = as.factor(year), .before = station_id)
  
  #split data by year
  df_split <- split(df2, df2$year)

  #assign names to split data
  names_df_split <- vector()
 for(i in 1:length(df_split)){
    names_df_split[i] <- as.character(df_split[[i]]$year[1])
 }
  names(df_split) <- names_df_split
  
  #define function to select the raster stack for the relevant year to extract from and extract the covariates
  extract_fun <- function(df_split, year){
    stack <- cov_stack_list[[paste0("cov_stack_", as.character(year))]]
    covs <- terra::extract(stack, df_split, bind = TRUE)
    
    return(covs)
  }
  
  #apply the extraction function to each year
  df_split_covs <- list()
  for(i in 1:length(df_split)){
    names <- names(df_split)
    df_split_covs[[i]] <- extract_fun(df_split[[i]], names[i])
  }
  
  names(df_split_covs) <- names_df_split <- names
  
  remove_year_names <- function(df_split_covs){
    old_names <- names(df_split_covs)
    new_names <- vector()
    for(i in 1:length(old_names)){
      if(str_detect(old_names[i], coll("annual_20"))){
        new_names[i] <- substr(old_names[i], 1, nchar(old_names[i]) - 5)
      } else{
        new_names[i] <- old_names[i]
      }
    }
    
    names(df_split_covs) <- new_names
    
    return(df_split_covs)
  }
  
  #apply renaming function
  covs_renamed <- map(df_split_covs, remove_year_names)
  
  #bind all years together 
  df_covs_final <- covs_renamed %>% 
    map(as.data.frame) %>% 
    bind_rows() %>% 
    select(-f)
  
  return(df_covs_final)
  
}


#apply annual extraction function
occ_dat_all_covs <- extract_annual_covs(sf_pts2, cov_stacks)

#write_csv(occ_dat_all_covs, "data/Camera_Data/master/ocp_onp_occ_dat_annual_covs_11-21-23.csv")


