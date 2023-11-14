#### OCP iSSF Module_02: Prepare Screened GPS Data in amt ####

# Author: Read Barbee

# Date:2023-07-12 
#Last updated: 2023-11-13

# Purpose: Prepare screened GPS data in amt to fit iSSFs.

# Inputs:
#   •	Screened GPS Locations (corrected for error and habitat bias)
#   •	Deployment list
#   •	Disperser list
#   •	Covariate Stack
#
# Outputs:
#   •	amt_step dataframe with fields for all necessary covariates


################################ Libraries #################################
library(tidyverse)
library(terra)
library(amt)
library(janitor)
library(doParallel)
library(sf)
library(mapview)



################################ User-Defined Parameters #################################

project_crs <- "epsg:5070" #NAD83 Albers Equal Area Projection

resample_int_hours <- 2 #also try 3 and 6

tolerance_mins <- 15

rand_steps <- 10

cov_folder_path <- "data/Habitat_Covariates/annual_covariates/annual_stacks"

#plan(multisession, workers = 8)

#########################################################################
##
## 1. Import and format location data
##
##########################################################################
# Mountain lion location data 
locs_screened <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_10-02-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = with_tz(date_time_local, tzone = "US/Pacific"))

get_dupes(locs_screened, deployment_id, date_time_utc)

# # Mountain lion deployments 
# deployments <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Deployments/collar_deployments_master_7-11-2023.csv")
# 
# # Dispersal inventory from Teams 
# dispersals <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_10-02-2023.csv")


#########################################################################
##
## 2. Import and format covariate data
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

#import static stack
static_stack <- rast("data/Habitat_Covariates/annual_covariates/static_stack.tif")

#date_ranges
# years_ndvi_evi_precip <- 2010:2023
# years_npp_gpp_landcover <- 2010:2022
# years_veg_hii <- 2010:2020
# years_osm <- 2014:2023

# #add time stamps to all covariate layers for extraction
# cov_stacks_ts <- list()
# for(i in 1:length(cov_stacks)){
#   if(str_starts(names(cov_stacks)[i], coll("dist"))){
#     start_year <- 2014
#   } else{
#     start_year <- 2010
#   }
#   
#   years <- start_year:(start_year + length(names(cov_stacks[[i]]))-1)
#   
#   timestamps <- vector()
#   for(j in 1:length(names(cov_stacks[[i]]))){
#     timestamps[j] <- paste0(years[j], "-01-01 00:00:00")
#   }
#   
#   terra::time(cov_stacks[[i]]) <- ymd_hms(timestamps)
#   
#   cov_stacks_ts[[i]] <- cov_stacks[[i]]
# }
# 
# names(cov_stacks_ts) <- names(cov_stacks)

#########################################################################
##
## 3. Convert data to amt tracks
##
##########################################################################

#Create track for each animal

multi_track <- function(d){
  make_track(d, lon_utm, lat_utm, date_time_utc, check_duplicates = TRUE, all_cols = TRUE, crs = 5070) 
  #%>% transform_coords(crs_to = 5070)
}

#Nest locations by animal_id 
locs_nested <- locs_screened %>% 
  nest_by(animal_id, sex, dispersal_status) 


#create tracks
locs_nested$tracks <-map(locs_nested$data, multi_track)


#subset only disperser locations
# dispersers <- locs_nested %>% 
#   filter(dispersal_status=="disperser")

#########################################################################
##
## 4. Check for duplicate locations in tracks
##
##########################################################################
#check for duplicate points in tracks
bind_rows(locs_nested$tracks) %>% get_dupes(deployment_id, t_)


#########################################################################
##
## 5. Examine Sampling Rates
##
##########################################################################

#inspect sampling rate for each individual

locs_nested$sr <-  map(locs_nested$tracks, summarize_sampling_rate)

locs_nested$median_sr <-  map(locs_nested$tracks, function(x){
  summ <- summarize_sampling_rate(x)
  med <- round(summ$median)
  return(med)
})

#sampling rate statistics by individual
sampling_rates <- locs_nested %>% 
  select(-c(data, tracks, median_sr)) %>% 
  unnest(cols = c(sr)) 

summary(sampling_rates)


#########################################################################
##
## 6. Identify dispersal events (only need to do this once and save to csv)
##
##########################################################################
# 
# # Map track
# filt <- dispersers %>% 
#   select(data) %>% 
#   unnest(cols=data) %>% 
#   filter(animal_id=="Archie") %>% 
#   #filter(date_time_utc>= ymd("2020/03/15")) %>% 
#   sf::st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070) %>% ungroup()
# 
# mapview(filt)
# 
# #calculate net squared displacement from first location over time
# dispersers$tracks <- map(dispersers$tracks, add_nsd)
# 
# #plot nsd over time to identify dispersal date 
# nsd_plot <- ggplot(data = dispersers %>% 
#                  filter(animal_id=="Archie") %>% 
#                  select(tracks) %>% 
#                  unnest(cols=tracks), 
#                #%>% filter(t_>= ymd("2020/03/15")), 
#                aes(x = t_, y = nsd_)) +
#   geom_point()
# 
# 
# plotly::ggplotly(nsd_plot)
# 
# #recheck Hana, Lady, Lolli
# #Lolli's track not trimmed correctly. Also doesn't seem like much of a dispersal
# #Harder to determine dispersal date for females because they don't disperse as far
# 

#########################################################################
##
## 7. Identify optimal resampling interval
##
##########################################################################

#function to calculate bursts for each individual---increment by hours to find optimum
test_fun1 <- function(x) {
  x %>% 
    amt::track_resample(rate = hours(2), tolerance = minutes(15)) %>% #resample to 2 hours
    amt::filter_min_n_burst() #divide into bursts of min 3 pts
}

test <- locs_nested

test$steps <- map(locs_nested$tracks, test_fun1)


test %>% 
  filter(nrow(steps)<10) %>% #some individuals with low numbers of steps cause the amt functions to fail
  select(steps)

#make list of individuals to remove (because they can't be resampled at the specified interval)
indiv_to_remove <- test %>% 
  filter(nrow(steps)<10) %>% 
  pull(animal_id)


#1 hour removes 67 individuals (72 with imp points)
#2 hour removes 35 individuals (34 with imp points)
#3 hours removes 29 individuals (35 with imp points)
#4 hours removes 32 individuals (31 with imp points)
#5 hours removes 64 individuals (62 with imp points)
#6 hours removes 10 individuals (8 with imp points)

#6 hours retains the most individuals
# 2 hours may offer the best balance between data resolution and data loss


#########################################################################
##
## 8. Resample Tracks/Generate Steps
##
#######################################################################

#boundary polygon for random steps
op_poly <- st_read("data/Habitat_Covariates/study_area_polys/ocp_study_area_poly_wa_only_10-30-23.shp")

#op_poly_buffer <- st_buffer(op_poly, 200)

## Create non-habitat zones where random steps shouldnt be generated (i.e. water bodies)
#Import water body polygons
water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/OP_water_polys_v1.shp") %>% st_transform(crs = 5070)

water_polys_cropped <-  water_polys %>%
  filter(WB_PERIOD_ =="PER" ) %>% #only include permanent water bodies with area > 100 m2
  filter(SHAPEAREA >= 100) %>%
  sf::st_crop(op_poly) #crop to study area

water_polys_mask <- water_polys_cropped %>%
  sf::st_union() %>% #dissolve polygons into single vector mask layer
  st_sf()



#function to convert locatons to steps--turned off water intersection filter
steps_calc <- function(x) {
  part1 <- x %>% 
    amt::track_resample(rate = hours(resample_int_hours), tolerance = minutes(tolerance_mins)) %>% #resample 
    amt::filter_min_n_burst() %>% #divide into bursts of min 3 pts
    amt::steps_by_burst() %>%  #convert to steps
    dplyr::filter(sl_ >=100) %>% 
    dplyr::filter(sl_<=20000) #remove stationary and unreasonable step lengths
  blank <- list()
  
  if(nrow(part1) < 3){ return(blank)}
  
  else{
  
  #generate 10 times the amount of desired available steps
  part2 <- part1 %>% 
    amt::random_steps(n_control = rand_steps*10) %>% #generate random steps per used step
    #amt::extract_covariates(cov_stack, where = "end") %>% #extract covariates at end of step
    amt::time_of_day(include.crepuscule = FALSE) %>% #calculate time of day for each step--not working
    mutate(unique_step = paste(burst_,step_id_,sep="_"), .before=step_id_) %>%  #add column for unique step id }
    mutate(log_sl_ = log(sl_),    
           cos_ta_ = cos(ta_))
  
  #determine whether each random step intersects a water body and/or is within the study area
  part3 <- part2 %>% 
    st_as_sf(coords = c("x2_", "y2_"), crs = 5070, remove=FALSE) %>% 
    mutate(intersects_water = lengths(st_intersects(., water_polys_mask)) > 0) %>% 
    mutate(intersects_study_area = lengths(st_intersects(., op_poly)) > 0) %>% 
    dplyr::select(unique_step, t1_, t2_, x2_, y2_, intersects_water, intersects_study_area) %>% 
    as.data.frame() %>% 
    select(-geometry)
  
  #join intersection status to original step frame
  part4 <- part2 %>% 
    left_join(part3, by=join_by(unique_step, t1_, t2_, x2_, y2_))
  
  #split used and random steps to filter the random steps by intersection
  tf_split <- split(part4, part4$case_)
  
  #Keep only random steps within the study area that don't intersect water and randomly sample the desired number from those remaining
  tf_split$`FALSE` <- tf_split$`FALSE` %>% 
    filter(intersects_study_area == TRUE) %>% 
    #filter(intersects_water == FALSE) %>% 
    slice_sample(n = rand_steps, by = unique_step)
  
  #put the used and random steps back together in a single frame
  final <- bind_rows(tf_split) %>% arrange(t2_, desc(case_))
  
  return(final)
  }
}

#remove individuals not compatible with the sampling interval
amt_locs <- locs_nested %>% 
  filter(!(animal_id %in% indiv_to_remove)) %>% 
  filter(!(animal_id %in% c("Butch", "Crash", "Freya", "Hera", "Otook_Tom", "Promise")))
  # %>% filter(!(animal_id %in% c("Belle", "Cato")))


#map the step function over each individual and append as a nested dataframe
amt_locs$steps <- map(amt_locs$tracks, steps_calc)

#map resulting points and check for errors
# bind_rows(amt_locs$steps) %>% 
#   #filter(intersects_water==TRUE) %>% 
#   #filter(case_==TRUE) %>% 
#   st_as_sf(coords=c("x2_", "y2_"), crs = 5070) %>% 
#   mapview::mapview()

#select the relevant columns and unnest the steps column
amt_steps <- amt_locs %>% 
  select(animal_id:dispersal_status,
         steps) %>% 
  filter(length(steps)>0) #%>%  # removes any individuals without enough steps to fit distributions (should only be 1)
  #unnest(cols=c(steps)) %>% 
  #relocate(unique_step, .after = step_id_) %>% 
  #ungroup()
  
  # Check global step length and turn angle distributions
  amt_steps %>%
    unnest(cols=c(steps)) %>% 
    #filter(case_==TRUE) %>%
    filter(sl_ >100 & sl_<20000 ) %>%
    pull(sl_) %>%
    hist()
  
  amt_steps %>% 
    unnest(cols=c(steps)) %>%  
    pull(ta_) %>% 
    hist()
  
  #########################################################################
  ##
  ## 11. Extract covariate values
  ##
  #######################################################################
  
  amt_steps_all_covs <- amt_steps
  
  #static covariates
  amt_steps_all_covs$steps <- map(amt_steps_all_covs$steps, extract_covariates, covariates = static_stack)

  #annual covariates (make sure to adjust column indices)
  extract_annual_covs <- function(steps, cov_stack_list){
    #create year column
    steps2 <- steps %>%
      mutate(year = as.factor(year(t1_)), .before = t1_)
    
    #split data by year
    steps_split <- split(steps2, steps2$year)
    
    #define function to select the raster stack for the relevant year to extract from and extract the covariates
    extract_fun <- function(steps_split, year){
      stack <- cov_stack_list[[paste0("cov_stack_", as.character(year))]]
      covs <- steps_split %>% extract_covariates(stack)
      
      return(covs)
    }
    
    #apply the extraction function to each year
    steps_split_covs <- list()
    for(i in 1:length(steps_split)){
      names <- names(steps_split)
      steps_split_covs[[i]] <- extract_fun(steps_split[[i]], names[i])
    }
    
    names(steps_split_covs) <- names
    
    #define function to remove the year suffixes from covariate columns 
    # remove_year_names <- function(steps_split_covs, start_name, end_name){
    #   names_to_change <-  names(steps_split_covs)[start_name:end_name]
    #   names(steps_split_covs)[start_name:end_name]<- substr(names_to_change, 1, nchar(names_to_change) - 5)
    #   return(steps_split_covs)
    #   
    # }
    
    remove_year_names <- function(steps_split_covs){
      old_names <- names(steps_split_covs)
      new_names <- vector()
      for(i in 1:length(old_names)){
        if(str_detect(old_names[i], coll("20"))){
          new_names[i] <- substr(old_names[i], 1, nchar(old_names[i]) - 5)
        } else{
          new_names[i] <- old_names[i]
        }
      }
      
      names(steps_split_covs) <- new_names
      
      return(steps_split_covs)
    }
    
    #apply renaming function
    covs_renamed <- map(steps_split_covs, remove_year_names)
    
    #bind all years together 
    steps_covs_final <- bind_rows(covs_renamed)
    
    return(steps_covs_final)
    
  }
  

  amt_steps_all_covs$steps <- map(amt_steps_all_covs$steps, extract_annual_covs, cov_stack_list = cov_stacks)
  
  
  #unnest to view missing values
  #test <- amt_steps_all_covs %>% unnest(cols=c(steps))
  
  
#########################################################################
##
## 9. Extract covariates
##
##########################################################################
  
  #REDUNDANT
  
  #amt_steps$steps <- map(amt_steps$steps, extract_covariates, covariates = static_stack)
  
  # test_steps <- amt_steps$steps[[15]]
  # test_steps2 <- test_steps
  # 
  # extract_annual_covs <- function(steps){
  #   
  #   stack_names <- names(cov_stacks_ts)
  #   
  #   for(i in 1:length(cov_stacks_ts)){
  #     steps <- extract_covariates_var_time(steps,  covariates = cov_stacks_ts[[i]],
  #                                          when = "any", 
  #                                          where = "end",
  #                                          max_time = years(4),
  #                                          name_covar = stack_names[i])
  #   }
  #   return(steps)
  # }
  # 
  # 
  # step_map_test <- map(amt_steps$steps, extract_annual_covs)
  # 
  # 
  # 
  # for(i in 1:length(cov_stacks_ts)){
  #   if(max(time(cov_stacks_ts[[i]])))
  #   test_steps2 <- extract_covariates_var_time(test_steps2,  covariates = cov_stacks_ts[[i]],
  #                                              when = "any", 
  #                                              where = "end",
  #                                              max_time = years(4),
  #                                              name_covar = stack_names[i])
  # }
  
#########################################################################
##
## 9. Create column indicating which steps were during active dispersals
##
##########################################################################
  #Import dispersal dates
  dispersal_dates <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_reduced_10-02-2023.csv") %>% 
    mutate(disp_date_nsd = mdy(disp_date_nsd)) %>% select(-sex)
  
  #Join dispersal dates to locs_nested by animal name
 amt_steps_all_covs <- amt_steps_all_covs %>% 
    left_join(dispersal_dates, by=join_by(animal_id))
  
  #trim dispersal tracks to dispersal dates
  # for (i in 1:nrow(locs_nested)){
  #   if(!is.na(locs_nested$disp_date_nsd[i])){
  #     locs_nested$tracks[[i]] <- locs_nested$tracks[[i]] %>% filter(t_ >= locs_nested$disp_date_nsd[i])
  #   }
  # }
  
  #Create column identifying locations during active dispersal ***
  for (i in 1:nrow(amt_steps_all_covs)){
    if(!is.na(amt_steps_all_covs$disp_date_nsd[i])){
      amt_steps_all_covs$steps[[i]] <- amt_steps_all_covs$steps[[i]] %>%
        mutate(dispersing = case_when(t2_ >= amt_steps_all_covs$disp_date_nsd[i] ~ TRUE,
                                      t2_ < amt_steps_all_covs$disp_date_nsd[i] ~ FALSE))
    }
  }
  
  #make sure it worked
  # locs_nested %>% filter(dispersal_status=="disperser") %>% 
  #   select(animal_id, tracks, disp_date_nsd) %>% 
  #   unnest(cols=tracks) %>% 
  #   group_by(animal_id) %>% 
  #   summarize(min = min(t_))
  # 
 
 #clean up
 amt_steps_all_covs <- amt_steps_all_covs %>% 
 unnest(cols=c(steps)) %>% 
 relocate(unique_step, .after = step_id_) %>% 
 ungroup()
 
#########################################################################
##
## 10. Remove individuals with too little data for inference
##
##########################################################################

#Remove individuals with < 30steps and/or < 20 days of steps: removes 8 individuals. not sure yet if this is necessary

removal_list <- amt_steps_all_covs %>% 
  group_by(animal_id) %>%
  mutate(date_range=interval(start=min(t1_), end=max(t2_)),
         step_days = as.duration(date_range)/ddays(1)) %>%
  distinct(unique_step, .keep_all = TRUE) %>%
  summarize(n_steps=n(), step_days = round(first(step_days), 0)) %>%
  #filter(n_steps<=30) %>% 
  #filter(step_days <= 30)
  filter(n_steps<30 | step_days <= 20) %>%
pull(animal_id)

amt_steps_all_covs <- amt_steps_all_covs %>%
  filter(!(animal_id %in% removal_list))


#########################################################################
##
## 11. Covariate transformations
##
##########################################################################

### Aspect ###

#convert aspect from degrees to radians and calculate northing and easting variables with cos and sin transformations respectively
# steps_unscaled <- amt_steps %>% 
#   mutate(aspect_rad = (pi*aspect)/180, .after=aspect) %>%
#   mutate(northing = cos(aspect_rad),
#          easting = sin(aspect_rad), .after=aspect_rad) %>% 
#   rename(aspect_deg = aspect)
  

### USFS Land Cover ###

#round mean values for landuse and landcover to nearest integer
steps_unscaled <- amt_steps_all_covs %>% 
  mutate(land_cover_usfs_annual = round(land_cover_usfs_annual),
         land_use_usfs_annual = round(land_use_usfs_annual))

#define land cover categories
steps_unscaled <- steps_unscaled %>% 
  mutate(land_cover_usfs_annual_str = case_when(land_cover_usfs_annual == 1 ~ "trees",
                                            land_cover_usfs_annual== 2 ~ "tall_trees_shrubs",
                                            land_cover_usfs_annual == 3 ~ "tree_shrub_mix",
                                            land_cover_usfs_annual == 4 ~ "gfh_tree_mix",
                                            land_cover_usfs_annual == 5 ~ "barren_tree_mix",
                                            land_cover_usfs_annual == 6 ~ "tall_shrubs",
                                            land_cover_usfs_annual == 7 ~ "shrubs",
                                            land_cover_usfs_annual == 8 ~ "gfh_shrub_mix",
                                            land_cover_usfs_annual == 9 ~ "barren_shrub_mix",
                                            land_cover_usfs_annual == 10 ~ "gfh",
                                            land_cover_usfs_annual == 11 ~ "barren_gfh_mix",
                                            land_cover_usfs_annual == 12 ~ "barren_impervious",
                                            land_cover_usfs_annual == 13 ~ "snow_ice",
                                            land_cover_usfs_annual == 14 ~ "water",
                                            land_cover_usfs_annual == 15 ~ NA_character_),
         .after=land_cover_usfs_annual)


### USFS Land Use ###

#define land use categories
steps_unscaled <- steps_unscaled %>% 
  mutate(land_use_usfs_annual_str = case_when(land_use_usfs_annual == 1 ~ "agriculture",
                                              land_use_usfs_annual == 2 ~ "developed",
                                              land_use_usfs_annual == 3 ~ "forest",
                                              land_use_usfs_annual == 4 ~ "non_forest_wetland",
                                              land_use_usfs_annual == 5 ~ "other",
                                              land_use_usfs_annual == 6 ~ "rangeland_pasture",
                                              land_use_usfs_annual == 7 ~ NA_character_),
         .after=land_use_usfs_annual)


######### Add Season Covariates ###

#hunting season (deer and elk): Sep 15 - Nov 15 (WDFW)
#calving season (deer and elk): May 15 - July 1 (WDFW)
#wet season: October - April (en.climate data.org)
#dry season: May - September (en.climate data.org)

# steps_unscaled <- steps_unscaled %>% 
#   mutate(season = case_when(month(t2_) >= 10 | month(t2_) <= 4 ~ "wet",
#                             month(t2_) < 10 & month(t2_) > 4 ~ "dry"),
#          hunting_season = ifelse(yday(t2_) %in% c(258:319), "yes", "no"),
#          calving_season = ifelse(yday(t2_) %in% c(135:182), "yes", "no"))


# check to make sure categorization worked correctly
# filter(year(t2_)== 2020) %>% 
#   group_by(season) %>% summarize(first = min(t2_), last = max(t2_))



#########################################################################
##
## 12. Create final step data frame for export
##
##########################################################################

## Resume here: reorder, rename, round, format, etc.


#add unique id and rearrange fields to final format
steps_final <- steps_unscaled %>% 
  mutate(unique_id = 1:nrow(steps_unscaled), .before= animal_id) %>% 
  select(unique_id:burst_,
         step_id_,
         unique_step,
         case_,
         x1_:sl_,
         log_sl_,
         ta_,
         cos_ta_,
         t1_:dt_,
         gpp,
         npp,
         ndvi,
         evi,
         tree_cover_hansen,
         perc_tree_cover:perc_nonveg,
         land_use_usfs,
         land_cover_usfs,
         precip,
         dist_water,
         elevation:aspect,
         northing,
         easting,
         tpi,
         tri,
         roads_hii,
         popdens_hii,
         landuse_hii,
         infra_hii,
         rails_hii,
         power_hii,
         season:calving_season,
         dispersing:disp_qual
         ) %>% 
  mutate(across(roads_hii:power_hii, round))


#scale covariates (optional)
# amt_steps_scaled <- amt_steps %>%
#   mutate(across(elev_start:landuse_hii_end, scale)) %>%
#   mutate(across(elev_start:landuse_hii_end, as.numeric))


#write_csv(steps_final, "data/Location_Data/Steps/2h_steps_unscaled_no_imp_10-02-2023.csv")




################################ GRAVEYARD #################################


# extract_annual_covs <- function(steps, cov_stack_list){
#   steps2 <- steps %>%
#     mutate(year = as.factor(year(t1_)), .before = t1_)
#   
#   steps_split <- split(steps2, steps2$year)
#   
#   extract_fun <- function(steps_split, year){
#     stack <- cov_stack_list[[paste0("cov_stack_", as.character(year))]]
#     covs <- steps_split %>% extract_covariates(stack)
#     
#     return(covs)
#   }
#   
#   steps_split_covs <- list()
#   for(i in 1:length(steps_split)){
#     names <- names(steps_split)
#     steps_split_covs[[i]] <- extract_fun(steps_split[[i]], names[i])
#   }
#   
#   remove_year_names <- function(steps_split_covs, start_name, end_name){
#     names_to_change <-  names(steps_split_covs)[start_name:end_name]
#     names(steps_split_covs)[start_name:end_name]<- substr(names_to_change, 1, nchar(names_to_change) - 5)
#     return(steps_split_covs)
#     
#   }
#   
#   covs_renamed <- map(steps_split_covs, remove_year_names, start_name = 18, end_name = 38)
#   
#   steps_covs_final <- bind_rows(covs_renamed)
#   
#   return(steps_covs_final)
#   
# } #Not working 
# 
# amt_steps_all_covs <- map(amt_steps$steps, extract_annual_covs, cov_stack_list = cov_stacks)
# 
# test <- bind_rows(steps_test)



#add timestamps to annual stacks

# add_timestamp_to_stack <- function(stack){
#   img_names <- names(stack)
#   year <- str_sub(first(img_names), start = -4)
#   ts <- ymd_hms(paste0(year, "-01-01 00:00:00"))
#   terra::time(stack) <- rep(ts, length(img_names))
#   
#   return(stack)
# }
# cov_stacks_annual <- cov_stacks
# cov_stacks_annual[["static_stack"]] <- NULL
# cov_stacks_annual <- map(cov_stacks_annual, add_timestamp_to_stack)
# 
# 
# annual_rasts <- list()
# for(i in 1:length(names(cov_stacks_annual))){
#   stack <- cov_stacks_annual[[i]]
#   
#   band_rasts <- list()
#   for(j in 1:length(names(stack))){
#     band_rasts[[names(stack)[j]]] <- stack[[j]]
#   }
#   
# annual_rasts[[names(cov_stacks_annual[i])]] <- band_rasts
# }
# 
# 
# cov_stacks_combined <- rast(unlist(annual_rasts))
# 
# test <-extract_covariates_var_time(steps2, cov_stacks_combined, max_time = years(1), when = "after")

#


# #remove flagged duplicates from dataframe
# locs_nested_no_dupes <- locs_nested %>% 
#   select(-c(data, sr)) %>% 
#   unnest(cols=c(tracks)) %>% 
#   filter(duplicate_==FALSE) %>% 
#   ungroup() %>%
#   nest_by(animal_id, sex, dispersal_status, median_sr) %>% 
#   unnest(cols=median_sr) %>% 
#   rename(tracks=data)

# tally <- vector()
# for (i in 2:length(sims_marked$imp_status)){
#   tally[[1]] <- 0
#   if(sims_marked$imp_status[[i-1]] == sims_marked$imp_status[[i]]){
#     tally[i] <- 1
#   } else{
#     tally[i] <- 0
#   }
# }

## splitgap method
# source("r_scripts/utility_functions.R")
# split_gaps <- split_at_gap(ctmm_dat[[i]] %>% mutate(animal_id == ctmm_dat[[i]]@info$identity), max_gap = 60*24)


# #summary(ctmm_fit)
# source("r_scripts/utility_functions.R")
# 
# ctmm_dat[[1]]$animal_id <-  ctmm_dat[[1]]@info$identity
# split_gaps <- split_at_gap(data=ctmm_dat[[1]], max_gap = 60*24)
# 
# track_list <- unique(split_gaps$ID_old)
# 
# for (i in 1:length(track_list)){
# ctmm_sims <- list()
# data <-  split_gaps %>% filter(ID_old == track_list[[i]])
# ctmm_sims[[i]] <- ctmm::simulate(object = mod, data = data, complete=TRUE)
# }
# 
# full_track <- bind_rows(ctmm_sims)

#remove steps with endpoints in water
# sf <- part2 %>% 
#   st_as_sf(coords=c("x2_", "y2_"), crs=5070, remove = FALSE)
# 
# water_polys_cropped <- sf::st_buffer(test_sf, dist = 10000) %>%
#   sf::st_union() %>%
#   sf::st_convex_hull() %>%
#   sf::st_intersection(water_polys_filtered) %>%
#   st_collection_extract('POLYGON') %>%
#   st_union() %>% 
#   st_sf()
# 
# sf2 <- sf %>% 
#   mutate(in_water = st_intersects(sf, water_polys_cropped, sparse=FALSE)) %>% 
#   filter(in_water==FALSE) %>% 
#   as_tibble() %>% 
#   select(-c(geometry, in_water))

# steps_calc <- function(x) {
#   part1 <- x %>% 
#     amt::track_resample(rate = hours(resample_int_hours), tolerance = minutes(tolerance_mins)) %>% #resample 
#     amt::filter_min_n_burst() %>% #divide into bursts of min 3 pts
#     amt::steps_by_burst() %>%  #convert to steps
#     dplyr::filter(sl_ >=100) %>% 
#     dplyr::filter(sl_<=20000) #remove stationary and unreasonable step lengths
#   
#   return(part1)
# }
# 
# steps_sample <- function(x) {
#   part2 <- x %>% 
#     amt::random_steps(n_control = rand_steps) %>% #generate random steps per used step
#     amt::extract_covariates(cov_stack, where = "end") %>% #extract covariates at end of step
#     amt::time_of_day(include.crepuscule = FALSE) %>% #calculate time of day for each step--not working
#     mutate(unique_step = paste(burst_,step_id_,sep="_"), .before=step_id_) %>%  #add column for unique step id }
#     mutate(log_sl_ = log(sl_),    
#            cos_ta_ = cos(ta_))
#   
#   return(part2)
# }


# #join identified dispersal dates with dis
# dispersers <- dispersers %>% 
#   select(-c(data)) %>% 
#   left_join(dispersal_dates, by=join_by(animal_id))
# 
# tracks_trimmed <- list()
# 
# for (i in 1:nrow(dispersers)){
#   tracks_trimmed[[i]] <- dispersers$tracks[[i]] %>% filter(t_ >= dispersers$disp_date_nsd[[i]])
# }
# 
# dispersers$tracks <- tracks_trimmed




#filter steps to remove outliers and stationary steps
# steps_unscaled <- amt_steps %>% 
#   filter(sl_ >100 & sl_<20000 ) 
#%>% terra::na.omit()


#make sure no data are missing after omitting steps with missing raster values
#plot_missing(steps_unscaled)