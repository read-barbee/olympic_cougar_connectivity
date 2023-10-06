#### OCP iSSF Module_02: Prepare Screened GPS Data in amt ####

# Author: Read Barbee

# Date:2023-07-12 
#Last updated: 2023-10-02

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

rand_steps <- 20

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

cov_stack <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2_asp.tif")

#rename covariate bands
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
                      "dist_water",
                      "northing",
                      "easting")

#inspect covariate stack--plotting takes a long time. Parallel doesn't seem to be working
# plan(multicore)
# terra::plot(cov_stack, parallel = TRUE)

#########################################################################
##
## 3. Additional Data Cleaning
##
##########################################################################

#Nest locations by animal_id 
locs_nested <- locs_screened %>% 
  nest_by(animal_id, sex, dispersal_status) 



#########################################################################
##
## 4. Convert data to amt tracks
##
##########################################################################

#Create track for each animal

multi_track <- function(d){
  make_track(d, lon_utm, lat_utm, date_time_utc, check_duplicates = TRUE, all_cols = TRUE, crs = 5070) 
  #%>% transform_coords(crs_to = 5070)
}

#create tracks
locs_nested$tracks <-map(locs_nested$data, multi_track)


#subset only disperser locations
dispersers <- locs_nested %>% 
  filter(dispersal_status=="disperser")

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
## 6. Remove duplicate locations (could do this in screening section)
##
##########################################################################
#check for duplicate points in tracks
bind_rows(locs_nested$tracks) %>% get_dupes(deployment_id, t_)


#initialize clusters for parallel computing
# n.clusters = 8
# my.cluster <- parallel::makeCluster(n.clusters, type = "PSOCK")
# doParallel::registerDoParallel(cl = my.cluster, cores = ncores)
# 
# #flag duplicates in tracks as any locations within 15 minutes of each other
# locs_nested$tracks <- foreach(i = 1:nrow(locs_nested), .packages = c("amt", "dplyr"))%dopar%{
#   locs_nested$tracks[[i]] <-locs_nested$tracks[[i]] %>% 
#     flag_duplicates(gamma = minutes(45), DOP = "dop") %>% 
#     filter(duplicate_==FALSE) %>% 
#     flag_duplicates(gamma = minutes(45), DOP = "dop") #flag twice to catch mutiple dupes in a row
# }
# 
# #examine which points were flagged as duplicates
# locs_nested$tracks  %>% filter(duplicate_ ==TRUE)
# 
# #remove duplicates from tracks
# locs_nested$tracks <- foreach (i = 1:nrow(locs_nested), .packages=c("dplyr"))%dopar%{
#   locs_nested$tracks[[i]] <- locs_nested$tracks[[i]] %>% filter(duplicate_==FALSE) 
# }
# #stop parallel computing cluster
# parallel::stopCluster(cl = my.cluster)

#Alternative: flag duplicates in tracks as any locations less than median sampling interval by 15 min
# for (i in 1:nrow(locs_nested)){
#   locs_nested$tracks[[i]] <-locs_nested$tracks[[i]] %>% flag_duplicates(gamma = hours(locs_nested$median_sr[[i]]) - minutes(15))
# }

#########################################################################
##
## 7. Identify dispersal events (only need to do this once and save to csv)
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
## 8. Identify optimal resampling interval
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
## 8. Resample Tracks/Generate Steps/Extract covariate values
##
#######################################################################

#function to convert locatons to steps and extract covariate values
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
  part2 <- part1 %>% 
    amt::random_steps(n_control = rand_steps) %>% #generate random steps per used step
    amt::extract_covariates(cov_stack, where = "end") %>% #extract covariates at end of step
    amt::time_of_day(include.crepuscule = FALSE) %>% #calculate time of day for each step--not working
    mutate(unique_step = paste(burst_,step_id_,sep="_"), .before=step_id_) %>%  #add column for unique step id }
    mutate(log_sl_ = log(sl_),    
           cos_ta_ = cos(ta_))
  
  return(part2)
  }
}

#remove individuals not compatible with the sampling interval
amt_locs <- locs_nested %>% 
  filter(!(animal_id %in% indiv_to_remove))
  # %>% filter(!(animal_id %in% c("Belle", "Cato")))

#troubleshooting
#amt_locs %>% summarize(count = nrow(tracks)) %>% filter(count<10)

#map the step function over each individual and append as a nested dataframe
amt_locs$steps <- map(amt_locs$tracks, steps_calc)



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
## 9. Create column indicating which steps were during active dispersals
##
##########################################################################
  #Import dispersal dates
  dispersal_dates <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_reduced_10-02-2023.csv") %>% 
    mutate(disp_date_nsd = mdy(disp_date_nsd)) %>% select(-sex)
  
  #Join dispersal dates to locs_nested by animal name
 amt_steps <- amt_steps %>% 
    left_join(dispersal_dates, by=join_by(animal_id))
  
  #trim dispersal tracks to dispersal dates
  # for (i in 1:nrow(locs_nested)){
  #   if(!is.na(locs_nested$disp_date_nsd[i])){
  #     locs_nested$tracks[[i]] <- locs_nested$tracks[[i]] %>% filter(t_ >= locs_nested$disp_date_nsd[i])
  #   }
  # }
  
  #Create column identifying locations during active dispersal ***
  for (i in 1:nrow(amt_steps)){
    if(!is.na(amt_steps$disp_date_nsd[i])){
      amt_steps$steps[[i]] <- amt_steps$steps[[i]] %>%
        mutate(dispersing = case_when(t2_ >= amt_steps$disp_date_nsd[i] ~ TRUE,
                                      t2_ < amt_steps$disp_date_nsd[i] ~ FALSE))
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
 amt_steps <- amt_steps %>% 
 unnest(cols=c(steps)) %>% 
 relocate(unique_step, .after = step_id_) %>% 
 ungroup()
 
#########################################################################
##
## 10. Remove individuals with too little data for inference
##
##########################################################################

#Remove individuals with < 30steps and/or < 20 days of steps: removes 8 individuals. not sure yet if this is necessary

removal_list <- amt_steps %>% 
  group_by(animal_id) %>%
  mutate(date_range=interval(start=min(t1_), end=max(t2_)),
         step_days = as.duration(date_range)/ddays(1)) %>%
  distinct(unique_step, .keep_all = TRUE) %>%
  summarize(n_steps=n(), step_days = round(first(step_days), 0)) %>%
  #filter(n_steps<=30) %>% 
  #filter(step_days <= 30)
  filter(n_steps<30 | step_days <= 20) %>%
pull(animal_id)

amt_steps <- amt_steps %>%
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
steps_unscaled <- amt_steps %>% 
  mutate(land_cover_usfs = round(land_cover_usfs),
         land_use_usfs = round(land_use_usfs))

#define land cover categories
steps_unscaled <- steps_unscaled %>% 
  mutate(land_cover_usfs = case_when(land_cover_usfs == 1 ~ "trees",
                                     land_cover_usfs == 2 ~ "tall_trees_shrubs",
                                     land_cover_usfs == 3 ~ "tree_shrub_mix",
                                     land_cover_usfs == 4 ~ "gfh_tree_mix",
                                     land_cover_usfs == 5 ~ "barren_tree_mix",
                                     land_cover_usfs == 6 ~ "tall_shrubs",
                                     land_cover_usfs == 7 ~ "shrubs",
                                     land_cover_usfs == 8 ~ "gfh_shrub_mix",
                                     land_cover_usfs == 9 ~ "barren_shrub_mix",
                                     land_cover_usfs == 10 ~ "gfh",
                                     land_cover_usfs == 11 ~ "barren_gfh_mix",
                                     land_cover_usfs == 12 ~ "barren_impervious",
                                     land_cover_usfs == 13 ~ "snow_ice",
                                     land_cover_usfs == 14 ~ "water",
                                     land_cover_usfs == 15 ~ NA_character_))


### USFS Land Use ###

#define land use categories
steps_unscaled <- steps_unscaled %>% 
  mutate(land_use_usfs = case_when(land_use_usfs == 1 ~ "agriculture",
                                   land_use_usfs == 2 ~ "developed",
                                   land_use_usfs == 3 ~ "forest",
                                   land_use_usfs == 4 ~ "non_forest_wetland",
                                   land_use_usfs == 5 ~ "other",
                                   land_use_usfs == 6 ~ "rangeland_pasture",
                                   land_use_usfs == 7 ~ NA_character_))


######### Add Season Covariates ###

#hunting season (deer and elk): Sep 15 - Nov 15 (WDFW)
#calving season (deer and elk): May 15 - July 1 (WDFW)
#wet season: October - April (en.climate data.org)
#dry season: May - September (en.climate data.org)

steps_unscaled <- steps_unscaled %>% 
  mutate(season = case_when(month(t2_) >= 10 | month(t2_) <= 4 ~ "wet",
                            month(t2_) < 10 & month(t2_) > 4 ~ "dry"),
         hunting_season = ifelse(yday(t2_) %in% c(258:319), "yes", "no"),
         calving_season = ifelse(yday(t2_) %in% c(135:182), "yes", "no"))


# check to make sure categorization worked correctly
# filter(year(t2_)== 2020) %>% 
#   group_by(season) %>% summarize(first = min(t2_), last = max(t2_))



#########################################################################
##
## 12. Create final step data frame for export
##
##########################################################################
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