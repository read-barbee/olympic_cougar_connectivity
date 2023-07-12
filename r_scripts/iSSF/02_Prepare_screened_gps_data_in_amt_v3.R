#### OCP iSSF Module_02: Prepare Screened GPS Data in amt ####

# Author: Read Barbee

# Date:2023-06-02 

# Purpose: Prepare screened GPS data in amt to fit iSSFs.

# Inputs:
#   •	Screened GPS Locations (corrected for error and habitat bias)
#   •	Deployment list
#   •	Disperser list
#   •	Covariate Stack
#
# Outputs:
#   •	amt_step dataframe with fields for all necessary covariates

#Steps
# •	Import location and covariate data
# •	Set time zones
# •	Formatting: Nest locations by individual and join columns for sex and dispersal status 
# •	Generate amt tracks for each individual
# •	Examine sampling rates and identify resampling interval to retain the most individuals
# •	Remove individuals not compatible with that resampling interval
# •	Resample tracks, generate used and random steps, and extract covariate values for end of each step
# •	Filter steps to remove outliers and stationary steps (100 m < sl_ < 20,000 m )
# •	Remove individuals with < 30steps and/or < 30 days of steps????***
# •	Transform aspect to continuous easting and northing values
# •	Define USFS landcover and landuse categories
# •	Derive covariates for climatic season, hunting season, and calving season
# •	Round all floating hii values to integers
# •	reorder columns and export to csv

#Total individuals: 126 (7/11/2023)

################################ Libraries #################################
library(tidyverse)
library(terra)
library(lubridate)
library(amt)
library(janitor)
library(DataExplorer)
library(furrr)
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
# Mountain lion location data (May 2023)
locs_screened <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_7-11-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = with_tz(date_time_local, tzone = "US/Pacific"))

locs_raw <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_master_7-11-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = with_tz(date_time_local, tzone = "US/Pacific"))

# Mountain lion deployments  (September 2022)
deployments <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Deployments/collar_deployments_master_7-11-2023.csv")

# Dispersal inventory from Teams (3-13-2023)
dispersals <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_5-11-2023.csv")

#########################################################################
##
## 2. Import and format covariate data
##
##########################################################################

cov_stack <- rast("data/Habitat_Covariates/puma_cov_stack_v2/tifs/puma_cov_stack_v2.tif")

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
                      "dist_water")

#inspect covariate stack--plotting takes a long time. Parallel doesn't seem to be working
# plan(multicore)
# terra::plot(cov_stack, parallel = TRUE)

#########################################################################
##
## 3. Additional Data Cleaning
##
##########################################################################

#Filter deploment list to get one row per individual
first_deps <- deployments %>% 
  clean_names() %>% 
  distinct(name, .keep_all = TRUE) %>% 
  mutate(name = trimws(name))

#Get vector of disperser names
disperser_names <- dispersals %>% 
  clean_names() %>% 
  filter(!is.na(name)) %>% 
  distinct(name) %>% 
  pull(name)


#Add dispersal status to deployment list
dem_cats <- first_deps %>% 
  mutate(dispersal_status = case_when(name %in% disperser_names == TRUE ~ "disperser",
                                      name %in% disperser_names == FALSE ~ "resident")) %>% 
  select(name,
         sex,
         dispersal_status) %>% 
  rename(animal_id=name)


#Nest locations by animal_id and add columns for sex and dispersal status
locs_nested <- locs_screened %>% 
  nest_by(animal_id) %>% 
  left_join(dem_cats, by = join_by(animal_id)) %>% 
  select(animal_id, sex, dispersal_status, data)

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
## 6. Examine Sampling Rates
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
  select(-c(data, tracks)) %>% 
  unnest(cols = c(sr)) 

#summary(sampling_rates)


#########################################################################
##
## 6. Remove duplicate locations (could do this in screening section)
##
##########################################################################

#initialize clusters for parallel computing
n.clusters = 8
my.cluster <- parallel::makeCluster(n.clusters, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster, cores = ncores)

#flag duplicates in tracks as any locations within 15 minutes of each other
locs_nested$tracks <- foreach(i = 1:nrow(locs_nested), .packages = c("amt", "dplyr"))%dopar%{
  locs_nested$tracks[[i]] <-locs_nested$tracks[[i]] %>% 
    flag_duplicates(gamma = minutes(45), DOP = "dop") %>% 
    filter(duplicate_==FALSE) %>% 
    flag_duplicates(gamma = minutes(45), DOP = "dop") #flag twice to catch mutiple dupes in a row
}

#remove flagged duplicates from dataframe
locs_nested_no_dupes <- locs_nested %>% 
  select(-c(data, sr)) %>% 
  unnest(cols=c(tracks)) %>% 
  filter(duplicate_==FALSE) %>% 
  ungroup() %>% 
  nest_by(animal_id, sex, dispersal_status, median_sr) %>% 
  unnest(cols=median_sr)

#flag duplicates in tracks as any locations less than median sampling interval by 15 min
# for (i in 1:nrow(locs_nested)){
#   locs_nested$tracks[[i]] <-locs_nested$tracks[[i]] %>% flag_duplicates(gamma = hours(locs_nested$median_sr[[i]]) - minutes(15))
# }

#########################################################################
##
## 5. Reconstruct tracks with ctmm
##
##########################################################################
library(ctmm)

#create telemetry object to fit ctmm models
ctmm_dat <- locs_nested_no_dupes %>% 
  unnest(cols=c(data)) %>%
  mutate(t_ = round_date(t_, unit="hour")) %>% 
  dplyr::select(animal_id, collar_id, t_, lat_wgs84,lon_wgs84) %>% 
  rename(individual.local.identifier = animal_id,
         tag.local.identifier = collar_id,
         timestamp = t_,
         location.long = lon_wgs84,
         location.lat = lat_wgs84) %>% 
  as.telemetry(datum = "+proj=longlat +datum=WGS84 +no_defs +type=crs",
               projection = "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" )



#for loop to fit ctmm model and impute points for each individual (~ 40 min with 8 cores)

system.time(ctmm_sims <- foreach(i = 1:length(ctmm_dat), .packages = c("ctmm", "dplyr")) %dopar%{
  #length(ctmm_dat)
  #initialize list for simulations
  ctmm_sims <- list()
  
  #calculate initial values for ctmm model
  guess <- ctmm.guess(ctmm_dat[[i]], interactive = FALSE)
  
  #perform ctmm model selection
  ctmm_fits<- ctmm.select(ctmm_dat[[i]], guess, verbose=TRUE, cores=8)
  
  #select top fitted ctmm model
  mod <- ctmm_fits[[1]]
  
  #simulate missing points based on the model and the observed data (simulation runs through observed pts)
  ctmm_sim <- ctmm::simulate(object = mod, data = ctmm_dat[[i]], complete=TRUE)

  #classify points as observed or imputed and assign unique group id to consecutive locs of same type
  sims_marked <- ctmm_sim %>% as_tibble() %>% 
    mutate(imp_status = case_when(ctmm_sim$timestamp %in% ctmm_dat[[i]]$timestamp ~ "observed", .default = "imputed")) %>% 
    mutate(group = data.table::rleid(imp_status))
  
  #count the number of consecutive values in each group
  group_counts <- sims_marked %>% group_by(group) %>% count()
  
  max_streak <- 24 / locs_nested_no_dupes$median_sr[[i]]
  
  #join the group counts to the original data frame and filter out groups of consecutive imputed points spanning 24 hr or more 
  sims_marked <- sims_marked %>% 
    left_join(group_counts, by=join_by(group)) %>% 
    filter(imp_status =="observed" | n < max_streak) %>% 
    mutate(animal_id = ctmm_dat[[i]]@info$identity, .before = t) %>% 
    select(-c(group, n))
  
  ctmm_sims[[i]] <- sims_marked
})

#stop paralllel computing cluster
parallel::stopCluster(cl = my.cluster)

#bind simulated frames together
ctmm_sims_unrouted <- bind_rows(ctmm_sims)

#note: simulated frame has 112 fewer observed points than the input locs_screened dataset...

#write_csv(ctmm_sims_unrouted, "data/Location_Data/imputed_paths/ctmm_sim_imputed_paths_unrouted_7-11-2023.csv")

## map imputed tracks for quality control
og_telem_sf <- ctmm_dat[[1]] %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% 
  mutate(type= "og")

sim_test_sf <- ctmm_sims[[36]] %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% 
  mutate(type ="sim")

mapview::mapview(sim_test_sf, zcol="imp_status")

mapview::mapview(bind_rows(sim_test_sf, og_telem_sf), zcol="type")


#Exclude Comet(23) and Cato(20). (too few points)


#might need this later: adds empty rows for all possible time steps before imputation
# all_timesteps <- seq(min(data$date_time_utc), max(data$date_time_utc), by = "2 hours")
# 
# full_trajectory <- data.frame(
#   date_time_utc = all_timesteps)
# 
# data_full <- left_join(full_trajectory, data, by=join_by(date_time_utc))

#########################################################################
##
## 5. Reroute imputed paths around major water bodies
##
##########################################################################

#Import water body polygons
water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp") %>% st_transform(crs = 5070)

#filter water body polygons to only include permanent water bodies
water_polys_filtered <- water_polys %>% filter(WB_PERIOD_ =="PER" ) #|WB_PERIOD_ =="INT"


#For loop to reroute all the imputed paths from ctmm (~ 17 min with 8 clusters)
system.time(sims_rerouted <- foreach(i = 1:length(ctmm_sims), .packages = c("dplyr", "sf", "pathroutr")) %dopar% {
  
  #name the iteration? Not sure what this does
  cat(i," ")
  
  #convert ctmm sim object to sf_object for pathroutr
  sim_sf <- ctmm_sims[[i]] %>% 
    as_tibble() %>% 
    st_as_sf(coords = c("x", "y"), crs = 5070, remove = FALSE)
  
  #Filter water body polygons to retain those within convex hull of all of individual's points
  water_polys_cropped <- sf::st_buffer(sim_sf, dist = 10000) %>%
    sf::st_union() %>%
    sf::st_convex_hull() %>%
    sf::st_intersection(water_polys_filtered) %>%
    st_collection_extract('POLYGON') %>%
    st_union() %>% 
    st_sf()
  
  #create buffer around barrier objects as visgraph for rerouting function (connects all verticies of barrier polygon with Delaunay triangle mesh and removes any edges that cross the barrier). Essentially it creates a roadmap of traversible terrain
  visgraph <- pathroutr::prt_visgraph(water_polys_cropped, buffer = 15)
  
  sims_rerouted <- list()
  path <- pathroutr::prt_trim(sim_sf, water_polys_cropped)
  sims_rerouted[[i]] <- pathroutr::prt_reroute(path, water_polys_cropped, visgraph) %>% 
    pathroutr::prt_update_points(path) 
})

#inspect rerouted paths
mapview::mapview(sims_rerouted[[16]], zcol = "imp_status")

#bind for loop iterations into single dataframe
sims_rerouted_full <- bind_rows(sims_rerouted)

#write_csv(sims_rerouted_full, "data/Location_Data/imputed_paths/ctmm_sim_imputed_paths_rerouted_7-11-2023.csv" )

#########################################################################
##
## 5. Identify dispersal events
##
##########################################################################

# Map track
filt <- dispersers %>% 
  select(data) %>% 
  unnest(cols=data) %>% 
  filter(animal_id=="Al") %>% 
  #filter(date_time_utc>= ymd("2020/03/15")) %>% 
  sf::st_as_sf(coords = c("lon_utm", "lat_utm"), crs = 5070) %>% ungroup()

mapview(filt)

#calculate net squared displacement from first location over time
dispersers$tracks <- map(dispersers$tracks, add_nsd)

#plot nsd over time to identify dispersal date 
nsd_plot <- ggplot(data = dispersers %>% 
                 filter(animal_id=="Al") %>% 
                 select(tracks) %>% 
                 unnest(cols=tracks), 
               #%>% filter(t_>= ymd("2020/03/15")), 
               aes(x = t_, y = nsd_)) +
  geom_point()


plotly::ggplotly(nsd_plot)

#recheck Hana, Lady, Lolli
#Lolli's track not trimmed correctly. Also doesn't seem like much of a dispersal
#Harder to determine dispersal date for females because they don't disperse as far


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

locs_nested_test <- locs_nested

locs_nested_test$steps <- map(locs_nested$tracks, test_fun1)


locs_nested_test %>% 
  filter(nrow(steps)<10) %>% #some individuals with low numbers of steps cause the amt functions to fail
  select(steps)

#make list of individuals to remove (because they can't be resampled at the specified interval)
indiv_to_remove <- locs_nested_test %>% 
  filter(nrow(steps)<10) %>% 
  pull(animal_id)


#1 hour removes 67 individuals 
#2 hour removes 35 individuals 
#3 hours removes 29 individuals 
#4 hours removes 32 individuals 
#5 hours removes 64 individuals 
#6 hours removes 10 individuals 

#6 hours retains the most individuals
# 2 hours may offer the best balance between data resolution and data loss


#########################################################################
##
## 8. Resample Tracks/Generate Steps/Extract covariate values
##
##########################################################################

# #Import water body polygons
# water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation/DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp") %>% st_transform(crs = 5070)
# 
# water_polys_filtered <- water_polys %>% filter(WB_PERIOD_ =="PER" ) #|WB_PERIOD_ =="INT"

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
  
  #hist(amt_steps$ta_)
  
#########################################################################
##
## 9. Create column indicating which steps were during active dispersals
##
##########################################################################
  #Import dispersal dates
  dispersal_dates <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_reduced_7-05-2023.csv") %>% 
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

#Remove individuals with < 30steps and/or < 30 days of steps: removes 4 individuals. not sure yet if this is necessary

removal_list <- amt_steps %>%
  group_by(animal_id) %>%
  mutate(date_range=interval(start=min(t1_), end=max(t2_)),
         step_days = as.duration(date_range)/ddays(1)) %>%
  distinct(unique_step, .keep_all = TRUE) %>%
  summarize(n_steps=n(), step_days = round(first(step_days), 0)) %>%
  #filter(n_steps<=30) %>% 
  #filter(step_days <= 30)
  filter(n_steps<30 | step_days <= 30) %>%
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
steps_unscaled <- amt_steps %>% 
  mutate(aspect_rad = (pi*aspect)/180, .after=aspect) %>%
  mutate(northing = cos(aspect_rad),
         easting = sin(aspect_rad), .after=aspect_rad) %>% 
  rename(aspect_deg = aspect)
  

### USFS Land Cover ###

#round mean values for landuse and landcover to nearest integer
steps_unscaled <- steps_unscaled %>% 
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
         land_cover_usfs,
         precip,
         dist_water,
         elevation:tpi,
         land_use_usfs,
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


#write_csv(steps_final, "data/Location_Data/Steps/2h_steps_unscaled_7-07-2023.csv")




################################ GRAVEYARD #################################

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