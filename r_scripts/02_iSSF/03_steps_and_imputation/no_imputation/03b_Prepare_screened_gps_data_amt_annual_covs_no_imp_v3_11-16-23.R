#### OCP iSSF Module_02: Prepare Screened GPS Data in amt ####

# Author: Read Barbee

# Date:2023-07-12 
#Last updated: 2023-11-16

# Purpose: Prepare screened GPS data in amt to fit iSSFs.

# Inputs:
#   •	Screened GPS Locations (corrected for error and habitat bias)
#   •	Deployment list
#   •	Disperser list
#   •	Covariate Stack
#
# Outputs:
#   •	amt_step dataframe with fields for all necessary covariates

#consider parallelizing some of the loops

################################ Libraries #################################
library(tidyverse)
library(terra)
library(amt)
library(janitor)
#library(doParallel)
library(sf)
library(mapview)



################################ User-Defined Parameters #################################

project_crs <- 5070 #NAD83 Albers Equal Area Projection

resample_int_hours <- 2 #interval to resample tracks

tolerance_mins <- 15

rand_steps <- 10

study_area_buffer_dist <- 500 #meters. 500m allows for walking out on beaches at low tide

cov_folder_path <- "/Users/tb201494/Desktop/annual_cov_stacks_1km_buffer"

#plan(multisession, workers = 8)

#########################################################################
##
## 1. Import and format Barrier Polygons
##
##########################################################################

## Study area boundary ##
op_poly <- st_read("data/Habitat_Covariates/study_area_polys/ocp_study_area_poly_wa_only_10-30-23.shp")
op_poly_buffer <- st_buffer(op_poly, study_area_buffer_dist)

## Water Polygons ##
water_polys <- st_read("data/Habitat_Covariates/washington_water_polygons/op_water_polys_with_salt.shp") %>% st_transform(crs = project_crs)

## Freshwater only ##
water_polys_cropped <-  water_polys %>%
  filter(WB_PERIOD_ =="PER" | OBJECTID != 1) %>% #only include permanent water bodies with area > 100 m2
  filter(SHAPEAREA >= 100) %>%
  sf::st_crop(op_poly) #crop to study area
## Freshwater only (dissolved) ##
water_polys_mask <- water_polys_cropped %>%
  sf::st_union() %>% #dissolve polygons into single vector mask layer
  st_sf()

## Saltwater only ##
ocean <- water_polys %>% filter(OBJECTID == 1)



## only for dispersal date identification (one time use) ##
# # Mountain lion deployments 
# deployments <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Deployments/collar_deployments_master_7-11-2023.csv")
# 
# # Dispersal inventory from Teams 
# dispersals <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_10-02-2023.csv")

#########################################################################
##
## 2. Import and format location data
##
##########################################################################
# Mountain lion location data 
locs_screened <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_10-02-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = with_tz(date_time_local, tzone = "US/Pacific"))

#check for duplicate locations
get_dupes(locs_screened, deployment_id, date_time_utc)

#remove_locations outside of study area if not already done in screening stage (should remove 2 pts)
locs_screened <- locs_screened %>% 
  st_as_sf(coords = c("lon_utm", "lat_utm"), crs = project_crs, remove=FALSE) %>% 
  mutate(intersects_study_area = lengths(st_intersects(., op_poly_buffer)) > 0) %>% 
  filter(intersects_study_area==TRUE) %>% 
  as_tibble() %>% 
  select(-c(intersects_study_area, geometry))

#inspect points
# pts_out <- locs_screened2 %>% as.data.frame() %>% select(-geometry)
# write_csv(pts_out, "pts_out.csv")

#########################################################################
##
## 3. (Optional) Remove locations in water (could do this during screening phase)
##
##########################################################################

#CAUTION: may remove many points not actually in water because of polygon error

# locs_screened2 <- locs_screened %>%
#   st_as_sf(coords=c("lon_utm", "lat_utm"), crs=project_crs, remove=FALSE) %>%
#   mutate(intersects_freshwater = lengths(st_intersects(., water_polys_mask)) > 0 ) %>%
#   mutate(intersects_ocean = lengths(st_intersects(., ocean)) > 0 ) %>%
#   mutate(intersects_water = case_when(intersects_ocean==FALSE & intersects_freshwater==FALSE ~FALSE,
#                                       .default=TRUE)) %>%
#   mutate(intersects_study_area = lengths(st_intersects(., op_poly_buffer)) > 0) %>%
#   filter(intersects_water==TRUE & intersects_study_area==TRUE)
# 
# mapview(list(locs_screened2, water_polys_mask, ocean))
  


#########################################################################
##
## 4. Import and format covariate data
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
static_stack <- rast("/Users/tb201494/Desktop/1km_buffer/static_stack_1km_buffer_single_tpi.tif")
mtpi <- rast("/Users/tb201494/Desktop/1km_buffer/static/mtpi_1km_buffer.tif")

static_stack$mtpi <- mtpi

static_stack <- static_stack[[c("elevation", "slope", "aspect", "aspect_northness", "aspect_eastness", "tpi", "mtpi", "tri", "distance_water")]]

#########################################################################
##
## 5. Convert locations to amt tracks
##
##########################################################################

#Create track for each animal

multi_track <- function(d){
  make_track(d, lon_utm, lat_utm, date_time_utc, check_duplicates = TRUE, all_cols = TRUE, crs = project_crs) 
  #%>% transform_coords(crs_to = project_crs)
}

#Nest locations by animal_id 
locs_nested <- locs_screened %>% 
  nest_by(animal_id, sex, dispersal_status) 


#create tracks
locs_nested$tracks <-map(locs_nested$data, multi_track)

#check for duplicate points in tracks
bind_rows(locs_nested$tracks) %>% get_dupes(deployment_id, t_)

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
  select(-c(data, tracks, median_sr)) %>% 
  unnest(cols = c(sr)) 

summary(sampling_rates)


#########################################################################
##
## 7. Identify dispersal events (ONLY NEED TO DO THIS ONCE and save to csv)
##
##########################################################################
# 
# # Map track
# filt <- dispersers %>% 
#   select(data) %>% 
#   unnest(cols=data) %>% 
#   filter(animal_id=="Archie") %>% 
#   #filter(date_time_utc>= ymd("2020/03/15")) %>% 
#   sf::st_as_sf(coords = c("lon_utm", "lat_utm"), crs = project_crs) %>% ungroup()
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
## 8. Identify optimal resampling interval. Remove individuals that can't be resampled
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
## 9. Resample Tracks/Generate Steps
##
#######################################################################

#function to convert locatons to steps. Won't generate random steps in water or outside the study area
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
      amt::random_steps(n_control = rand_steps*5) %>% #generate random steps per used step
      #amt::extract_covariates(cov_stack, where = "end") %>% #extract covariates at end of step
      amt::time_of_day(include.crepuscule = FALSE) %>% #calculate time of day for each step--not working
      mutate(unique_step = paste(burst_,step_id_,sep="_"), .before=step_id_) %>%  #add column for unique step id }
      mutate(log_sl_ = log(sl_),    
             cos_ta_ = cos(ta_))
    
    #determine whether each random step intersects a water body and/or is within the study area. this is the bottleneck
    part3 <- part2 %>% 
      st_as_sf(coords = c("x2_", "y2_"), crs = project_crs, remove=FALSE) %>% 
      mutate(intersects_freshwater = lengths(st_intersects(., water_polys_mask)) > 0 ) %>% #this is the bottleneck (water polys mask. works with just ocean or just water polys, but not both together)
      mutate(intersects_ocean = lengths(st_intersects(., ocean)) > 0 ) %>%
      mutate(intersects_water = case_when(intersects_ocean==FALSE & intersects_freshwater==FALSE ~FALSE,
                                          .default=TRUE)) %>%
      mutate(intersects_study_area = lengths(st_intersects(., op_poly_buffer)) > 0) %>%
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
      filter(intersects_water == FALSE) %>% 
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


#map the step function over each individual and append as a nested dataframe
amt_locs$steps <- map(amt_locs$tracks, steps_calc)

#select the relevant columns and unnest the steps column
amt_steps <- amt_locs %>% 
  select(animal_id:dispersal_status,
         steps) %>% 
  filter(length(steps)>0) # removes any individuals without enough steps to fit distributions (should only be 1)

#map resulting points and check for errors
# bind_rows(amt_steps$steps) %>%
#   filter(intersects_water==TRUE) %>%
#   #filter(intersects_study_area==FALSE) %>%
#   filter(case_==TRUE) %>%
#   st_as_sf(coords=c("x2_", "y2_"), crs = project_crs) %>%
#   mapview::mapview()


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
## 10. Extract covariate values
##
#######################################################################

amt_steps_all_covs <- amt_steps

#static covariates
amt_steps_all_covs$steps <- map(amt_steps_all_covs$steps, extract_covariates, covariates = static_stack)

#annual covariates (make sure to adjust column indices if more are added)
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

#########################################################################
##
## 11. Interpolate values in water-masked extent of GEE products that are actually on land
##
#########################################################################

#points near bodies of water don't have hii values, npp, gpp, or perc veg values because of coarse water-masking in GEE products. Solve by interpolation

#list of covariates with coarse water masking to interpolate values for
covs_to_interp_300m <- c("perc_nonveg", "perc_nontree_veg", "perc_tree_cov", "hii_annual", "roads_hii_annual", "infra_hii_annual", "landuse_hii_annual", "popdens_hii_annual", "mtpi")

#Interpolate precip, npp and gpp with 600 m radius b/c of coarser scale
covs_to_interp_600m <- c("precip_annual", "npp_annual", "gpp_annual")

#function to interpolate values from specified buffer radius around points
interp_cov_vals <- function(cov, df, buffer_size){

#prep data and add na status column
col_name <- cov 
column <- df[[col_name]] 

#bypass function if there are no NA values
if(sum(is.na(column)) == 0){
  return(df)
} else if(sum(is.na(column)) > 0){
df2 <- df %>% 
  mutate(na_status = is.na(!!sym(col_name))) 

#split data by na_status
na_split <- split(df2, df2$na_status)

#select the na values for interpolation
na_only <- na_split$`TRUE` %>% 
  st_as_sf(coords=c("x2_", "y2_"), crs = project_crs, remove = FALSE) %>% 
  st_buffer(buffer_size)

#split na_values by year
na_year <- split(na_only, na_only$year)

#interpolate missing values as average of cells within buffer radius for each year
for(i in 1:length(na_year)){
 if(nrow(na_year[[i]])>0){
  year <- names(na_year)[i]
  stack <- cov_stacks[[paste0("cov_stack_", as.character(year))]]
  rast <- stack[[str_detect(names(stack), coll(col_name))]]
  
  #extract values
  na_year[[i]][[col_name]] <- terra::extract(rast, na_year[[i]], fun=function(x){mean(x, na.rm=T)})[,2]
  #replace NaN values with NA
  na_year[[i]][[col_name]] <- ifelse(is.nan(na_year[[i]][[col_name]]), NA, na_year[[i]][[col_name]])
  }
}

#substitute interpolated values and recombine datafarme
na_only_sub <- bind_rows(na_year) %>% as.data.frame() %>% select(-geometry)
na_split$`TRUE` <- na_only_sub
out <- bind_rows(na_split) %>% arrange(t2_, desc(case_)) #animal_id

return(out)
}
}

#function to map interpolation function across all covariates
cov_interp_map <- function(steps_interp){
#Interpolate missing values at 300m radius ~13 min
for(i in 1:length(covs_to_interp_300m)){
  steps_interp <- interp_cov_vals(covs_to_interp_300m[i], df = steps_interp, buffer_size = 300)
print(paste0(i, "/", length(covs_to_interp_300m)))
}

#Interpolate missing values at 600m radius
for(i in 1:length(covs_to_interp_600m)){
  steps_interp<- interp_cov_vals(covs_to_interp_600m[i], df = steps_interp, buffer_size = 600)
  print(paste0(i, "/", length(covs_to_interp_600m)))
}

#re-nest step data
# steps_interp<- steps_interp %>% 
#   select(-na_status) #%>% 
  #nest_by(animal_id, sex, dispersal_status) %>% 
  #rename(steps = data)
return(steps_interp)
}

steps_interp <- amt_steps_all_covs

#map interpolation across all indidivual step frames ~ 15 min
steps_interp$steps <- map(amt_steps_all_covs$steps, cov_interp_map)

#alternative for loop to get more detailed progress
# system.time(
# for(i in 1:length(amt_steps_all_covs$steps)){
#   steps_interp$steps[i] <- cov_interp_map(amt_steps_all_covs$steps[i])
#   print(paste0("blah_", i, "/", length(amt_steps_all_covs$steps)))
#   
# })


#Check remaining NA points
steps_interp %>% unnest(cols=c(steps)) %>% select(-na_status) %>% DataExplorer::plot_missing()

#########################################################################
##
## 12. Create column indicating which steps were during active dispersals
##
##########################################################################
#Import dispersal dates
dispersal_dates <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_reduced_10-02-2023.csv") %>% 
  mutate(disp_date_nsd = mdy(disp_date_nsd)) %>% select(-sex)

#Join dispersal dates to locs_nested by animal name
amt_steps_all_covs_interp <- steps_interp %>% 
  left_join(dispersal_dates, by=join_by(animal_id))

#Create column identifying locations during active dispersal ***
for (i in 1:nrow(amt_steps_all_covs_interp)){
  if(!is.na(amt_steps_all_covs_interp$disp_date_nsd[i])){
    amt_steps_all_covs_interp$steps[[i]] <- amt_steps_all_covs_interp$steps[[i]] %>%
      mutate(dispersing = case_when(t2_ >= amt_steps_all_covs_interp$disp_date_nsd[i] ~ TRUE,
                                    t2_ < amt_steps_all_covs_interp$disp_date_nsd[i] ~ FALSE))
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
amt_steps_all_covs_interp <- amt_steps_all_covs_interp %>% 
  unnest(cols=c(steps)) %>% 
  relocate(unique_step, .after = step_id_) %>% 
  ungroup()

#########################################################################
##
## 13. Remove individuals with too little data for inference
##
##########################################################################

#Remove individuals with < 30steps and/or < 20 days of steps: removes 8 individuals. not sure yet if this is necessary

removal_list <- amt_steps_all_covs_interp %>% 
  group_by(animal_id) %>%
  mutate(date_range=interval(start=min(t1_), end=max(t2_)),
         step_days = as.duration(date_range)/ddays(1)) %>%
  distinct(unique_step, .keep_all = TRUE) %>%
  summarize(n_steps=n(), step_days = round(first(step_days), 0)) %>%
  #filter(n_steps<=30) %>% 
  #filter(step_days <= 30)
  filter(n_steps<30 | step_days <= 20) %>%
  pull(animal_id)

amt_steps_all_covs_interp <- amt_steps_all_covs_interp %>%
  filter(!(animal_id %in% removal_list))

#########################################################################
##
## 11. Covariate transformations
##
##########################################################################

### USFS Land Cover ###
#round mean values for landuse and landcover to nearest integer
steps_unscaled <- amt_steps_all_covs_interp %>% 
  mutate(land_cover_usfs_annual = round(land_cover_usfs_annual),
         land_use_usfs_annual = round(land_use_usfs_annual))


#########################################################################
##
## 12. Create final step data frame for export
##
##########################################################################
#add unique id and rearrange fields to final format
steps_final <- steps_unscaled %>% 
  mutate(unique_id = 1:nrow(steps_unscaled), .before= animal_id) %>% 
  select(-na_status)
  


#inspect points
steps_final %>% filter(is.na(popdens_hii_annual)) %>% filter(case_==FALSE) %>% 
  st_as_sf(coords=c("x2_", "y2_"), crs=project_crs, remove=FALSE) %>% mapview()



#write_csv(steps_final, "data/Location_Data/Steps/2h_steps_unscaled_no_imp_annual_cov_11-16-2023.csv")




