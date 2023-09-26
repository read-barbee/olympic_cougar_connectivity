#### Reconstruct and Reroute Paths with CTMM ####

# Author: Read Barbee

# Date:2023-07-12 

# Purpose:


################################ Libraries #################################
library(tidyverse)
library(terra)
library(amt)
library(janitor)
#library(DataExplorer)
#library(furrr)
library(doParallel)
library(sf)
library(mapview)

#########################################################################
##
## 1. Import and format location data
##
##########################################################################
# Mountain lion location data (May 2023)
locs_screened <- read_csv("data/Location_Data/Source_Files/locations_master/gps_locs_dop_screened_7-11-2023.csv", col_types = list(fix_type = col_character())) %>% 
  mutate(date_time_local = with_tz(date_time_local, tzone = "US/Pacific"))

# Mountain lion deployments  
deployments <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Deployments/collar_deployments_master_7-11-2023.csv")

# Dispersal inventory from Teams
dispersals <- read_csv("data/Location_Data/Metadata/From_Teams/Formatted_for_R/Dispersals/dispersals_master_5-11-2023.csv")

#########################################################################
##
## 2. Additional Data Cleaning
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
## 3. Convert data to amt tracks
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
## 4. Examine Sampling Rates
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

#summary(sampling_rates) not working


#########################################################################
##
## 5. Remove duplicate locations (could do this in screening section)
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
## 6. Reconstruct tracks with ctmm (~40 min 8 clusters)
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
#parallel::stopCluster(cl = my.cluster)

#bind simulated frames together
ctmm_sims_unrouted <- bind_rows(ctmm_sims)

####note: simulated frame has 112 fewer observed points than the input locs_screened dataset...
#rerouted frame has 115 fewer points

#write_csv(ctmm_sims_unrouted, "data/Location_Data/imputed_paths/ctmm_sim_imputed_paths_unrouted_7-11-2023.csv")

## map imputed tracks for quality control
og_telem_sf <- ctmm_dat[[1]] %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% 
  mutate(type= "og")

sim_test_sf <- ctmm_sims[[1]] %>% as_tibble() %>% st_as_sf(coords = c("x", "y"), crs = 5070) %>% 
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
## 7. Reroute imputed paths around major water bodies (~ 20 min 8 clusters)
##
##### WARNING: Pathroutr only updates GEOMETRIES ###
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

### WARNING: Pathroutr only updates GEOMETRIES ###

#stop paralllel computing cluster
parallel::stopCluster(cl = my.cluster)

#inspect rerouted paths
mapview::mapview(sims_rerouted[[1]], zcol = "imp_status")

#bind for loop iterations into single dataframe
sims_rerouted_full <- bind_rows(sims_rerouted)

rerouted_coords <- sims_rerouted_full %>% st_coordinates()

sims_rerouted_full_df <- sims_rerouted_full %>% 
  as_tibble() %>%
  mutate(x_rerouted = rerouted_coords[,1],
         y_rerouted = rerouted_coords[,2], .after=t) %>% 
  rename(x_old = x,
         y_old = y,
         longitude_old = longitude,
         latitude_old = latitude) %>% 
  relocate(timestamp, .before=t) %>% 
  select(-geometry)

#write_csv(sims_rerouted_full_df, "data/Location_Data/imputed_paths/ctmm_sim_imputed_paths_rerouted_7-11-2023.csv" )
