library(tidyverse)
library(sf)
library(janitor)

vec_dat <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Raw_Data/Vectronics/Vectronic_complete_download_raw_2023-05-16.csv") %>% clean_names()

lotek_dat <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Raw_Data/Lotek/Lotek_complete_download_raw_2023-05-16.CSV") %>% clean_names()

mort_meta <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Metadata/From_Teams/north_op_cougar_mortalities_v2_6-29-2023.csv") %>% clean_names()


mort_meta <- mort_meta %>% filter(!is.na(mort_date) & !is.na(lat_wgs84)) %>% 
  mutate(mort_date = mdy(mort_date),
         visit_date = mdy(visit_date))

mort_meta_vec <- mort_meta %>% filter(collar_brand=="Vectronic")
mort_meta_lotek <- mort_meta %>% filter(collar_brand=="Lotek")



################################ Vectronics #################################

sample_frames_vec <- list()

for (i in 1:length(mort_meta_vec$collar_id)){
  if(is.na(mort_meta_vec$visit_date[i]==TRUE) | mort_meta_vec$visit_date[i] == mort_meta_vec$mort_date[i]){
  sample_frames_vec[[i]] <- vec_dat %>% filter(collar_id == mort_meta_vec$collar_id[i]) %>% 
    filter(acq_time_utc >= mort_meta_vec$mort_date[i] & acq_time_utc <= (mort_meta_vec$mort_date[i] + day(1) ))
  }
  else{
    sample_frames_vec[[i]] <-vec_dat %>% filter(collar_id == mort_meta_vec$collar_id[i]) %>% 
      filter(acq_time_utc >= mort_meta_vec$mort_date[i] & acq_time_utc <= mort_meta_vec$visit_date[i])
  }
}

sample_frames_vec<- sample_frames_vec %>%
  keep(~ nrow(.) > 0)

vec_frames_sf <- map(sample_frames_vec, st_as_sf, coords = c("longitude_deg", "latitude_deg"), crs=4326, na.fail=FALSE)

mapview::mapview(vec_frames_sf[[1]])



################################ lotek #################################

sample_frames_lotek <- list()

for (i in 1:length(mort_meta_lotek$collar_id)){
  if(is.na(mort_meta_lotek$visit_date[i]==TRUE) | mort_meta_lotek$visit_date[i] == mort_meta_lotek$mort_date[i]){
    sample_frames_lotek[[i]] <- lotek_dat %>% filter(device_id == mort_meta_lotek$collar_id[i]) %>% 
      filter(date_time_gmt >= mort_meta_lotek$mort_date[i] & date_time_gmt <= (mort_meta_lotek$mort_date[i] + day(1) ))
  }
  else{
    sample_frames_lotek[[i]] <-lotek_dat %>% filter(device_id == mort_meta_lotek$collar_id[i]) %>% 
      filter(date_time_gmt >= mort_meta_lotek$mort_date[i] & date_time_gmt <= mort_meta_lotek$visit_date[i])
  }
}

sample_frames_lotek<- sample_frames_lotek %>%
  keep(~ nrow(.) > 0)

lotek_frames_sf <- map(sample_frames_lotek, st_as_sf, coords = c("longitude", "latitude"), crs=4326, na.fail=FALSE)

mapview::mapview(lotek_frames_sf)









## Vec data by mort status
# morts <- vec_dat %>% filter(mortality_status =="Mortality No Radius")
# 
# morts_sf <- morts %>% 
#   select(-animal) %>% 
#   na.omit() %>%  
#   mutate(collar_id = as.factor(collar_id)) %>% 
#   st_as_sf(coords=c("longitude_deg", "latitude_deg"), crs=4326)
# 
# mapview::mapview(morts_sf, zcol="collar_id")



### Append collar model data to master locations

# collar_mods <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Metadata/From_Teams/Formatted_for_R/collar_model_list_partial.csv") %>% clean_names()
# 
# locs <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Source_Files/locations_master/gps_locs_master_6-23-2023.csv", col_types = list(fix_type = col_character()))
# collar_mods <- collar_mods %>% select(collar_id, collar_brand, collar_type, satellite_system) %>% 
#   mutate(collar_id = as.character(collar_id))
# 
# locs2 <- locs %>% left_join(collar_mods, by=join_by(collar_id))
# 
# #write_csv(locs2, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location_Data/Source_Files/locations_master/gps_locs_master_6-29-2023.csv")





################################ GRAVEYARD #################################

# temp<- mort_meta2 %>% select(cougar_id, 
#                                     date_time_mort_signal_received,
#                                     estimated_date_time_of_death, 
#                                     visit_date, 
#                                     carcass_latitude, 
#                                     carcass_longitude, 
#                                     general_location,
#                                     cause_of_death) %>% 
#   rename(animal_id = cougar_id,
#          lat_wgs84 = carcass_latitude,
#          lon_wgs84 = carcass_longitude) %>% 
#   filter(!(animal_id %in% c("Kitten (Uncollared)", "Cougar (Uncollared)", "Cougar  (Uncollared)", "Didi"))) %>% 
#   filter(!(animal_id %in% unique(mort_meta$animal_id))) %>% 
#   filter(animal_id != "Damir") %>% 
#   select(-date_time_mort_signal_received) %>% 
#   rename(mort_date = estimated_date_time_of_death) %>% 
#   mutate(mort_date = ymd(floor_date(mort_date, unit = "day")),
#          visit_date = ymd(floor_date(visit_date, unit = "day")))