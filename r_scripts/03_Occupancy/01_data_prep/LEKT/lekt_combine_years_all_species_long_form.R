#### LEKT Cam det histories ####

# Author: Read Barbee

# Date:2024-10-16 

# Purpose:

################################ libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)

#########################################################################
##
## 1. Combine detections for all species across years
##
##########################################################################
#species detections
lekt_det_2019 <- read_csv("data/Camera_Data/all_species/LEKT/2019/SNA_20190703_20200214_2022.09.19_21.23_metadata_tbl.csv") %>% 
  clean_names() %>% 
  select(station,cameras, date_time_original, species) %>% 
  rename(timestamp = date_time_original,
         station_id = station,
         camera_id = cameras) %>% 
  mutate(timestamp = as.character(mdy_hm(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2019", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)


#species detections
lekt_det_2020 <- read_csv("data/Camera_Data/all_species/LEKT/2020/S3071_20200513_20201210_2022.12.12_22.30_metadata_tbl.V2.csv") %>% 
  clean_names() %>% 
  select(station, cameras, date_time_original, species) %>% 
  rename(timestamp = date_time_original,
         station_id = station,
         camera_id = cameras) %>% 
  filter(is.na(station_id) == F ) %>% 
  mutate(timestamp = as.character(mdy_hm(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2020", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)



lekt_det_2021 <- read_csv("data/Camera_Data/2021/LEKT_2021/S3071_20210419_20211209_2022.09.12_20.04_metadata_tbl.csv", col_types = cols(DateTimeOriginal = col_character())) %>% 
  clean_names() %>% 
  select(station, cameras, date_time_original, species) %>% 
  rename(timestamp = date_time_original,
         station_id = station,
         camera_id = cameras) %>% 
  filter(!is.na(station_id)) %>% 
  mutate(timestamp = as.character(ymd_hms(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2021", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)

lekt_det_2022 <- read_csv("data/Camera_Data/all_species/LEKT/2022/2022 OCP Camera Data.csv") %>% 
  clean_names() %>% 
  select(station_num, timestamp, common_name) %>% 
  rename(station_id = station_num,
         species = common_name) %>%
  mutate(camera_id = "Camera1", .after = station_id) %>% 
  filter(!is.na(station_id)) %>% 
  mutate(station_id = str_remove_all(station_id, "Check.*")) %>% 
  mutate(timestamp = as.character(mdy_hm(timestamp))) %>%
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F,
                       too_few = "align_start") %>% 
  mutate(dep_year = "2022", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)


lekt_det_2023 <- read_csv("data/Camera_Data/all_species/LEKT/2023/2023.Data_key_spp_time_date_shifts_complete.csv",
                          col_types = cols(timestamp = col_character())) %>% 
  rename(station_id = Station,
         species = common_name) %>% 
  select(station_id, timestamp, species) %>% 
  filter(!is.na(station_id)) %>% 
  separate_wider_delim(cols = timestamp, 
                       delim = " ", 
                       names = c("date", "time"), 
                       cols_remove = F) %>% 
  mutate(station_id = paste0("Station", station_id),
         camera_id = "Camera1") %>% 
  mutate(dep_year = "2023", 
         act_year = year(timestamp),
         .before = station_id) %>% 
  select(station_id, camera_id, dep_year, act_year, species, timestamp, date, time)


lekt_dets_all <- bind_rows(lekt_det_2019,
                           lekt_det_2020,
                           lekt_det_2021,
                           lekt_det_2022,
                           lekt_det_2023)%>% 
  mutate(grid = "LEKT", .before = dep_year) %>% 
  mutate(station_id = str_replace_all(station_id, 
                                      pattern = coll("Station"),
                                      replacement =coll("LEKT_"))) %>% 
  unite("deployment_id", station_id, dep_year, sep="_", remove=F) %>% 
  unite("station_year", station_id, act_year, sep="_", remove=F) %>% 
  select(grid, deployment_id, station_year, station_id, camera_id, dep_year, act_year, species, timestamp, date, time)

#write_csv(lekt_dets_all, "data/Camera_Data/lekt_det_all_species_2019_2023.csv")

###############################################################################

#2. Import viewshed files

###############################################################################

#2019 viewsheds have data but have not been calculated

#2020
vs_2020 <- read_csv("data/Camera_Data/all_species/LEKT/2020/lekt_2020_viewsheds_edited.csv") %>% 
  rename(station_id = `2020.Sta.Num`, 
         viewshed = Total_Area,
         exclude = Exclude) %>%
  select(station_id, viewshed, exclude) %>% 
  filter(exclude == FALSE) %>% 
  mutate(station_id = str_remove(station_id, "_.*")) %>% 
  select(-exclude)
  
vs_2021 <- read_csv("data/Camera_Data/all_species/LEKT/2021/lekt_2021_viewsheds_edited.csv") %>% 
  rename(station_id = `2021.Sta.Num`, 
         viewshed = Total_Area) %>%
  select(station_id, viewshed) %>%
  mutate(camera_id = case_when(str_detect(station_id, coll("b")) ~ "Camera2",
                               .default = "Camera1"), .after = station_id) %>% 
  mutate(station_id = str_remove(station_id, "a")) %>% 
  mutate(station_id = str_remove(station_id, "b")) %>% 
  mutate(station_id = paste0("Station", station_id))


vs_2022 <- read_csv("data/Camera_Data/all_species/LEKT/2022/lekt_2022_viewsheds_edited.csv") %>% 
  rename(station_id = cam,
         viewshed = area)

vs_2023 <- read_csv("data/Camera_Data/all_species/LEKT/2023/lekt_2023_viewsheds_edited.csv") %>% 
  rename(station_id = `2023.Sta.Num`,
         viewshed = Total_Area)

###############################################################################

#3. Format individual activity sheets

###############################################################################

#camera activity
lekt_act_2019 <- read_csv("data/Camera_Data/all_species/LEKT/2019/2019.Activity.Sheet.csv") %>% 
  select(-c(study, relative_location, camera_id, cam_num:pull_date)) %>%
  rename(station_id = station,
         camera_id = cameras) %>% 
  mutate(viewshed = NA,
         station_id = str_replace(station_id, coll("Station"), coll("LEKT"))) %>% 
  relocate(viewshed, .after = camera_id) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  # mutate(pull_date = case_when(pull_date == "Stolen" ~ NA,
  #                              .default = pull_date)) %>% 
  mutate(#set_up_date = as.character(mdy(set_up_date)),
         #pull_date = as.character(mdy(pull_date)),
         act_date = as.character(mdy(act_date))) %>% 
  relocate(camera_id, .after = station_id) %>% 
  relocate(longitude, .before = latitude) %>% 
  mutate(dep_year = "2019", 
         act_year = year(ymd(act_date)),
         .after = camera_id) %>%
  mutate(survey_id = "LEKT_2019", .before = station_id)



lekt_act_2020 <-read_csv("data/Camera_Data/all_species/LEKT/2020/2020.Camera.Activity.Sheet.IDS.copy.csv") %>% 
  select(-c(study, camera_id, set_up_date:elevation)) %>%
  rename(longitude = x,
         latitude = y,
         station_id = station,
         camera_id = cameras) %>% 
  left_join(vs_2020, by = "station_id") %>%
  relocate(viewshed, .after = camera_id) %>% 
  mutate(station_id = str_replace(station_id, coll("Station"), coll("LEKT"))) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>%
  mutate(#set_up_date = as.character(mdy(set_up_date)),
         #pull_date = as.character(mdy(pull_date)),
         act_date = as.character(dmy(act_date))) %>%
  relocate(camera_id, .after = station_id) %>% 
  mutate(dep_year = "2020", 
         act_year = year(ymd(act_date)),
         .after = camera_id) %>%
  mutate(survey_id = "LEKT_2020", .before = station_id)


lekt_act_2021 <-read_csv("data/Camera_Data/all_species/LEKT/2021/2021.Camera.Activity.Sheet.csv") %>% 
  select(-c(study, camera_id, elevation)) %>% 
  rename(longitude = x,
         latitude = y,
         station_id=station,
         camera_id = cameras) %>%
  left_join(vs_2021, by = c("station_id", "camera_id")) %>%
  relocate(viewshed, .after = camera_id) %>% 
  mutate(station_id = str_replace(station_id, coll("Station"), coll("LEKT"))) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>%
  mutate(#set_up_date = as.character(mdy(set_up_date)),
    #pull_date = as.character(mdy(pull_date)),
    act_date = as.character(mdy(act_date))) %>% 
  relocate(camera_id, .after = station_id) %>% 
  mutate(dep_year = "2021", 
         act_year = year(ymd(act_date)),
         .after = camera_id) %>%
  mutate(survey_id = "LEKT_2021", .before = station_id)



lekt_act_2022 <-read_csv("data/Camera_Data/all_species/LEKT/2022/2022.Activity.Sheet.NEW.csv") %>%
  clean_names() %>% 
  select(-c(study, camera_id, elevation)) %>% 
  rename(longitude = x,
         latitude = y,
         station_id=station,
         camera_id = cameras) %>%
  left_join(vs_2022, by = "station_id") %>%
  relocate(viewshed, .after = camera_id) %>% 
  mutate(station_id = str_replace(station_id, coll("Station"), coll("LEKT"))) %>% 
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(act_date = as.character(mdy(act_date))) %>% 
  relocate(camera_id, .after = station_id) %>% 
  mutate(dep_year = "2022", 
         act_year = year(ymd(act_date)),
         .after = camera_id) %>%
  mutate(survey_id = "LEKT_2022", .before = station_id)
 

lekt_act_2023 <- read_csv("data/Camera_Data/all_species/LEKT/2023/2023.Activity Sheet Only_edited.csv") %>% 
  clean_names() %>%
  select(-c(viewshed:end)) %>%
  rename(longitude = x,
         latitude = y,
         station_id=cam) %>%
  left_join(vs_2023, by = "station_id") %>%
  relocate(viewshed, .after = camera_id) %>% 
  mutate(camera_id = case_when(str_detect(camera_id, coll("_2")) ~ "Camera2",
                               .default = "Camera1")) %>% 
  mutate(camera_id = case_when(station_id == "33_2" ~ "Camera1",
                               station_id == "56_2" ~ "Camera2",
                               .default = camera_id)) %>%
  mutate(station_id = case_when(station_id != "33_2" ~ str_remove(station_id, coll("_2")),
                                .default = station_id)) %>%
  mutate(station_id = paste0("LEKT", station_id)) %>%
  pivot_longer(-c(station_id:viewshed), names_to = "act_date", values_to = "status") %>% 
  mutate(act_date = str_remove_all(act_date, coll("x"))) %>% 
  mutate(act_date = str_replace_all(act_date, coll("_"), coll("/"))) %>% 
  mutate(act_date = as.character(mdy(act_date))) %>% 
  relocate(camera_id, .after = station_id) %>% 
  mutate(dep_year = "2023", 
         act_year = year(ymd(act_date)),
         .after = camera_id) %>%
  mutate(survey_id = "LEKT_2023", .before = station_id)
 
    
    #check locations of rows with same station name but different coords
    # lekt_act_2023 %>% st_as_sf(coords = c(x = "longitude", y = "latitude"), crs = 4326) %>% filter(station_id == "33" | station_id == "33_2") %>% 
    # mapview::mapview()


###############################################################################

#4. Combine activity sheets in long and wide formats

###############################################################################

lekt_act_all <- bind_rows(lekt_act_2019,
                          lekt_act_2020,
                          lekt_act_2021,
                          lekt_act_2022,
                          lekt_act_2023)


lekt_act_all_final_long <- lekt_act_all %>% 
  mutate(camera_id = str_replace(camera_id, coll("Camera"), coll("CAM"))) %>% 
  #unite("station_year", station_id, act_year, sep = "_", remove = F) %>% 
  unite("deployment_id", survey_id, station_id, camera_id, sep="_", remove = F) %>% 
  relocate(survey_id, .before = deployment_id) %>% 
  select(-act_year)

#check for duplicates in deployment id and activity date
lekt_act_all_final_long  %>%
  select(deployment_id, act_date, status) %>% get_dupes(deployment_id, act_date) #%>% distinct(deployment_id)

# #camera_level_test
lekt_act_all_final_long  %>%
  select(deployment_id, act_date, status) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  column_to_rownames("deployment_id") %>%
  as.matrix() %>%
  camtrapR:::camopPlot(., lattice = TRUE)


lekt_act_all_final_wide <- lekt_act_all_final_long %>%
  pivot_wider(names_from = act_date, values_from = status)
  
lekt_act_all_final_wide %>% get_dupes(deployment_id)

################################ Graveyard #################################


# test <- bind_rows(lekt_act_2019,
#                   lekt_act_2020) %>% 
#   pivot_wider(names_from = "act_date", values_from = "status")
# 
# 
# #camera_level_test
# test %>% 
#   unite("station_cam", station_id, camera_id, sep = "_") %>% View()
#   column_to_rownames("station_cam") %>% 
#   select(-c(longitude:n_cams)) %>% 
#   as.matrix() %>% 
#   camtrapR:::camopPlot(., lattice = TRUE)
# 
# 
# #station_level_test
# test %>% 
#   group_by(station_id) %>% 
#   summarize(longitude = first(longitude),
#             latitude = first(latitude),
#             n_cams = n(),
#             across(-c(camera_id:latitude), sum_or_na)) %>% 
#   ungroup() %>% 
#   column_to_rownames("station_id") %>% 
#   select(-c(longitude:n_cams)) %>% 
#   as.matrix() %>% 
#   camtrapR:::camopPlot(., lattice = TRUE)
