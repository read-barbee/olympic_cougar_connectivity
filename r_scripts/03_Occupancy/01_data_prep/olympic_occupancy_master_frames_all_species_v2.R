#### Olympic Occupancy Master Frames All Species ####

# Author: Read Barbee

# Date:2024-10-21 

# Purpose:

################################ Helper functions #################################

fill_missing_dates <- function(data, mode){
  
  date_cols <- data %>% 
    select(matches("\\d")) %>% names()
  
  if(mode == "survey"){
    start_date <- min(ymd(date_cols))
    end_date <- max(ymd(date_cols))
    
  }else if (mode == "all" | is.null(mode)){
    start_year <- year(min(ymd(date_cols)))
    end_year <- year(max(ymd(date_cols)))
    
    start_date <- ymd(paste0(start_year, "/01/01"))
    end_date <- ymd(paste0(end_year, "/12/31"))
  }
  
  date_seq <- seq(start_date, end_date, by = "day") %>% as.character()
  
  for(i in date_seq){
    if(!(as.character(i) %in% date_cols))
      data <- data %>% 
        mutate(!!as.character(i) := NA)
  }
  
  out <- data %>%
    pivot_longer(-c(survey_id:dep_year), names_to = "name", values_to = "value") %>%
    mutate(name = ymd(name)) %>%
    arrange(name) %>%
    #mutate(name = as.character(name)) %>%
    pivot_wider(names_from = name, values_from = value)
  
  return(out)
}

## fill_jdays ##
#fill missing juilan days in each year for stacking the detection history
fill_jdays <- function(data){
  
  date_cols <- data %>%
    select(matches("\\d")) %>% names()
  
  for(i in 1:365){
    if(!(as.character(i) %in% date_cols))
      data <- data %>%
        mutate(!!as.character(i) := NA)
  }
  
  out <- data %>%
    pivot_longer(matches("\\d"), names_to = "name", values_to = "value") %>%
    mutate(name = as.numeric(name)) %>%
    arrange(name) %>%
    pivot_wider(names_from = name, values_from = value)
  
  return(out)
}

## fill_jdays ##
#fill missing juilan days in each year for stacking the detection history- not tested
fill_periods<- function(data, nper = 365){
  
  date_cols <- data %>%
    select(matches("\\d")) %>% names()
  
  for(i in 1:nper){
    if(!(as.character(i) %in% date_cols))
      data <- data %>%
        mutate(!!as.character(i) := NA)
  }
  
  out <- data %>%
    pivot_longer(matches("\\d"), names_to = "name", values_to = "value") %>% 
    mutate(name = as.numeric(name)) %>%
    arrange(name) %>%
    pivot_wider(names_from = name, values_from = value)
  
  return(out)
}

#stack wide-form detection history that spans multiple years
stack_detection <- function(det_hist, c_names){
  #change the column names back to dates
  
  det_hist_named <- det_hist %>% 
    rename_with(.cols= matches("\\d"), function(x){c_names}) 
  
  #pivot to long form
  pivot <- det_hist_named %>%
    as.data.frame() %>% 
    pivot_longer(matches("\\d"), names_to = "act_date", values_to = "detection") %>% 
    mutate(act_date = ymd(act_date),
           jday = yday(act_date)) %>% 
    mutate(act_year = as.factor(year(act_date))) 
  
  #split by activity year, pivot each year frame wider, and fill in missing julian days from 1-365
  split_fill <- pivot %>% 
    split(.$act_year) %>% 
    map(., function(x){
      x %>% select(-c(act_year, act_date)) %>% 
        pivot_wider(names_from = jday, values_from = detection)
    }) %>% 
    map(., fill_jdays)
  
  
  #for each year frame, filter only for deployments that are labelled with that year
  for(i in 1:length(split_fill)){
    year_name <- names(split_fill)[i]
    
    split_fill[[i]] <- split_fill[[i]] %>% 
      filter(str_detect(deployment_id, pattern = year_name))
  }
  
  #bind the individual annual frames together
  out <- split_fill %>% 
    bind_rows() #%>% 
}


nper = ceiling(length(period_start_dates)/7)

#stack wide-form detection history that spans multiple years- not working
stack_detection_non_daily <- function(det_hist, nper){
  #change the column names back to dates
  
  period_ints <- det_hist %>% select(matches("\\d")) %>% names()
  
  
  period_start_dates <-  det_hist %>%
    select(matches("\\d")) %>% 
    names() %>% 
    str_remove("_.*")
  
  period_end_dates <-  det_hist %>%
    select(matches("\\d")) %>% 
    names() %>% 
    str_remove(".*_")
  
  ints <- lubridate::interval(period_start_dates, period_end_dates)
  
  
  
  period_nums <- 1:length(period_start_dates) %>% as.character()
  
  #period_names <- 1:length(period_cols) %>% as.character()
  
  det_hist_named <- det_hist #%>% 
    #rename_with(.cols= matches("\\d"), function(x){period_start_dates}) 
  
  #pivot to long form
  pivot <- det_hist_named %>%
    as.data.frame() %>% 
    pivot_longer(matches("\\d"), names_to = "period", values_to = "detection") %>%
    # mutate(act_date = ymd(act_date),
    #        jday = yday(act_date)) %>% 
    mutate(act_year = as.factor(year(period))) #%>% 
    # pivot_wider(names_from = period, values_from = detection) %>% 
    # rename_with(.cols= matches("\\d"), function(x){period_nums}) %>% 
    # pivot_longer(matches("\\d"), names_to = "period", values_to = "detection") 

    # 
  #split by activity year, pivot each year frame wider
  split_fill <- pivot %>%
  split(.$act_year) %>% 
    map(., function(x){
      x %>% 
        #select(-act_year) %>%
        mutate(period = as.numeric(period)) %>% 
        mutate(period = ceiling(period/7)) %>% 
        mutate(period = as.character(period)) %>% 
        pivot_wider(names_from = period, values_from = detection) %>% View() 
        #fill_periods(nper = nper)
         
    }) 

  
  
  #for each year frame, filter only for deployments that are labelled with that year. Return NULL for years with no activity
  for(i in 1:length(split_fill)){
    year_name <- names(split_fill)[i]
    
    col_mat <- split_fill[[i]] %>% select(matches("\\d"))
    
    if(all(is.na(col_mat))){
      split_fill[[i]] <- NA
    } else{
    
    tmp <- split_fill[[i]] %>% 
      filter(str_detect(deployment_id, pattern = year_name))
    }
    
    if(nrow(tmp) == 0){
      split_fill[[i]] <- NA
    } else{
      split_fill[[i]] <- tmp
    }
  }
  
  filt <- split_fill[!is.na(split_fill)]
  
  #bind the individual annual frames together
  out <- filt %>% 
    bind_rows() 
}

## sum_or_na ##
#for grouped rows in wideform activity histories: take the sum of rows by column or return NA if all values are NA
sum_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA)  # Return NA if all values are NA
  } else {
    return(sum(x, na.rm = TRUE))  # Return the sum if at least one value is not NA
  }
}

max_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA)  # Return NA if all values are NA
  } else {
    return(max(x, na.rm = TRUE))  # Return the sum if at least one value is not NA
  }
}

################################ Libraries #################################
library(tidyverse)
library(janitor)
library(camtrapR)
library(terra)
library(sf)

###############################################################################

#1. Import activity files

###############################################################################

#LEKT
lekt_act <- read_csv("data/Camera_Data/all_species/LEKT/lekt_cam_act_2019_2023_long.csv") %>% 
  mutate(grid = "LEKT", .before = survey_id) %>% 
  relocate(viewshed, .after = pull_date)

#MAKAH
makah_act <- read_csv("data/Camera_Data/all_species/MAKAH/makah_cam_act_2021_long.csv") %>% 
  mutate(grid = "MAKAH", .before = survey_id) 

#PNPTC
pnptc_act <- read_csv("data/Camera_Data/all_species/PNPTC/pnptc_cam_act_2020_2022_long.csv") %>% 
  mutate(grid = "PNPTC", .before = survey_id)

#QUIN
quin_act <- read_csv("data/Camera_Data/all_species/QUIN/quin_cam_act_2022_long.csv") %>% 
  mutate(grid = "QUIN", .before = survey_id)

#SKOK
skok_act <- read_csv("data/Camera_Data/all_species/SKOK/skok_act_all_species_2021_2023_long.csv") %>% 
  mutate(grid = "SKOK", .before = survey_id) %>% 
  relocate(viewshed, .after = pull_date)

#ONP
onp_act <- read_csv("data/Camera_Data/all_species/ONP/onp_cam_act_2013-2016_long.csv") %>% 
  mutate(grid = "ONP", .before = survey_id) %>% 
  relocate(viewshed, .after = pull_date)


###############################################################################

#1. Combine activity files

###############################################################################

#OCP only
ocp_act <- bind_rows(lekt_act,
                     makah_act,
                     pnptc_act,
                     quin_act,
                     skok_act)

#OCP with Olympic national park
all_act <- bind_rows(ocp_act,
                     onp_act)


###############################################################################

#3. Plot combined activity matrix

###############################################################################

ocp_act_wide <- ocp_act %>% 
  arrange(set_date, survey_id) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>%
  select(-grid) %>% 
  fill_missing_dates(mode = "survey") 


ocp_act_mat_wide <- ocp_act_wide %>% 
  select(-c(survey_id, station_id:dep_year)) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()

#plot camera operation matrix
camtrapR:::camopPlot(camOp = ocp_act_mat_wide, palette = "Heat", lattice = TRUE)


all_act_wide <- all_act %>% 
  arrange(set_date, survey_id) %>% 
  pivot_wider(names_from = act_date, values_from = status) %>% 
  select(-grid) %>% 
  fill_missing_dates(mode = "survey") 
  
all_act_mat_wide <- all_act_wide %>% 
  select(-c(survey_id, station_id:dep_year)) %>%
  column_to_rownames("deployment_id") %>% 
  as.matrix()

#plot camera operation matrix
camtrapR:::camopPlot(camOp = all_act_mat_wide, palette = "Heat", lattice = TRUE)


###############################################################################

#4. Import detection files

###############################################################################

#LEKT
lekt_det <- read_csv("data/Camera_Data/all_species/LEKT/lekt_detections_all_species_2019-2023.csv", 
                     col_types = cols(timestamp = col_character(), 
                                      date = col_character(), 
                                      time = col_character())) %>% 
  mutate(grid = "LEKT", .before = deployment_id) #%>%
  #mutate(timestamp = ymd_hms(timestamp)) %>% 
  #select(-c(date, time))

#MAKAH
makah_det <- read_csv("data/Camera_Data/all_species/MAKAH/makah_detections_all_species_2021.csv",
                      col_types = cols(timestamp = col_character(), 
                                       date = col_character(), 
                                       time = col_character())) %>% 
  mutate(grid = "MAKAH", .before = deployment_id) #%>% 
  #select(-c(date, time))

#PNPTC
pnptc_det <- read_csv("data/Camera_Data/all_species/PNPTC/pnptc_detections_all_species_2020-2022.csv",
                      col_types = cols(timestamp = col_character(), 
                                       date = col_character(), 
                                       time = col_character())) %>% 
  mutate(grid = "PNPTC", .before = deployment_id) #%>% 
  #select(-c(date, time))

#QUIN
quin_det <- read_csv("data/Camera_Data/all_species/QUIN/quin_detections_all_species_2022.csv",
                     col_types = cols(timestamp = col_character(), 
                                      date = col_character(), 
                                      time = col_character())) %>% 
  mutate(grid = "QUIN", .before = deployment_id) #%>% 
  #select(-c(date, time)) %>% 
 

#SKOK
skok_det <- read_csv("data/Camera_Data/all_species/SKOK/skok_detections_all_species_2019-2023.csv",
                     col_types = cols(timestamp = col_character(), 
                                      date = col_character(), 
                                      time = col_character())) %>% 
  mutate(grid = "SKOK", .before = deployment_id) #%>% 
  # mutate(timestamp = ymd_hms(timestamp)) #%>% 
  # mutate(across(c(timestamp, date, time), as.character)) %>% 
  # mutate(timestamp = case_when(is.na(timestamp) ~ paste0(date, " ", time),
  #                              .default = timestamp)) %>% 
  # mutate(timestamp = ymd_hms(timestamp)) %>% 
  #mutate(timestamp = as.character(timestamp)) %>% 
  #select(-c(date, time))

#ONP
onp_det <- read_csv("data/Camera_Data/all_species/ONP/onp_detections_all_species_2013-2016.csv",
                    col_types = cols(timestamp = col_character())) %>% 
  mutate(grid = "ONP", .before = deployment_id) %>%
  select(-station_id) %>% 
  separate_wider_delim(timestamp, delim = " " , names = c("date", "time"), cols_remove = F, too_few = "align_start") %>% relocate(timestamp, .before = date) #%>% 
  #select(-c(date, time)) #%>% 
 # mutate(timestamp = as.character(timestamp))


###############################################################################

#4. Combine detection files

###############################################################################

#OCP with Olympic national park
all_det <- bind_rows(lekt_det,
                     makah_det,
                     pnptc_det,
                     quin_det,
                     skok_det,
                     onp_det)


ocp_det <- bind_rows(lekt_det,
                     makah_det,
                     pnptc_det,
                     quin_det,
                     skok_det)

###############################################################################

#4. Standardize species names

###############################################################################

sp_names <- all_det %>% 
  filter(is.na(species)) %>% 
  distinct(species) %>% 
  pull(species)
  
#write_csv(tibble(old_name = sp_names), "data/Camera_Data/all_species/all_op_species_dictionary.csv")
  
sp_dictionary <- read_csv("data/Camera_Data/all_species/all_op_species_dictionary.csv")

all_det_standard <- all_det %>%
  rename(old_name = species) %>% 
  left_join(sp_dictionary, by = "old_name") %>%
  rename(common_name = standard_name) %>% 
  select(grid:act_year, common_name:age, timestamp, date, time)

#make sure all matches were successful
all_det_standard %>% filter(is.na(common_name))

ocp_det_standard <- all_det_standard %>% filter(grid!="ONP")
  
#look up species with camtrapR: clunky
#test <- checkSpeciesNames(sp_names, searchtype = "common")

#export standardized detection data
# write_csv(all_det_standard, "data/Camera_Data/all_species/species_detections_standardized_ocp_onp_2013-2023.csv")
# write_csv(ocp_det_standard , "data/Camera_Data/all_species/species_detections_standardized_ocp_2019-2023.csv")

###############################################################################

#3. Plot combined detection matrix for one species

###############################################################################


dh_test <- detectionHistory(all_det_standard,
                         "Black_tailed_deer",
                         all_act_mat_wide,
                         output = "binary",
                         stationCol = "deployment_id",
                         speciesCol = "common_name",
                         recordDateTimeCol = "timestamp",
                         recordDateTimeFormat = "ymd HMS",
                         occasionLength = 1,
                         includeEffort = F,
                         day1 = "survey",
                         datesAsOccasionNames = F)

det_hist <- dh_test$detection_history

colnames(det_hist) <- colnames( all_act_mat_wide)

camtrapR:::camopPlot(det_hist, lattice = TRUE)


all_det_standard2 <- all_det_standard %>% 
  mutate(common_name = case_when(common_name == "Odocoileus_sp" ~ "Black_tailed_deer",
                                 .default = common_name)) #%>% 
  #mutate(timestamp = ymd_hms(timestamp)) %>% filter(is.na(timestamp))

all_det_standard2 %>% filter(common_name == "Odocoileus_sp")


all_det_standard3 <- all_det_standard2 %>% 
  mutate(timestamp = paste0(date, " ", time)) #%>% 
  #mutate(timestamp = ymd_hms(timestamp))
  #mutate(timestamp = str_trim(timestamp)) %>% 
  #mutate(timestamp = str_replace(timestamp, "00:00:00", "00:00:01")) #%>% 
 # mutate(timestamp = ymd_hms(format(as.POSIXct(timestamp), format = "%Y-%m-%d %T"))) 


###############################################################################

#3. Make wide-form detection histories for species of interest

###############################################################################


species_list <- c("Puma", "Black_tailed_deer", "Roosevelt_elk", "American_black_bear")


#station_ids <- unique(stations$station_id)

#Species detection plot
# all_act_wide_plot <- all_act_wide %>% 
# filter(deployment_id != "PNPTC_2022_PNPTC_DNR16")
# 
# all_det_plot <- all_det_standard3 %>% 
#   filter(deployment_id != "PNPTC_2022_PNPTC_DNR16")
# 
# test <- detectionMaps(all_act_wide_plot,
#                       all_det_plot,
#                       Xcol = "longitude",
#                       Ycol = "latitude",
#                       speciesToShow = "Black_tailed_deer",
#                       stationCol = "deployment_id",
#                       speciesCol = "common_name",
#                       smallPoints = 3)


det_hists_wide <- list()

#effort_mats_wide <- list()

for (i in 1:length(species_list)){
  
  species_i <- species_list[i]
  
  cmr<- camtrapR::detectionHistory(recordTable = all_det_standard3, 
                                   camOp = all_act_mat_wide, 
                                   species = species_i, 
                                   output = "binary",
                                   stationCol = "deployment_id",
                                   speciesCol = "common_name",
                                   recordDateTimeCol = "timestamp",
                                   recordDateTimeFormat = "ymd HMS",
                                   occasionLength = 1,
                                   day1 = "survey",
                                   includeEffort = TRUE,#effort matrix is always the same
                                   scaleEffort = FALSE,
                                   timeZone = "UTC",
                                   datesAsOccasionNames = F)
  
  
  det_hist <- cmr$detection_history %>% 
    as.data.frame() %>% 
    rownames_to_column("deployment_id") %>% 
    as_tibble()
  
  #effort_mats_wide[[i]] <- cmr$effort
  
  #bind detection history back to metadata
  
  metadata <- all_act_wide %>% 
    select(deployment_id, !matches("\\d"))
  
  det_hists_wide[[i]] <- metadata %>% 
    left_join(det_hist, by = "deployment_id") %>% 
    rename_with(.cols= matches("\\d"), function(x){as.character(colnames(all_act_mat_wide))})
  
  print(species_i)
}

names(det_hists_wide) <- species_list
#names(effort_mats_wide) <- species_list


#format and export effort matrix (only one for all species)
effort_mat_j <- cmr$effort %>% 
  as.data.frame() %>% 
  rownames_to_column("deployment_id") %>% 
  as_tibble()

effort_mat <- metadata %>% 
  left_join(effort_mat_j, by = "deployment_id") #%>% 
  #rename_with(.cols= matches("\\d"), function(x){as.character(colnames(all_act_mat_wide))})

effort_mat %>% 
  select(matches("\\d")) %>% 
  as.matrix() %>% 
  raster::raster() %>% 
  plot()

#write_csv(effort_mat, "effort_mat_all_species_daily_ocp_onp_2013-2023_weekly.csv")

det_hists_wide[[1]] %>% 
  select(matches("\\d")) %>% 
  as.matrix() %>% 
  raster::raster() %>% 
  plot()

#export detection histories for later reference
# for(i in 1:length(det_hists_wide)){
# 
#   dh_name <- names(det_hists_wide)[i] %>% str_replace(coll(" "), coll("_"))
#   dh_out <- det_hists_wide[[i]] %>% as.data.frame()
# 
#   write_csv(dh_out, paste0(dh_name, "_det_hist_wide_2013-2023_weekly.csv"))
# 
#   print(i)
# }



# det_hist_test <- cmr$detection_history
# 
# colnames(det_hist_test) <- colnames(all_act_mat_wide)
# 
# camtrapR:::camopPlot(det_hist, lattice = TRUE)

###############################################################################

#7. Stack the detection history

###############################################################################


c_names <- colnames(all_act_mat_wide)

det_hists_stacked <- map(det_hists_wide, stack_detection, c_names = c_names)

stack_effort <- stack_detection(effort_mat, c_names)


#function to agreggate stacked detection and effort matricies to periods of a specified number of days
agg_stacked_occasions <- function(det_hist, ndays, type = "detection"){
  
  if(type == "detection"){
    fun = max_or_na
  } else if(type == "effort"){
    fun = sum_or_na
  }
  
  out <- det_hist %>% 
    pivot_longer(matches("\\d"), names_to = "period", values_to = "detection") %>% 
    mutate(period = as.numeric(period)) %>% 
    mutate(interval = (period - 1) %/% ndays + 1) %>% 
    group_by(deployment_id, interval) %>% 
    summarize(across(c(survey_id:dep_year), first),
              detection = fun(detection), .groups = "drop") %>% 
    ungroup() %>% 
    pivot_wider(names_from = interval, values_from = detection)
  
  return(out)
}



effort_stacked_weekly <- agg_stacked_occasions(stack_effort, 7, type = "effort")

effort_stacked_2_weeks <- agg_stacked_occasions(stack_effort, 14, type = "effort")



det_hists_stacked_weekly <- map(det_hists_stacked, 
                                agg_stacked_occasions,
                                ndays = 7,
                                type = "detection",
                                .progress = T)

det_hists_stacked_2_weeks <- map(det_hists_stacked, 
                                agg_stacked_occasions,
                                ndays = 14,
                                type = "detection",
                                .progress = T)


# write_csv(stack_effort, "effort_mat_all_species_daily_ocp_onp_2013-2023_stacked.csv")

#export stacked detection histories 
# for(i in 1:length(det_hists_stacked)){
# 
#   dh_name <- names(det_hists_stacked)[i] %>% str_replace(coll(" "), coll("_"))
#   dh_out <- det_hists_stacked[[i]] %>% as.data.frame()
# 
#   write_csv(dh_out, paste0(dh_name, "_det_hist_stacked_2013-2023.csv"))
# 
#   print(i)
# }




#filter(!all(is.na(c_across(-deployment_id))))


#format and plot using camtrapR

date_seq <- seq(mdy("01/01/2020"), mdy("12/31/2020"), by = "day") %>% as.character()

plot_mat <- stack_effort %>% 
  select(deployment_id, matches("\\d")) %>% 
  rename_with(.cols= matches("\\d"), function(x){as.character(date_seq)}) %>% 
  column_to_rownames("deployment_id") %>% 
  as.matrix()


camtrapR:::camopPlot(plot_mat, lattice = TRUE)


###############################################################################

#7. Stack the detection history v2

###############################################################################


#  Source functions to make occupancy surveys of different lenghts
# source("/Users/tb201494/Library/CloudStorage/Box-Box/Panthera Tiger Program/Panthera_Tekai/Analysis_materials_Chin/Occupancy/example1/make_int_fun.R")
# 
# 
# #  Sampling occasion specs
# survey_period <- "week" 
# number_of_periods <- "1"  # how many days you want each survey to be
# 
# # Note date column should be date format:
# # $ date   : Date, format: "2009-05-02" "2009-08-10" "2009-05-14" ...
# 
# 
# #  Format encounter history 
# #   Observations
# obs_int <- det_hists_wide[[1]] %>% 
#   fill_missing_dates(mode = "all")
#   pivot_longer(matches("\\d"), names_to = "act_date", values_to = "detection") %>% mutate(act_date = ymd(act_date)) %>% 
#   mutate(act_year = year(act_date)) %>% 
#   dplyr::select(deployment_id, act_year, act_date, detection) %>% 
#   as_tibble()
# 
# #compute the new sampling intervals
# obs_tmp1 <- make_interval(obs_int, act_date, survey_period, number_of_periods) 
# 
# 
# #stagger the intervals properly for grouping
# obs_tmp2 <- obs_tmp1 %>% 
#   dplyr::group_by(deployment_id, act_year) %>% 
#   dplyr::mutate(survey_interval = interval - min(interval) + 1) %>%
#   dplyr::ungroup() 
# 
# #lump data by the new intervals-takes way too long with this structure
# obs_tmp3 <- obs_tmp2 %>%
#   dplyr::select(deployment_id, act_year, detection, survey_interval) %>% 
#   group_by(survey_interval, deployment_id, act_year) %>% 
#   mutate(detection = max(detection, na.rm = TRUE)) %>%
#   slice(1) %>%
#   as.data.frame() 






