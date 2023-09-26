##### Lotek Data Formatting v2 #####


#Author: Read Barbee

#Date: 2022-09-30

#Purpose: format complete location download from Lotek for combination with vectronic web data and collar downloads. Get everything into Movebank format.

##Steps:
#1. Import all collar locations from Lotek and deployment list from local file
#1b (optional): generate summary counts of collars and deployments to compare between web and local list
#2. Filter Lotek data for collars in deployment list
#3. Filter out retrieved collars? (could do this via matching after all data are combined) Note: figure out if retrieved collars contain all data on site and more or if there is data on the site that's not on the collars

#4. Trim each track by deployment dates and assign to individual based on deployment list

#5. Combine all deployments in to single dataframe


library(tidyverse)
library(lubridate)
library(janitor)
library(collar)


################################## Step 1: Data Import ###########################

#1. import csv from complete lotek site download
lotek_complete_raw <- read_csv("data/Location Data/Raw Data/Lotek/Lotek_complete_download_raw_2022-09-27.csv") %>% clean_names()

#2. import cougar deployment list from local file and 
cougar_deployments <- read_csv("data/Location Data/OCP_Cougar_Deployments_9-30-22.csv")%>% clean_names()

#3. format dates
cougar_deployments <- cougar_deployments %>% 
  mutate(deployment_date = mdy(deployment_date),
         end_date = mdy(end_date))


################################## Step 2 & 3: Filter Lotek Data #######################

#1. Create vector of IDs for retrieved collars
retrieved_collars <- cougar_deployments %>% 
  filter(file_source == "Collar") %>% 
  select(collar_id) %>% 
  unique() %>% 
  unlist() %>% 
  as.vector()


#2. Filter Lotek data only for collars in deployment list and remove data for retrieved collars
lotek_web_filtered <- lotek_complete_raw %>% 
  filter(device_id %in% cougar_deployments$collar_id) %>% 
  filter(!(device_id %in% retrieved_collars))


#3. subset lotek collars from deployment list
cougar_deployments_lotek <- cougar_deployments %>% 
  filter(collar_brand == "Lotek")


#4. subset collar deployments not on the website
not_on_web<- setdiff(cougar_deployments_lotek$collar_id, lotek_web_filtered$device_id)


#5. remove collar deployments not on the website from deployment list for trimming
cougar_deployments_lotek_web <- cougar_deployments_lotek %>% 
  filter(!(collar_id %in% not_on_web))


##################### Step 4: Trim Tracks by Deployment Dates #######################

#Create function to extract deployments based on collar ID and deployment dates and append columns of animal_id and deployment_id (Working)

extract_deployments <- function(data, animal_id, collar_id, start_date, end_date) {
  if (is.na(end_date) == TRUE){ #don't filter by end date if it's not included
    trimmed_track <- data %>% 
      filter(device_id == collar_id) %>% 
      filter(date_time_gmt>= as.POSIXct(start_date)) %>% 
      mutate(animal_id = animal_id,
             deployment_id = paste(animal_id, "_", collar_id)) %>% 
      select(deployment_id, animal_id, everything(), -device_name)
  } else{
    trimmed_track <- data %>% 
      filter(device_id == collar_id) %>% 
      filter(date_time_gmt>= as.POSIXct(start_date) & date_time_gmt<=as.POSIXct(end_date)) %>% 
      mutate(animal_id = animal_id,
             deployment_id = paste(animal_id, "_", collar_id)) %>% 
      select(deployment_id, animal_id, everything(), -device_name)
  }
  
  return(trimmed_track)
}

#create empty list to fill with data frames of trimmed individual tracks
trimmed_deployments <- list()

#loop to generate trimmed tracks for each deployment and append them to list (working)
for (i in 1:nrow(cougar_deployments_lotek_web)) {
  trimmed_deployments[[i]]<- extract_deployments(
    lotek_web_filtered,
    cougar_deployments_lotek_web$name[[i]],
    cougar_deployments_lotek_web$collar_id[[i]],
    cougar_deployments_lotek_web$deployment_date[[i]],
    cougar_deployments_lotek_web$end_date[[i]])
}



####################Step 5: Combine Trimmed Tracks to Single DF #######################

lotek_web_final <- bind_rows(trimmed_deployments) %>% 
  mutate(deployment_id = str_replace_all(deployment_id, fixed(" "), ""))


#write_csv(lotek_web_final, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/lotek_web_trimmed_combined_2022-10-03.csv")


##################### Step 6: Convert Data Frame to Standard Format ############

#Standard fields:
#deployment_id
#animal_id
#collar_id
#date_time_gmt
#latitude
#longitude
#altitude_m
#dop
#fix_type
#temp_c
#main_v
#back_v

lotek_web_formatted <- lotek_web_final %>% 
  select(-date_time_local) %>% 
  rename(collar_id = device_id,
         fix_type = fix_status,
         altitude_m = altitude
         ) %>% 
  relocate(dop, .after = altitude_m)


#write_csv(lotek_web_formatted, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/lotek_web_final_2022-10-03.csv")



#Future Improvements:
#1. Add functionality to scrape data from lotek automatically

