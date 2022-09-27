##### Lotek Data Formatting #####


#Author: Read Barbee

#Date: 2022-09-27

#Purpose: format complete location download from Lotek for combination with vectronic web data and collar downloads. Get everything into Movebank format.


library(tidyverse)
library(lubridate)
library(janitor)


#import complete lotek csv
lotek_complete_raw <- read_csv("data/Location Data/Raw Data/Lotek/Lotek_complete_download_raw_2022-09-27.csv") %>% clean_names()

#import cougar subject metadata
cougar_info <- read_csv("data/Location Data/OCP_Cougar_Collar_Data_Files.csv") %>% 
  clean_names()


#filter just for lotek collars
cougar_info_lotek <- cougar_info %>%
  na_if("Active") %>% 
  mutate(deployment_date = mdy(deployment_date),
         end_date = mdy(end_date)
         ) %>% 
  filter(collar_brand=="Lotek")

#get list of unique collar IDs
collar_ids <- unique(cougar_info_lotek$collar_id)


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
deployments <- list()

#loop to generate trimmed tracks for each daployment and append them to list (working)
for (i in 1:nrow(cougar_info_lotek)) {
  deployments[[i]]<- extract_deployments2(
    lotek_complete_raw,
    cougar_info_lotek$name[[i]],
    cougar_info_lotek$collar_id[[i]],
    cougar_info_lotek$deployment_date[[i]],
    cougar_info_lotek$end_date[[i]])
}


## Next Steps:

#1. Filter out list elements with 0 rows (either remove deployments not on the lotek site before creating the list, or remove them afterwards by filtering the list)

#2. Inspect list elements to make sure they filtered correctly

#3. Delist and combine all list elements into a single data frame

#4. (Optional) format unified dataframe to movebank standards

#5. Repeat for Vectronics

#6. Combine into single "web" dataframe

#7. Format and combine tracks from downloaded collars


#Decision: remove individuals with downloaded collars from deployment list before extracting, or combine all and remove duplicates? Probably the former.




