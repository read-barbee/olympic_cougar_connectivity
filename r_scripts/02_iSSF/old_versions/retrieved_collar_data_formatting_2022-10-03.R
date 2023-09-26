##### Collar Download Data Formatting #####


#Author: Read Barbee

#Date: 2022-10-03

#Purpose: Filter for collared dates and combine and format files from retrieved collar downloads. Get everything into Movebank format.

library(tidyverse)
library(lubridate)
library(janitor)
library(collar)


#load all collar download files into single dataframe
collar_downloads <-
  list.files(
    path="/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Raw Data/Collar Downloads/Lotek/formatted", 
    pattern = "*.csv",
    full.names = TRUE) %>% 
  map_df(~read_csv(.)) %>% 
  clean_names() %>% 
  mutate(gmt_time = mdy_hm(gmt_time))

#Import cougar deployment list for date filtering
cougar_deployments <- read_csv("data/Location Data/OCP_Cougar_Deployments_9-30-22.csv")%>% clean_names() %>% 
  mutate(deployment_date = mdy(deployment_date),
         end_date = mdy(end_date))


#filter deployments to subset downloaded collars only
downloaded_deployments<- cougar_deployments %>% 
  filter(file_source == "Collar")


#Filter tracks by deployment times
#Create function to extract deployments based on collar ID and deployment dates and append columns of animal_id and deployment_id (Working)

extract_deployments <- function(data, animal_id, collar_id, start_date, end_date) {
  if (is.na(end_date) == TRUE){ #don't filter by end date if it's not included
    trimmed_track <- data %>% 
      filter(collar_id == collar_id) %>% 
      filter(gmt_time>= as.POSIXct(start_date)) %>% 
      mutate(animal_id = animal_id,
             deployment_id = paste(animal_id, "_", collar_id)) %>% 
      select(deployment_id, animal_id, everything())
  } else{
    trimmed_track <- data %>% 
      filter(collar_id == collar_id) %>% 
      filter(gmt_time>= as.POSIXct(start_date) & gmt_time<=as.POSIXct(end_date)) %>% 
      mutate(animal_id = animal_id,
             deployment_id = paste(animal_id, "_", collar_id)) %>% 
      select(deployment_id, animal_id, everything())
  }
  
  return(trimmed_track)
}

#create empty list to fill with data frames of trimmed individual tracks
trimmed_deployments <- list()

#loop to generate trimmed tracks for each deployment and append them to list (working)
for (i in 1:nrow(downloaded_deployments)) {
  trimmed_deployments[[i]]<- extract_deployments(
    collar_downloads,
    downloaded_deployments$name[[i]],
    downloaded_deployments$collar_id[[i]],
    downloaded_deployments$deployment_date[[i]],
    downloaded_deployments$end_date[[i]])
}


#recombine trimmed tracks into single dataframe

retrieved_collars_final <- bind_rows(trimmed_deployments) %>% 
  mutate(deployment_id = str_replace_all(deployment_id, fixed(" "), ""))


#write_csv(retrieved_collars_final, "/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/retrieved_collars_trimmed_combined_2022-10-03.csv")


#####################Convert Data Frame to Standard Format ############

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

retrieved_collars_formatted <- retrieved_collars_final %>% 
  select(-c(duration, cause_of_fix)) %>% 
  rename(date_time_gmt = gmt_time,
         altitude_m = altitude,
         temp_c = temperature,
         main_v = voltage,
         fix_type = satellites) %>% 
  relocate(c(dop, fix_type), .after= altitude_m) %>% 
  mutate(fix_type = as.character(fix_type))
  


#write_csv(retrieved_collars_formatted,"//Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/retrieved_collars_final_2022-10-03.csv")



### Next Steps: 

#Trim and concatenate data from historical lotek downloads (not cuurently on site)

#Trim and concatenate data from historical vectronics downloads (not currently on site)

#Standardize formatting across all files




#Combine all files into master location list
