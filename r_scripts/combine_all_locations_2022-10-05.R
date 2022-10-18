##### Create Master Location Data Frame from Source Files #####


#Author: Read Barbee

#Creation Date: 2022-10-05
#Updated: 2022-10-18

#Purpose: Combine all source files into single master location dataframe and check for completeness

library(tidyverse)
library(lubridate)
library(janitor)


################## Load Source Files ####################

lotek_web <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/lotek_web_final_2022-10-04.csv") %>% 
  mutate(collar_id=as.character(collar_id))

vec_web <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/vec_web_final_2022-10-04.csv") %>% 
  mutate(collar_id=as.character(collar_id))

retrieved <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/retrieved_collars_final_2022-10-05.csv") %>% 
  mutate(collar_id=as.character(collar_id))

historic <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/Source Files/Formatted/hist_downloads_final_2022-10-17.csv") %>% 
  mutate(sats_used =  str_replace_all(sats_used, "N/A", NA_character_),
         sats_used = as.numeric(sats_used)) %>% 
  mutate(deployment_id= str_replace_all(deployment_id, "Tswift_38019", "TSwift_38019")) %>% 
  mutate(deployment_id = str_replace_all(deployment_id, " ", ""))


all_locations <- bind_rows(lotek_web, vec_web, retrieved, historic)

all_locations <-all_locations %>% 
  mutate(deployment_id = str_replace_all(deployment_id, " ", "")) %>% 
  mutate(animal_id = str_replace_all(animal_id, " ", ""))

unique(all_locations$deployment_id)

deployments <- read_csv("/Users/tb201494/Library/CloudStorage/Box-Box/olympic_cougar_connectivity/data/Location Data/OCP_Cougar_Deployments_9-30-22.csv") %>% 
  clean_names() %>% 
  mutate(deployment_id = paste0(name,"_",collar_id),
         deployment_date=mdy(deployment_date),
         end_date=mdy(end_date))

setdiff(all_locations$deployment_id, deployments$deployment_id)
setdiff(deployments$deployment_id, all_locations$deployment_id)

##Only 143 unique deployments because Rue_87529 is duplicated in deployment file

deployments %>% 
  mutate(file_source = case_when(deployment_id %in% lotek_web$deployment_id ~ "lotek_web", deployment_id %in% vec_web$deployment_id ~ "vectronic_web", deployment_id %in% retrieved$deployment_id ~ "collar", deployment_id %in% historic$deployment_id ~ "historic_download", TRUE ~ "missing")) %>% 
  View()


#Write to csv

#write_csv(all_locations, "data/Location Data/Source Files/all_locations_2022-10-17.csv")


### trim all tracks based on deployment dates
#Create function to extract deployments based on collar ID and deployment dates and append columns of animal_id and deployment_id (Working)

extract_deployments <- function(data_p, animal_id_p, collar_id_p, start_date_p, end_date_p) {
  if (is.na(end_date_p) == TRUE){ #don't filter by end date if it's not included
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(date_time_gmt>= as.POSIXct(start_date_p))
  } else{
    trimmed_track <- data_p %>% 
      filter(collar_id == collar_id_p) %>% 
      filter(date_time_gmt>= as.POSIXct(start_date_p) & date_time_gmt<=as.POSIXct(end_date_p)) 
  }
  
  return(trimmed_track)
}
#create empty list to fill with data frames of trimmed individual tracks
trimmed_deployments <- list()

#loop to generate trimmed tracks for each deployment and append them to list (working)
for (i in 1:nrow(deployments)) {
  trimmed_deployments[[i]]<- extract_deployments(
    data_p = all_locations,
    animal_id_p = deployments$name[[i]],
    collar_id_p = deployments$collar_id[[i]],
    start_date_p = deployments$deployment_date[[i]],
    end_date_p = deployments$end_date[[i]])
}

all_locations_trimmed <- bind_rows(trimmed_deployments)

#trimming removed 193 locations (probably all from historic downloads....but should there be that many?)


#write_csv(all_locations_trimmed, "data/Location Data/Source Files/all_locations_trimmed_2022-10-18.csv")





