##### Lotek Data Formatting #####


#Author: Read Barbee

#Date: 2022-09-27

#Purpose: format complete location download from Lotek for combination with vectronic web data and collar downloads. Get everything into Movebank format.


library(tidyverse)
library(lubridate)
library(janitor)
library(collar)


#Total OCP deployments to date: 144
#Total OCP collars deployed to date: 116

#Deployments with data from retrieved collars: 12

#Lotek
#Total OCP Lotek deployments to date: 121
#OCP lotek cougar collars to date: 102

#OCP cougar collars on Lotek site: 102


#Total OCP Lotek deployments to date: 121

#Current active OCP Lotek deployments: 33
#Bobcat collars: 10
#Patagonia collars: 7

#Total collars on Lotek site: 146


##Functions to create matrices summarizing deployment and collar counts from cougar_info

#could make these functions simpler by using group_by()

###########Optional Data Scraping #######################


lotek_user = "Elbroch"
lotek_pass = "olympicOCP123"

lotek_login(lotek_user, lotek_pass) #Not working currently

lotek_dat <- fetch_lotek_positions(device_id=selected_ids_lotek, start_date=start, end_date=end)

############################ Data Import ############################


#import complete lotek csv
lotek_complete_raw <- read_csv("data/Location Data/Raw Data/Lotek/Lotek_complete_download_raw_2022-09-27.csv") %>% clean_names()

#import cougar subject metadata
cougar_info <- read_csv("data/Location Data/OCP_Cougar_Collar_Data_Files.csv") %>% 
  clean_names()



############################ Explore Data ############################
deployment_counts <- function(deployment_list){
  total_deployments <- nrow(deployment_list)
  active_deployments <- deployment_list %>% 
    filter(end_date=="Active") %>% 
    nrow()
  complete_deployments <- deployment_list %>% 
    filter(end_date!="Active") %>% 
    nrow()
  
  total_deployments_lotek <- deployment_list %>% 
    filter(collar_brand=="Lotek") %>% 
    nrow()
  active_deployments_lotek <- deployment_list %>% 
    filter(collar_brand=="Lotek") %>% 
    filter(end_date=="Active") %>% 
    nrow()
  complete_deployments_lotek <- deployment_list %>% 
    filter(collar_brand=="Lotek") %>% 
    filter(end_date!="Active") %>% 
    nrow()
  
  total_deployments_vec <- deployment_list %>% 
    filter(collar_brand=="Vectronic") %>% 
    nrow()
  active_deployments_vec <- deployment_list %>% 
    filter(collar_brand=="Vectronic") %>% 
    filter(end_date=="Active") %>% 
    nrow()
  complete_deployments_vec <- deployment_list %>% 
    filter(collar_brand=="Vectronic") %>% 
    filter(end_date!="Active") %>% 
    nrow()
  
  lotek <- c(active_deployments_lotek, complete_deployments_lotek, total_deployments_lotek)
  vec <- c(active_deployments_vec, complete_deployments_vec, total_deployments_vec)
  total <- c(active_deployments, complete_deployments, total_deployments)
  
  summary_dat <- rbind(lotek, vec, total)
  
  colnames(summary_dat) <-  c("Active", "Complete", "Total")
  rownames(summary_dat)<- c("Lotek", "Vectronic", "All")
  
  return(summary_dat)
  
}  

collar_counts <- function(deployment_list){
  
  total_collars <- length(unique(deployment_list$collar_id))
  
  active_collars <- deployment_list %>% 
    filter(end_date=="Active") %>% 
    select(collar_id) %>% 
    unique() %>% 
    nrow()
  
  used_collars <- deployment_list %>% 
    filter(end_date!="Active") %>% 
    select(collar_id) %>% 
    unique() %>% 
    nrow()
  
  total_collars_lotek <- deployment_list %>% 
    filter(collar_brand=="Lotek") %>% 
    select(collar_id) %>% 
    unique() %>% 
    nrow()
  
  active_collars_lotek <- deployment_list %>% 
    filter(end_date=="Active") %>% 
    filter(collar_brand=="Lotek") %>% 
    select(collar_id) %>% 
    unique() %>% 
    nrow()
  
  inactive_collars_lotek <- deployment_list %>% 
    filter(end_date!="Active") %>% 
    filter(collar_brand=="Lotek") %>% 
    select(collar_id) %>% 
    unique() %>% 
    nrow()
  
  
  total_collars_vec <- deployment_list %>% 
    filter(collar_brand=="Vectronic") %>% 
    select(collar_id) %>% 
    unique() %>% 
    nrow()
  
  active_collars_vec <- deployment_list %>% 
    filter(end_date=="Active") %>% 
    filter(collar_brand=="Vectronic") %>% 
    select(collar_id) %>% 
    unique() %>% 
    nrow()
  
  inactive_collars_vec <- deployment_list %>% 
    filter(end_date!="Active") %>% 
    filter(collar_brand=="Vectronic") %>% 
    select(collar_id) %>% 
    unique() %>% 
    nrow()
  
  retrieved_collars_lotek <- deployment_list %>% 
    filter(file_source=="Collar") %>% 
    filter(collar_brand=="Lotek") %>% 
    nrow()
  
  retrieved_collars_vec <- deployment_list %>% 
    filter(file_source=="Collar") %>% 
    filter(collar_brand=="Vectronic") %>% 
    nrow()
  
  retrieved_collars_total <- deployment_list %>% 
    filter(file_source=="Collar") %>% 
    nrow()
  
  lotek <- c(active_collars_lotek, inactive_collars_lotek, total_collars_lotek)
  vec <- c(active_collars_vec, inactive_collars_vec, total_collars_vec)
  total <- c(active_collars, used_collars, total_collars)
  retrieved <- c(retrieved_collars_lotek, retrieved_collars_vec, retrieved_collars_total)
  
  summary_dat <- rbind(lotek, vec, total)
  summary_dat <- cbind(summary_dat, retrieved)
  
  colnames(summary_dat) <-  c("Active", "Complete", "Total","Retrieved")
  rownames(summary_dat)<- c("Lotek", "Vectronic", "All")
  
  return(summary_dat)
}

dep_counts <- deployment_summary()


web_inventory <- function(web_download){
  lotek_ids_all <- cougar_info %>% 
    filter(collar_brand=="Lotek") %>% 
    select(collar_id) %>% 
    unique() %>% 
    unlist() %>% 
    as.vector()
  
  lotek_ids_web <- unique(lotek_complete_raw$device_id) 
}

#determine which lotek collar IDs are not on the website
lotek_ids_all <- cougar_info %>% 
  filter(collar_brand=="Lotek") %>% 
  select(collar_id) %>% 
  unique() %>% 
  unlist() %>% 
  as.vector()

lotek_ids_web <- unique(lotek_complete_raw$device_id) 

ids_missing_from_web <- setdiff(lotek_ids_all, lotek_ids_web)

#Note: 44 collars are currently on the lotek site that are not part of the Olympic Cougar Project and/or are currently undeployed and/or are bobcat collars

#check names of individuals missing from site
indiv_missing_from_web <- cougar_info %>% 
  filter(collar_id %in% ids_missing_from_web) %>% 
  select(name, collar_id)

#all collars missing from web are old Skokomish or Quinault collars. Incorporate from downloaded file

#Create vector of collar_ids of retrieved collars with downloaded data
retrieved_collar_downloads <- cougar_info %>% 
  filter(collar_brand=="Lotek" & file_source == "Collar") %>% 
  select(collar_id) %>% 
  unique() %>% 
  unlist() %>% 
  as.vector()

#filter just for lotek collars included in the current web server and not in the list of retrieved collar downloads
cougar_info_lotek_web <- cougar_info %>%
  na_if("Active") %>% 
  mutate(deployment_date = mdy(deployment_date),
         end_date = mdy(end_date)
         ) %>% 
  filter(collar_brand=="Lotek") %>% 
  filter(!(collar_id %in% retrieved_collar_downloads)) %>% 
  filter(collar_id %in% lotek_complete_raw$device_id)
  

#get list of unique collar IDs
collar_ids <- unique(cougar_info_lotek_web$collar_id)



######################## Extract Deployments ####################################

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

#loop to generate trimmed tracks for each deployment and append them to list (working)
for (i in 1:nrow(cougar_info_lotek_web)) {
  deployments[[i]]<- extract_deployments(
    lotek_complete_raw,
    cougar_info_lotek_web$name[[i]],
    cougar_info_lotek_web$collar_id[[i]],
    cougar_info_lotek_web$deployment_date[[i]],
    cougar_info_lotek_web$end_date[[i]])
}


## Next Steps:


#2. Inspect list elements to make sure they filtered correctly

#3. Delist and combine all list elements into a single data frame

#4. (Optional) format unified dataframe to movebank standards

#5. Repeat for Vectronics

#6. Combine into single "web" dataframe

#7. Format and combine tracks from downloaded collars


#Decision: remove individuals with downloaded collars from deployment list before extracting, or combine all and remove duplicates? Probably the former.







